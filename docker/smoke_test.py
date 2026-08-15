#!/usr/bin/env python
"""
Prove that both halves of this image actually run.

The container ships two entry points -- an R one (`/opt/run_dkcc.R`) and a Python
one (`import devkidcc`) -- and until now CI only proved the image *built*. It
built green five times while both paths were broken end to end: the Python path
died on a missing UMAP reduction, the R path died inside SeuratDisk's h5ad
converter. Neither is visible from a successful `docker build`.

So: synthesise a small counts matrix over the packaged reference genes, push it
through both paths, and require classification columns to come back. No network,
no data download, no fixture in git -- the gene list is already inside the
wheel.

The synthetic matrix is noise, so the labels it produces are meaningless and
deliberately not asserted on. What is asserted is that the pipeline completes and
returns `DKCC`/`LineageID` for every cell, which is exactly the property that was
silently false.

Run it yourself against a pulled image:

    singularity exec dkcc.sif python /opt/smoke_test.py
    docker run --rm ghcr.io/kidneyregeneration/dkcc:latest python /opt/smoke_test.py
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Enough cells for DKCC()'s internal neighbour graph and UMAP to be well posed,
# and small enough that the whole check runs in a couple of minutes on a CI box.
N_CELLS = 200

# A subset of the reference genes rather than all ~10k: scPred only needs enough
# overlap to align, and the runtime is dominated by the gene count.
N_GENES = 4000

# Columns DKCC() is contracted to add. These are what downstream code reads.
REQUIRED_COLUMNS = ("DKCC", "LineageID")

R_SCRIPT = Path("/opt/run_dkcc.R")


def reference_genes() -> list[str]:
    """
    The gene list shipped inside the installed wrapper. Sourcing it from the
    package rather than a checked-in copy means the smoke test cannot drift away
    from the models the image actually contains.
    """
    import devkidcc

    path = Path(devkidcc.__file__).parent / "data" / "reference_genes.txt"
    if not path.exists():
        raise SystemExit(f"FAIL: packaged reference gene list missing at {path}")

    genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(genes) < N_GENES:
        raise SystemExit(f"FAIL: only {len(genes)} reference genes, need {N_GENES}")
    return genes[:N_GENES]


def synthetic_adata(genes: list[str]) -> ad.AnnData:
    """
    A negative-binomial-ish counts matrix: Poisson draws over a lognormal mean,
    which gives the mean-variance relationship and sparsity of real UMI data
    without pretending to carry biology.
    """
    rng = np.random.default_rng(0)
    means = rng.lognormal(mean=-1.0, sigma=1.5, size=(1, len(genes)))
    counts = rng.poisson(np.repeat(means, N_CELLS, axis=0)).astype(np.float32)

    # Guarantee no empty cell: a zero library size divides by zero in
    # LogNormalize and the failure would look like a code bug, not a data one.
    empty = counts.sum(axis=1) == 0
    if empty.any():
        counts[empty, 0] = 1.0

    obs = pd.DataFrame(index=[f"cell_{i:04d}" for i in range(N_CELLS)])
    obs["orig.ident"] = "smoke"
    var = pd.DataFrame(index=pd.Index(genes, name=None))
    return ad.AnnData(X=counts, obs=obs, var=var)


def check_columns(obs: pd.DataFrame, label: str) -> None:
    missing = [c for c in REQUIRED_COLUMNS if c not in obs.columns]
    if missing:
        raise SystemExit(f"FAIL [{label}]: missing column(s) {missing}")

    for col in REQUIRED_COLUMNS:
        n_null = obs[col].isna().sum()
        if n_null:
            raise SystemExit(f"FAIL [{label}]: {col} is null for {n_null} cells")

    counts = obs["LineageID"].value_counts()
    print(f"  [{label}] LineageID: {counts.to_dict()}")


def check_python_path(adata: ad.AnnData) -> None:
    print("Python path: devkidcc.classify_kidney_cells()")
    import devkidcc

    print(f"  wrapper: {devkidcc.__file__}")
    result = devkidcc.classify_kidney_cells(adata, copy=True, verbose=False)

    # The wrapper projects onto the reference genes before handing off to R; the
    # object it returns must still be the caller's, at full width.
    if result.n_vars != adata.n_vars:
        raise SystemExit(
            f"FAIL [python]: returned {result.n_vars} genes, input had {adata.n_vars} "
            "-- the reference-gene projection leaked into the caller's object"
        )
    check_columns(result.obs, "python")


def check_r_path(adata: ad.AnnData, workdir: Path) -> None:
    print(f"R path: Rscript {R_SCRIPT}")
    if not R_SCRIPT.exists():
        raise SystemExit(f"FAIL [r]: {R_SCRIPT} not in the image")

    src = workdir / "smoke_in.h5ad"
    dst = workdir / "smoke_out.h5ad"
    adata.write_h5ad(src)

    proc = subprocess.run(
        ["Rscript", str(R_SCRIPT), "--input", str(src), "--output", str(dst)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not dst.exists():
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"FAIL [r]: run_dkcc.R exited {proc.returncode}")

    check_columns(ad.read_h5ad(dst).obs, "r")


def main() -> int:
    print(f"DevKidCC container smoke test -- {N_CELLS} cells x {N_GENES} genes\n")
    genes = reference_genes()
    adata = synthetic_adata(genes)

    with tempfile.TemporaryDirectory(prefix="dkcc_smoke_") as tmp:
        workdir = Path(tmp)
        check_python_path(adata.copy())
        print()
        check_r_path(adata.copy(), workdir)

    print("\nPASS: both the R and Python entry points classified the input.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
