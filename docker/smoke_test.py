#!/usr/bin/env python
"""
Prove that both halves of this image actually run.

The container ships two entry points -- Python for h5ad, R for Seurat-native
files -- and until now CI only proved the image *built*. It built green five
times while both paths were broken end to end: the Python path died on a missing
UMAP reduction, the R path died inside SeuratDisk's h5ad converter. Neither is
visible from a successful `docker build`.

So: synthesise a small counts matrix over the packaged reference genes and push
it through every way a user can reach the classifier --

    1. the Python API,      devkidcc.classify_kidney_cells()
    2. an h5ad via /opt/dkcc, which must route it to run_dkcc.py
    3. an .rds  via /opt/dkcc, which must route it to run_dkcc.R

-- requiring classification columns back from each. No network, no data
download, no fixture in git: the gene list is already inside the wheel and the
.rds is built here by Seurat itself.

Checks 2 and 3 are as much about the routing as the classifying. h5ad used to be
read by R calling back into Python through reticulate, and the h5ad-shaped bugs
all lived in that round trip; a regression that quietly sent h5ad back to R would
otherwise look identical to a passing run.

The synthetic matrix is noise, so the labels it produces are meaningless and
deliberately not asserted on. What is asserted is that each path completes and
returns `DKCC`/`LineageID` for every cell, which is exactly the property that was
silently false.

Run it yourself against a pulled image:

    singularity exec dkcc.sif python /opt/smoke_test.py
    docker run --rm ghcr.io/kidneyregeneration/dkcc:latest python /opt/smoke_test.py
"""

from __future__ import annotations

import subprocess
import sys
import traceback
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

DISPATCH = Path("/opt/dkcc")
PY_SCRIPT = Path("/opt/run_dkcc.py")
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
        raise RuntimeError(f"FAIL: packaged reference gene list missing at {path}")

    genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(genes) < N_GENES:
        raise RuntimeError(f"FAIL: only {len(genes)} reference genes, need {N_GENES}")
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
        raise RuntimeError(f"FAIL [{label}]: missing column(s) {missing}")

    for col in REQUIRED_COLUMNS:
        n_null = obs[col].isna().sum()
        if n_null:
            raise RuntimeError(f"FAIL [{label}]: {col} is null for {n_null} cells")

    counts = obs["LineageID"].value_counts()
    print(f"  [{label}] LineageID: {counts.to_dict()}")


def run(cmd: list[str], label: str) -> subprocess.CompletedProcess:
    """
    Run a subprocess, and on any failure -- not just a non-zero exit -- put its
    own output in the log. These scripts can exit 0 having written a file with no
    labels in it, and swallowing the log in that case costs a full container
    rebuild to learn nothing.
    """
    print(f"  $ {' '.join(cmd)}")
    return subprocess.run(cmd, capture_output=True, text=True)


def require_routed_to(proc: subprocess.CompletedProcess, expected: str, label: str) -> None:
    """
    The dispatcher announces which side it picked. Assert on it: sending h5ad
    back to R would still classify (R used to convert it via reticulate) and so
    would pass every other check here while reintroducing the exact round trip
    this split removed.
    """
    if expected not in proc.stdout:
        raise RuntimeError(
            f"FAIL [{label}]: expected the dispatcher to route to the {expected} "
            "entry point; it did not say so"
        )
    print(f"  [{label}] routed to the {expected} entry point")


def check_python_api(adata: ad.AnnData) -> None:
    print("1. Python API: devkidcc.classify_kidney_cells()")
    import devkidcc

    print(f"  wrapper: {devkidcc.__file__}")
    result = devkidcc.classify_kidney_cells(adata, copy=True, verbose=False)

    # The wrapper projects onto the reference genes before handing off to R; the
    # object it returns must still be the caller's, at full width.
    if result.n_vars != adata.n_vars:
        raise RuntimeError(
            f"FAIL [python-api]: returned {result.n_vars} genes, input had "
            f"{adata.n_vars} -- the reference-gene projection leaked into the "
            "caller's object"
        )
    check_columns(result.obs, "python-api")


def check_h5ad_path(adata: ad.AnnData, workdir: Path) -> None:
    print(f"2. h5ad via {DISPATCH} (must route to Python)")
    for path in (DISPATCH, PY_SCRIPT):
        if not path.exists():
            raise RuntimeError(f"FAIL [h5ad]: {path} not in the image")

    src = workdir / "smoke_in.h5ad"
    dst = workdir / "smoke_out.h5ad"
    adata.write_h5ad(src)

    proc = run([str(DISPATCH), "--input", str(src), "--output", str(dst)], "h5ad")
    try:
        require_routed_to(proc, "Python", "h5ad")
        if proc.returncode != 0:
            raise RuntimeError(f"FAIL [h5ad]: dispatcher exited {proc.returncode}")
        if not dst.exists():
            raise RuntimeError(f"FAIL [h5ad]: no output written at {dst}")
        check_columns(ad.read_h5ad(dst).obs, "h5ad")
    except Exception:
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise


# Built by Seurat rather than converted from the h5ad: converting is precisely
# what this image no longer does, and a fixture that needed the conversion
# libraries to exist would defeat the point of removing them.
BUILD_RDS = r"""
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
counts <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
seu <- CreateSeuratObject(counts = counts, project = "smoke")
saveRDS(seu, args[2])
cat("built", args[2], "with", ncol(seu), "cells x", nrow(seu), "genes\n")
"""

# The R object's metadata comes back out as a CSV so the assertions stay in one
# place, in Python, for every path.
DUMP_OBS = r"""
args <- commandArgs(trailingOnly = TRUE)
seu <- readRDS(args[1])
write.csv(seu[[]], args[2])
"""


def check_r_path(adata: ad.AnnData, workdir: Path) -> None:
    print(f"3. .rds via {DISPATCH} (must route to R)")
    if not R_SCRIPT.exists():
        raise RuntimeError(f"FAIL [rds]: {R_SCRIPT} not in the image")

    counts_csv = workdir / "smoke_counts.csv"
    src = workdir / "smoke_in.rds"
    dst = workdir / "smoke_out.rds"
    obs_csv = workdir / "smoke_out_obs.csv"

    # Seurat wants genes as rows.
    pd.DataFrame(
        adata.X.T, index=adata.var_names, columns=adata.obs_names
    ).to_csv(counts_csv)

    build = run(["Rscript", "-e", BUILD_RDS, str(counts_csv), str(src)], "rds")
    if build.returncode != 0 or not src.exists():
        sys.stdout.write(build.stdout)
        sys.stderr.write(build.stderr)
        raise RuntimeError("FAIL [rds]: could not build the .rds fixture")
    print(f"  {build.stdout.strip().splitlines()[-1]}")

    proc = run([str(DISPATCH), "--input", str(src), "--output", str(dst)], "rds")
    try:
        require_routed_to(proc, "R", "rds")
        if proc.returncode != 0:
            raise RuntimeError(f"FAIL [rds]: dispatcher exited {proc.returncode}")
        if not dst.exists():
            raise RuntimeError(f"FAIL [rds]: no output written at {dst}")
    except Exception:
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise

    dump = run(["Rscript", "-e", DUMP_OBS, str(dst), str(obs_csv)], "rds")
    if dump.returncode != 0 or not obs_csv.exists():
        sys.stdout.write(dump.stdout)
        sys.stderr.write(dump.stderr)
        raise RuntimeError("FAIL [rds]: could not read metadata back out of the output")

    obs = pd.read_csv(obs_csv, index_col=0)
    print(f"  [rds] metadata written: {list(obs.columns)}")
    check_columns(obs, "rds")


def main() -> int:
    print(f"DevKidCC container smoke test -- {N_CELLS} cells x {N_GENES} genes\n")
    genes = reference_genes()
    adata = synthetic_adata(genes)

    # Every path runs even when an earlier one fails. Each CI round trip costs a
    # container build, so one run should report everything that is broken rather
    # than the first thing that is.
    failures: list[str] = []
    with tempfile.TemporaryDirectory(prefix="dkcc_smoke_") as tmp:
        workdir = Path(tmp)
        for label, check in (
            ("python-api", lambda: check_python_api(adata.copy())),
            ("h5ad", lambda: check_h5ad_path(adata.copy(), workdir)),
            ("rds", lambda: check_r_path(adata.copy(), workdir)),
        ):
            try:
                check()
            except Exception:
                traceback.print_exc()
                failures.append(label)
            print()

    if failures:
        print(f"FAIL: {', '.join(failures)} path(s) did not classify the input.")
        return 1

    print("PASS: the Python API, h5ad routing and .rds routing all classified the input.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
