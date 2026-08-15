#!/usr/bin/env python
"""Build the small test inputs the other examples run on.

Two 600-cell subsets, one organoid and one fetal kidney, written in both
formats so the same cells can go down both routes:

    <name>.h5ad         AnnData, for the Python route
    <name>.mtx          MatrixMarket + <name>.genes.txt + <name>.cells.txt,
                        which build_rds_from_mtx.R turns into a Seurat .rds

600 cells is deliberate. It classifies in about a minute, fits in a couple
of GB, and is still large enough that the two routes disagreeing would show
up -- the point of the example set is the comparison, not the biology.

The subsets are drawn with a fixed seed so re-running this reproduces the
same cells, and therefore the same labels, on any machine.

Both inputs must be raw counts: DKCC's scPred models expect to normalise
the data themselves, and log-normalised input silently produces nonsense
rather than an error. The *_qc.h5ad files carry raw counts; the
*_processed.h5ad ones next to them do not.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.io
import scipy.sparse as sp

SEED = 42
N_CELLS = 600

SOURCES = {
    "organoid_howden": "/data/homeserver/data/datasets/Howden_2019_Organoids/processed/Howden_2019_Organoids_qc.h5ad",
    "fetal_menon": "/data/homeserver/data/datasets/Menon_2018_HFK/processed/Menon_2018_HFK_qc.h5ad",
}


def subset(src: Path, n: int, seed: int) -> ad.AnnData:
    a = ad.read_h5ad(src)
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(a.n_obs, n, replace=False))
    out = a[idx].copy()
    x = out.X
    if x.max() < 30:
        raise SystemExit(
            f"{src} looks log-normalised (max {x.max():.2f}); DKCC needs raw counts"
        )
    return out


def write(a: ad.AnnData, stem: Path) -> None:
    a.write_h5ad(stem.with_suffix(".h5ad"))
    # Seurat wants genes x cells, so the matrix goes out transposed.
    m = sp.csr_matrix(a.X).T.tocoo()
    scipy.io.mmwrite(str(stem) + ".mtx", m)
    stem.with_suffix(".genes.txt").write_text("\n".join(a.var_names) + "\n")
    stem.with_suffix(".cells.txt").write_text("\n".join(a.obs_names) + "\n")
    print(f"{stem.name}: {a.n_obs} cells x {a.n_vars} genes")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("."), help="output directory")
    p.add_argument("--n-cells", type=int, default=N_CELLS)
    p.add_argument("--seed", type=int, default=SEED)
    args = p.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    for name, src in SOURCES.items():
        src = Path(src)
        if not src.exists():
            print(f"skipping {name}: {src} not present on this machine")
            continue
        write(subset(src, args.n_cells, args.seed), args.out / name)

    print("\nNow build the .rds companions:")
    for name in SOURCES:
        s = args.out / name
        print(f"  Rscript build_rds_from_mtx.R {s}.mtx {s}.genes.txt {s}.cells.txt {s}.rds")


if __name__ == "__main__":
    main()
