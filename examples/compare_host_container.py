#!/usr/bin/env python
"""Compare two runs of the same data made in different environments.

    compare_host_container.py <dir-a> <dir-b> [dataset ...]

`compare_routes.py` answers "do the R and Python routes agree with each other?".
This answers the other half: "does the container give the same answer as the
machine the package was developed on?" -- which is what an external user pulling
the image is really asking.

Both directories must hold `<dataset>_py.h5ad` written by `classify_h5ad.py`
from the same input, so the labels line up by barcode and a join is meaningful.

The container pins its own R against whatever the host has, so near-identical is
the pass condition and bit-identical would be luck. (In the v0.5.1 validation run
it came back identical anyway, across R 4.4.3 vs 4.5.3.)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from sklearn.metrics import cohen_kappa_score

DEFAULT_DATASETS = ("organoid_howden", "fetal_menon")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dir_a", type=Path, help="results directory from one environment")
    p.add_argument("dir_b", type=Path, help="results directory from the other")
    p.add_argument("datasets", nargs="*", default=list(DEFAULT_DATASETS))
    p.add_argument("--label-a", default="A")
    p.add_argument("--label-b", default="B")
    p.add_argument("--out", type=Path, help="write the summary table here as CSV")
    args = p.parse_args()

    rows = []
    for dataset in args.datasets or DEFAULT_DATASETS:
        fa = args.dir_a / f"{dataset}_py.h5ad"
        fb = args.dir_b / f"{dataset}_py.h5ad"
        if not (fa.exists() and fb.exists()):
            print(f"skipping {dataset}: need both {fa} and {fb}")
            continue

        a_ad, b_ad = ad.read_h5ad(fa), ad.read_h5ad(fb)
        shared = a_ad.obs_names.intersection(b_ad.obs_names)
        print(f"\n=== {dataset} ===")
        print(f"{args.label_a} {a_ad.n_obs} cells, {args.label_b} {b_ad.n_obs} cells, "
              f"shared {len(shared)}")
        if len(shared) == 0:
            print("  no shared barcodes -- these are not the same cells")
            continue

        for col in ("LineageID", "DKCC"):
            a = a_ad.obs.loc[shared, col].astype(str).values
            b = b_ad.obs.loc[shared, col].astype(str).values
            agree = float((a == b).mean())
            kappa = float(cohen_kappa_score(a, b))
            print(f"{col:10s} agreement {agree * 100:6.2f}%   kappa {kappa:.4f}")

            diff = pd.DataFrame({args.label_a: a, args.label_b: b})
            diff = diff[diff[args.label_a] != diff[args.label_b]]
            counts = diff.groupby([args.label_a, args.label_b]).size()
            for (x, y), n in counts.sort_values(ascending=False).head(5).items():
                print(f"    {x} -> {y}: {n}")

            rows.append(dict(dataset=dataset, column=col, n_cells=len(shared),
                             agreement=round(agree, 4), kappa=round(kappa, 4)))

    if args.out and rows:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
