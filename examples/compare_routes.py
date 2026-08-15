#!/usr/bin/env python
"""
Compare the two entry points cell by cell.

    python compare_routes.py --h5ad sample_DKCC.h5ad --labels sample_labels.csv

The container routes by format -- .h5ad is read by Python, .rds by R -- so the
obvious question is whether the same cells get the same answer either way. Build
both inputs from one source matrix, classify each through its own route, and
this reports the agreement.

They should agree closely but need not agree perfectly. Both routes end in the
same DevKidCC::DKCC(), but the Python wrapper projects onto the reference genes
before handing counts to R, which changes the library size that Seurat's
LogNormalize divides by. The wrapper ships the full-matrix library sizes across
in metadata to correct for exactly this; the residual is what that correction
does not reach, plus KNN smoothing, which is neighbourhood-dependent and so
sensitive to small upstream differences.

A high nineties agreement is the expected result. A large disagreement, or one
concentrated in a single lineage, means the two routes are no longer doing the
same thing -- which is the regression this script exists to catch.
"""

from __future__ import annotations

import argparse
import sys

import anndata as ad
import pandas as pd


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--h5ad", required=True,
                   help="Classified .h5ad from the Python route")
    p.add_argument("--labels", required=True,
                   help="CSV from export_labels.R for the R route")
    p.add_argument("--out", help="Where to write a per-column agreement CSV")
    p.add_argument("--min-agreement", type=float, default=0.0,
                   help="Exit non-zero if LineageID agreement falls below this "
                        "fraction (default: %(default)s, i.e. report only)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    py = ad.read_h5ad(args.h5ad).obs
    r = pd.read_csv(args.labels).set_index("cell")

    shared = py.index.intersection(r.index)
    print(f"cells: python {len(py)}, R {len(r)}, shared {len(shared)}")
    if len(shared) == 0:
        print("ERROR: no cell names in common -- are these the same input?",
              file=sys.stderr)
        return 1

    rows = []
    for column in ("LineageID", "DKCC"):
        a = py.loc[shared, column].astype(str)
        b = r.loc[shared, column].astype(str)
        agree = float((a.values == b.values).mean())

        try:
            from sklearn.metrics import cohen_kappa_score
            kappa = float(cohen_kappa_score(a.values, b.values))
        except ImportError:
            kappa = float("nan")

        rows.append({"column": column, "n_cells": len(shared),
                     "agreement": round(agree, 4), "kappa": round(kappa, 4)})
        print(f"{column:10s} agreement {agree:7.2%}   kappa {kappa:.4f}")

        # Where they differ matters more than how often. A handful of cells
        # moving between adjacent types is noise; a whole class relabelling is
        # not, and only the breakdown distinguishes them.
        diff = a[a.values != b.values]
        if len(diff):
            pairs = pd.Series(
                [f"{x} -> {y}" for x, y in zip(b[a.values != b.values], diff)]
            ).value_counts().head(5)
            print(f"  top disagreements (R -> Python):")
            for pair, n in pairs.items():
                print(f"    {pair:<34s} {n}")

    table = pd.DataFrame(rows)
    if args.out:
        table.to_csv(args.out, index=False)
        print(f"Wrote {args.out}")

    lineage_agreement = table.loc[table["column"] == "LineageID", "agreement"].iloc[0]
    if lineage_agreement < args.min_agreement:
        print(f"\nFAIL: LineageID agreement {lineage_agreement:.2%} is below the "
              f"{args.min_agreement:.2%} floor", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
