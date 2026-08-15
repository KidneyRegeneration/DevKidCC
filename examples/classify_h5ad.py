#!/usr/bin/env python
"""
Classify an AnnData object with DevKidCC, from Python.

This is the interactive counterpart to `/opt/run_dkcc.py`: the same call, written
out step by step so it can be read, adapted, or lifted into a methods section.
Run it as a script or paste the body into a notebook.

    python classify_h5ad.py --input organoid.h5ad --output organoid_DKCC.h5ad

What it demonstrates
--------------------
1. DevKidCC expects **raw counts**. The models were fit against counts and
   normalisation happens inside; handing it log-normalised values gives wrong
   answers silently rather than an error, so this script checks and refuses.
2. `classify_kidney_cells()` returns your object at its original gene width.
   The wrapper projects onto the ~10k reference genes internally to keep the
   handoff to R small, but that projection must not reach the caller -- an
   earlier version rebound the object to the projection and quietly returned a
   third of the genes.
3. The labels arrive as two `.obs` columns:
     LineageID -- Nephron / Stroma / UrEp / NPC / NPC-like / Endo / unassigned
     DKCC      -- the finer type within that lineage
   `unassigned` is a real answer, not a failure: scPred rejects a cell whose
   probability does not clear the threshold, and organoid data legitimately
   contains cells outside the reference.
"""

from __future__ import annotations

import argparse
import sys

import anndata as ad
import numpy as np
import pandas as pd

import devkidcc


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("-i", "--input", required=True, help="Input .h5ad (raw counts)")
    p.add_argument("-o", "--output", help="Where to write the classified .h5ad")
    p.add_argument("-s", "--summary", help="Where to write a label-count CSV")
    p.add_argument(
        "-t", "--threshold", type=float, default=0.7,
        help="scPred probability below which a cell is left unassigned "
             "(default: %(default)s)",
    )
    p.add_argument(
        "-k", "--knn-iter", type=int, default=None,
        help="Rounds of KNN smoothing over unassigned cells; 0 disables it "
             "(default: the package default)",
    )
    return p.parse_args(argv)


def looks_log_normalised(adata: ad.AnnData) -> bool:
    """
    Log-normalised data tops out around 8-10 and carries fractional values;
    counts are non-negative integers reaching the hundreds or thousands. The
    distinction matters more than it looks: Seurat's LogNormalize divides by a
    library size, so pre-normalised input is scored against the wrong scale.
    """
    x = adata.X
    sample = x[: min(200, x.shape[0])]
    sample = sample.toarray() if hasattr(sample, "toarray") else np.asarray(sample)
    if sample.size == 0:
        return False
    return bool(sample.max() < 50 and not np.allclose(sample, np.round(sample)))


def describe_input(adata: ad.AnnData) -> None:
    x = adata.X
    sample = x[: min(200, x.shape[0])]
    sample = sample.toarray() if hasattr(sample, "toarray") else np.asarray(sample)
    print(f"Input : {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"        X max {sample.max():.1f}, "
          f"{'sparse' if hasattr(x, 'toarray') else 'dense'}, dtype {x.dtype}")


def summarise(adata: ad.AnnData) -> pd.DataFrame:
    """One tidy table of both label columns -- the thing worth reporting."""
    rows = []
    for column in ("LineageID", "DKCC"):
        counts = adata.obs[column].value_counts()
        for label, n in counts.items():
            rows.append({
                "column": column,
                "label": label,
                "n_cells": int(n),
                "percent": round(100 * n / adata.n_obs, 2),
            })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    print(f"DevKidCC {getattr(devkidcc, '__version__', 'unknown')} "
          f"via {devkidcc.__file__}")
    adata = ad.read_h5ad(args.input)
    describe_input(adata)

    if looks_log_normalised(adata):
        print("\nERROR: this looks log-normalised, not raw counts.", file=sys.stderr)
        print("       DevKidCC normalises internally; pass the raw matrix "
              "(e.g. adata.raw or a *_qc.h5ad).", file=sys.stderr)
        return 1

    n_vars_in = adata.n_vars

    kwargs = {"threshold": args.threshold, "copy": True}
    if args.knn_iter is not None:
        kwargs["knn_iter"] = args.knn_iter
        kwargs["knn_smooth"] = args.knn_iter > 0

    print("\nClassifying (this shells out to R; a few minutes for large inputs)...")
    result = devkidcc.classify_kidney_cells(adata, **kwargs)

    # The projection is an implementation detail of the handoff. If it ever
    # reaches the caller, genes have been silently deleted from their data.
    assert result.n_vars == n_vars_in, (
        f"gene width changed: {n_vars_in} -> {result.n_vars}"
    )
    missing = [c for c in ("DKCC", "LineageID") if c not in result.obs.columns]
    if missing:
        print(f"ERROR: classification returned no {missing}", file=sys.stderr)
        return 1

    print(f"\nReturned {result.n_obs} cells x {result.n_vars} genes "
          "(gene width unchanged)\n")

    table = summarise(result)
    for column in ("LineageID", "DKCC"):
        print(f"{column}:")
        block = table[table["column"] == column]
        for _, r in block.iterrows():
            print(f"  {r['label']:<16s} {r['n_cells']:>6d}  {r['percent']:>6.2f}%")
        print()

    assigned = (result.obs["LineageID"] != "unassigned").sum()
    print(f"Assigned: {assigned}/{result.n_obs} "
          f"({100 * assigned / result.n_obs:.1f}%)")

    if args.output:
        result.write_h5ad(args.output)
        print(f"Wrote {args.output}")
    if args.summary:
        table.to_csv(args.summary, index=False)
        print(f"Wrote {args.summary}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
