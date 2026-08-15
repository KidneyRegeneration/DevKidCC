#!/usr/bin/env python
"""
Classify an h5ad file with DevKidCC.

h5ad is AnnData's format, and this image already contains anndata and the
`devkidcc` wrapper, so reading it in Python is a direct file read. The R entry
point used to do this instead by calling back into Python through reticulate --
R -> reticulate -> anndata -> Seurat, classify, then Seurat -> reticulate ->
anndata -> disk. Every h5ad-specific bug we have hit lived in that round trip:
reticulate provisioning its own interpreter, the resulting environment having no
anndata, and scCustomize's writer dropping the classification columns from obs.

None of those failures are possible from here. h5ad now belongs to Python and
R-native formats (.rds, .RData, .h5seurat, 10x .h5) belong to run_dkcc.R;
/opt/dkcc routes by extension. It also means an h5ad gets the same treatment
however it arrives -- through this CLI or through
`devkidcc.classify_kidney_cells()` directly -- rather than two entry points
running subtly different pipelines over the same file.

Usage:
    python /opt/run_dkcc.py --input sample.h5ad --output sample_DKCC.h5ad
"""

from __future__ import annotations

import argparse
import sys

import anndata as ad

import devkidcc

# Columns DKCC() is contracted to add. Downstream code reads these, and an h5ad
# written without them is the failure that is hardest to notice.
REQUIRED_COLUMNS = ("DKCC", "LineageID")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run DevKidCC classification on an h5ad file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input .h5ad file")
    parser.add_argument("-o", "--output", required=True, help="Output .h5ad file")
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.7,
        help="scPred assignment probability threshold (default: %(default)s)",
    )
    parser.add_argument(
        "-m", "--max-iter", type=int, default=1,
        help="DKCC max.iter (default: %(default)s)",
    )
    parser.add_argument(
        "-k", "--knn-iter", type=int, default=None,
        help="KNN smoothing iterations; 0 disables smoothing (default: DKCC's own)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    print(f"Reading input: {args.input}")
    adata = ad.read_h5ad(args.input)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    kwargs = {"threshold": args.threshold, "max_iter": args.max_iter, "copy": True}
    if args.knn_iter is not None:
        kwargs["knn_iter"] = args.knn_iter
        kwargs["knn_smooth"] = args.knn_iter > 0

    result = devkidcc.classify_kidney_cells(adata, **kwargs)

    # Same check the R side makes: an output file that converted cleanly but
    # carries no labels is worse than a crash, because nothing reports it.
    missing = [c for c in REQUIRED_COLUMNS if c not in result.obs.columns]
    if missing:
        print(f"ERROR: classification produced no {'/'.join(missing)} column",
              file=sys.stderr)
        return 1

    print(f"Saving output: {args.output}")
    result.write_h5ad(args.output)
    print(f"  LineageID: {result.obs['LineageID'].value_counts().to_dict()}")
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
