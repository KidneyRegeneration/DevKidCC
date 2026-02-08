"""
Test DevKidCC classification on real kidney data (Liao et al. dataset)
"""

import sys
import os
from pathlib import Path
import scanpy as sc
import pandas as pd

print("=" * 70)
print("DevKidCC Classification Test - Liao Kidney Dataset")
print("=" * 70)

# Path to the dataset
DATA_PATH = Path(r"c:\Users\sbwil\Documents\VisualStudioProjects\scRNAseq\Sean_adult-human-kidney_vs_D13p14\LiaoKidneyCombined.h5ad")

if not DATA_PATH.exists():
    print(f"\n[ERROR] Dataset not found at: {DATA_PATH}")
    print("Please update the DATA_PATH in the script.")
    sys.exit(1)

# Step 1: Load data
print("\n[1/3] Loading dataset...")
print(f"  Path: {DATA_PATH}")

try:
    adata = sc.read_h5ad(DATA_PATH)
    print(f"  [OK] Loaded successfully")
    print(f"  Shape: {adata.n_obs} cells x {adata.n_vars} genes")

    # Show some info about the dataset
    if adata.obs.columns.size > 0:
        print(f"  Metadata columns: {len(adata.obs.columns)}")
        print(f"    {', '.join(list(adata.obs.columns[:5]))}", end="")
        if len(adata.obs.columns) > 5:
            print(f", ... ({len(adata.obs.columns) - 5} more)")
        else:
            print()

except Exception as e:
    print(f"  [FAIL] Failed to load data: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 2: Initialize classifier
print("\n[2/3] Initializing DevKidCC classifier...")

try:
    from devkidcc import DevKidCCClassifier

    classifier = DevKidCCClassifier(verbose=True)
    print()

except Exception as e:
    print(f"  [FAIL] Failed to initialize classifier: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 3: Run classification
print("\n[3/3] Running DKCC classification on full dataset...")
print(f"  Processing {adata.n_obs} cells...")
print(f"  This will take several minutes...\n")

try:
    # Run classification
    adata_classified = classifier.classify(adata, copy=False)

    # Verification
    print("\n" + "=" * 70)
    print("CLASSIFICATION RESULTS")
    print("=" * 70)

    # Check for classification columns
    if 'LineageID' in adata_classified.obs.columns and 'DKCC' in adata_classified.obs.columns:
        print("\n[OK] Classification successful!")

        # Lineage distribution
        print("\n" + "-" * 70)
        print("Lineage Distribution:")
        print("-" * 70)
        lineage_counts = adata_classified.obs['LineageID'].value_counts().sort_index()
        for lineage, count in lineage_counts.items():
            pct = 100 * count / len(adata_classified)
            print(f"  {lineage:15s}: {count:6d} cells ({pct:5.1f}%)")

        # Cell type distribution
        print("\n" + "-" * 70)
        print("Cell Type Distribution (Top 15):")
        print("-" * 70)
        dkcc_counts = adata_classified.obs['DKCC'].value_counts()
        for i, (cell_type, count) in enumerate(dkcc_counts.head(15).items()):
            pct = 100 * count / len(adata_classified)
            print(f"  {i+1:2d}. {cell_type:15s}: {count:6d} cells ({pct:5.1f}%)")

        if len(dkcc_counts) > 15:
            remaining = len(dkcc_counts) - 15
            remaining_cells = dkcc_counts.iloc[15:].sum()
            pct = 100 * remaining_cells / len(adata_classified)
            print(f"  ... {remaining} more cell types: {remaining_cells:6d} cells ({pct:5.1f}%)")

        # Score columns
        score_cols = [c for c in adata_classified.obs.columns if 'scpred_' in c.lower()]
        if score_cols:
            print(f"\n[OK] Found {len(score_cols)} scpred score columns")

        # Summary statistics
        print("\n" + "-" * 70)
        print("Summary Statistics:")
        print("-" * 70)
        print(f"  Total cells:          {adata_classified.n_obs:,}")
        print(f"  Total genes:          {adata_classified.n_vars:,}")
        print(f"  Unique lineages:      {adata_classified.obs['LineageID'].nunique()}")
        print(f"  Unique cell types:    {adata_classified.obs['DKCC'].nunique()}")
        print(f"  Unassigned cells:     {(adata_classified.obs['DKCC'] == 'unassigned').sum():,}")

        # Save results
        output_path = DATA_PATH.parent / "LiaoKidneyCombined_DKCC_classified.h5ad"
        print(f"\n  Saving results to: {output_path}")
        adata_classified.write_h5ad(output_path)
        print(f"  [OK] Results saved!")

        # Also save a CSV of the classification results
        csv_path = DATA_PATH.parent / "LiaoKidneyCombined_DKCC_results.csv"
        print(f"\n  Saving classification table to: {csv_path}")

        # Select relevant columns
        result_cols = ['LineageID', 'DKCC']
        result_cols.extend([c for c in score_cols[:10]])  # Add first 10 score columns

        adata_classified.obs[result_cols].to_csv(csv_path)
        print(f"  [OK] CSV saved!")

        print("\n" + "=" * 70)
        print("TEST PASSED [OK]")
        print("=" * 70)
        print("\nClassified data saved to:")
        print(f"  - Full data: {output_path}")
        print(f"  - Results CSV: {csv_path}")
        print("\nYou can now:")
        print("  1. Load the classified data: adata = sc.read_h5ad('{0}')".format(output_path.name))
        print("  2. Visualize with scanpy (UMAP colored by DKCC)")
        print("  3. Use the classification results for downstream analysis")
        print("=" * 70)

    else:
        print("\n[FAIL] Classification columns not found!")
        print(f"Available columns: {list(adata_classified.obs.columns)}")
        sys.exit(1)

except KeyboardInterrupt:
    print("\n\n[INTERRUPTED] Classification cancelled by user")
    sys.exit(1)

except Exception as e:
    print(f"\n[FAIL] Classification failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
