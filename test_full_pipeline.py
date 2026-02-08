"""
Full end-to-end test of DevKidCC-python classification pipeline
This tests the Python wrapper calling the R DevKidCC package
"""

import sys
import tempfile
import os
from pathlib import Path

print("=" * 60)
print("DevKidCC-python Full Pipeline Test")
print("=" * 60)

# Step 1: Test imports
print("\n[1/5] Testing imports...")
try:
    from devkidcc import DevKidCCClassifier, classify_kidney_cells
    import scanpy as sc
    import anndata as ad
    import pandas as pd
    import numpy as np
    print("    [OK] All imports successful")
except Exception as e:
    print(f"    [FAIL] Import failed: {e}")
    sys.exit(1)

# Step 2: Initialize classifier (this will check/install R packages)
print("\n[2/5] Initializing DevKidCC classifier...")
print("    (This will check R package dependencies)\n")
try:
    classifier = DevKidCCClassifier(verbose=True, install_deps=True)
    print("\n    [OK] Classifier initialized")
except Exception as e:
    print(f"\n    [FAIL] Initialization failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 3: Load test data from R package
print("\n[3/5] Loading test data from R DevKidCC package...")
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    # Load the organoid dataset from R
    ro.r('''
        library(DevKidCC)
        data("organoid", package = "DevKidCC")
        test_data <- organoid
    ''')

    # Get dimensions
    n_cells = ro.r('ncol(test_data)')[0]
    n_genes = ro.r('nrow(test_data)')[0]

    print(f"    [OK] Loaded test dataset: {int(n_cells)} cells × {int(n_genes)} genes")

    # Create a temporary h5ad file
    temp_dir = tempfile.mkdtemp()
    h5ad_path = os.path.join(temp_dir, "test_data.h5ad")

    # Convert to h5seurat then h5ad
    ro.r(f'''
        library(Seurat)
        library(SeuratDisk)

        # Ensure we have a proper Seurat object
        if (!inherits(test_data, "Seurat")) {{
            test_data <- as.Seurat(test_data)
        }}

        # Update to latest Seurat version if needed
        test_data <- UpdateSeuratObject(test_data)

        # Save as h5seurat
        temp_file <- tempfile(fileext = ".h5seurat")
        SaveH5Seurat(test_data, filename = temp_file, overwrite = TRUE)

        # Convert to h5ad
        Convert(temp_file, dest = "h5ad", overwrite = TRUE)

        # Get the h5ad path
        h5ad_file <- sub("\\\\.h5seurat$", ".h5ad", temp_file)
    ''')

    h5ad_from_r = ro.r('h5ad_file')[0]

    # Read into Python
    adata = sc.read_h5ad(h5ad_from_r)
    print(f"    [OK] Converted to AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

except Exception as e:
    print(f"    [FAIL] Failed to load test data: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 4: Run classification
print("\n[4/5] Running DevKidCC classification...")
print("    (This may take several minutes)\n")
try:
    # Take a small subset for faster testing
    adata_subset = adata[:500].copy()  # Use first 500 cells
    print(f"    Using subset of {adata_subset.n_obs} cells for testing\n")

    # Run classification
    adata_classified = classifier.classify(adata_subset, copy=True)

    print("\n    [OK] Classification completed successfully")

except Exception as e:
    print(f"\n    [FAIL] Classification failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 5: Verify results
print("\n[5/5] Verifying results...")
try:
    # Check for expected columns
    expected_cols = ['LineageID', 'DKCC']
    missing_cols = [col for col in expected_cols if col not in adata_classified.obs.columns]

    if missing_cols:
        print(f"    [FAIL] Missing expected columns: {missing_cols}")
        sys.exit(1)

    print("    [OK] Found expected columns: LineageID, DKCC")

    # Show results
    print("\n    Classification Results:")
    print("    " + "-" * 50)

    if 'LineageID' in adata_classified.obs.columns:
        n_lineages = adata_classified.obs['LineageID'].nunique()
        print(f"    Lineages identified: {n_lineages}")
        print("\n    Lineage distribution:")
        lineage_counts = adata_classified.obs['LineageID'].value_counts()
        for lineage, count in lineage_counts.items():
            print(f"      {lineage}: {count}")

    if 'DKCC' in adata_classified.obs.columns:
        n_types = adata_classified.obs['DKCC'].nunique()
        print(f"\n    Cell types identified: {n_types}")
        print("\n    Top 5 cell types:")
        type_counts = adata_classified.obs['DKCC'].value_counts().head(5)
        for cell_type, count in type_counts.items():
            print(f"      {cell_type}: {count}")

    # Check for score columns
    score_cols = [c for c in adata_classified.obs.columns if 'scpred_' in c.lower()]
    if score_cols:
        print(f"\n    [OK] Found {len(score_cols)} scpred score columns")

    print("\n" + "=" * 60)
    print("TEST PASSED [OK]")
    print("=" * 60)
    print("\nThe DevKidCC-python wrapper is working correctly!")
    print("You can now use it to classify your own kidney scRNA-seq data.")
    print("\nExample usage:")
    print("    from devkidcc import classify_kidney_cells")
    print("    import scanpy as sc")
    print("    ")
    print("    adata = sc.read_h5ad('your_data.h5ad')")
    print("    adata = classify_kidney_cells(adata)")
    print("    print(adata.obs['DKCC'].value_counts())")
    print("=" * 60)

except Exception as e:
    print(f"    [FAIL] Verification failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
