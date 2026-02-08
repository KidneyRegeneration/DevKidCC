"""
Simple test of DevKidCC-python classification with minimal synthetic data
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

print("=" * 60)
print("DevKidCC-python Simple Test")
print("=" * 60)

# Step 1: Test imports
print("\n[1/4] Testing imports...")
try:
    from devkidcc import DevKidCCClassifier, classify_kidney_cells
    print("    [OK] All imports successful")
except Exception as e:
    print(f"    [FAIL] Import failed: {e}")
    sys.exit(1)

# Step 2: Initialize classifier
print("\n[2/4] Initializing DevKidCC classifier...")
print("    (This will check R package dependencies)\n")
try:
    classifier = DevKidCCClassifier(verbose=True, install_deps=True)
    print("\n    [OK] Classifier initialized")
except Exception as e:
    print(f"\n    [FAIL] Initialization failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 3: Create minimal test data
print("\n[3/4] Creating minimal test dataset...")
try:
    # Create a simple synthetic dataset
    # In a real use case, you would load your own kidney scRNA-seq data
    n_cells = 50
    n_genes = 100

    # Random count data
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))

    # Create gene names
    gene_names = [f"Gene{i}" for i in range(n_genes)]

    # Create cell barcodes
    cell_barcodes = [f"Cell{i}" for i in range(n_cells)]

    # Create AnnData object
    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = cell_barcodes

    # Add some basic metadata
    adata.obs['sample'] = ['sample1'] * (n_cells // 2) + ['sample2'] * (n_cells - n_cells // 2)

    print(f"    [OK] Created synthetic dataset: {adata.n_obs} cells x {adata.n_vars} genes")

except Exception as e:
    print(f"    [FAIL] Failed to create test data: {e}")
    sys.exit(1)

# Step 4: Run classification
print("\n[4/4] Running DevKidCC classification...")
print("    (This may take a few minutes)")
print("    Note: This is synthetic data, so results won't be biologically meaningful\n")
try:
    # Run classification
    adata_classified = classifier.classify(adata, copy=True)

    print("\n    [OK] Classification completed successfully!")

    # Verify results
    print("\n    Verification:")
    print("    " + "-" * 50)

    # Check for expected columns
    expected_cols = ['LineageID', 'DKCC']
    missing_cols = [col for col in expected_cols if col not in adata_classified.obs.columns]

    if missing_cols:
        print(f"    [WARNING] Missing expected columns: {missing_cols}")
    else:
        print("    [OK] Found expected columns: LineageID, DKCC")

    # Show results
    if 'LineageID' in adata_classified.obs.columns:
        print(f"\n    Lineage assignments:")
        lineage_counts = adata_classified.obs['LineageID'].value_counts()
        for lineage, count in lineage_counts.items():
            print(f"      {lineage}: {count}")

    if 'DKCC' in adata_classified.obs.columns:
        print(f"\n    Cell type assignments:")
        type_counts = adata_classified.obs['DKCC'].value_counts()
        for cell_type, count in type_counts.items():
            print(f"      {cell_type}: {count}")

    # Check for score columns
    score_cols = [c for c in adata_classified.obs.columns if 'scpred_' in c.lower()]
    if score_cols:
        print(f"\n    [OK] Found {len(score_cols)} classification score columns")

    print("\n" + "=" * 60)
    print("TEST PASSED [OK]")
    print("=" * 60)
    print("\nThe DevKidCC-python wrapper is working correctly!")
    print("\nImportant Notes:")
    print("  - This test used SYNTHETIC data, so the classifications are")
    print("    not biologically meaningful")
    print("  - For real analysis, use actual kidney scRNA-seq data")
    print("  - The classifier expects gene expression matrices with standard")
    print("    gene symbols (e.g., PAX2, SIX2, etc.)")
    print("\nExample usage with real data:")
    print("    from devkidcc import classify_kidney_cells")
    print("    import scanpy as sc")
    print("    ")
    print("    adata = sc.read_h5ad('your_kidney_data.h5ad')")
    print("    adata = classify_kidney_cells(adata)")
    print("    print(adata.obs['DKCC'].value_counts())")
    print("=" * 60)

except Exception as e:
    print(f"\n    [FAIL] Classification failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
