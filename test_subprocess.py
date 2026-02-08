"""
Test the subprocess-based DevKidCC classifier
"""

import sys
import numpy as np
import scanpy as sc
import anndata as ad

print("=" * 60)
print("DevKidCC Subprocess-Based Classifier Test")
print("=" * 60)

# Step 1: Test imports
print("\n[1/3] Testing imports...")
try:
    from devkidcc import DevKidCCClassifier, classify_kidney_cells
    print("    [OK] Imports successful")
except Exception as e:
    print(f"    [FAIL] Import failed: {e}")
    sys.exit(1)

# Step 2: Initialize classifier
print("\n[2/3] Initializing classifier...")
try:
    classifier = DevKidCCClassifier(verbose=True)
    print()  # Add spacing after verbose output
except Exception as e:
    print(f"    [FAIL] Initialization failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Step 3: Create test data and classify
print("\n[3/3] Testing classification with synthetic data...")
print("    NOTE: This uses synthetic data, so results won't be")
print("    biologically meaningful. Use real kidney data for actual analysis.\n")

try:
    # Create minimal synthetic dataset
    n_cells = 100
    n_genes = 200

    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))

    gene_names = [f"Gene{i}" for i in range(n_genes)]
    cell_barcodes = [f"Cell{i}" for i in range(n_cells)]

    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = cell_barcodes

    print(f"    Created synthetic dataset: {adata.n_obs} cells x {adata.n_vars} genes\n")

    # Classify
    adata_classified = classifier.classify(adata, copy=True)

    # Verify
    print("\n" + "=" * 60)
    print("TEST RESULTS")
    print("=" * 60)

    expected_cols = ['LineageID', 'DKCC']
    missing_cols = [col for col in expected_cols if col not in adata_classified.obs.columns]

    if missing_cols:
        print(f"[FAIL] Missing expected columns: {missing_cols}")
        sys.exit(1)

    print("[OK] Found expected classification columns")

    print("\nClassification Results:")
    if 'LineageID' in adata_classified.obs.columns:
        print("\n  Lineages:")
        for lineage, count in adata_classified.obs['LineageID'].value_counts().items():
            print(f"    {lineage}: {count}")

    if 'DKCC' in adata_classified.obs.columns:
        print("\n  Cell Types:")
        for cell_type, count in adata_classified.obs['DKCC'].value_counts().items():
            print(f"    {cell_type}: {count}")

    print("\n" + "=" * 60)
    print("TEST PASSED [OK]")
    print("=" * 60)
    print("\nThe subprocess-based DevKidCC wrapper is working!")
    print("\nFor real analysis, use actual kidney scRNA-seq data:")
    print("  adata = sc.read_h5ad('your_kidney_data.h5ad')")
    print("  adata = classify_kidney_cells(adata)")
    print("=" * 60)

except Exception as e:
    print(f"\n[FAIL] Classification failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
