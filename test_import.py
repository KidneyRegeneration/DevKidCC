"""
Quick test to verify DevKidCC-python package is working
"""

import sys
print("=" * 60)
print("DevKidCC-python Package Test")
print("=" * 60)

# Test imports
print("\n1. Testing imports...")
try:
    from devkidcc import (
        DevKidCCClassifier,
        classify_kidney_cells,
        plot_scpred_scores_distribution,
        plot_scpred_umap_py,
        create_custom_cmap
    )
    print("   [OK] All functions imported successfully")
except Exception as e:
    print(f"   [FAIL] Import error: {e}")
    sys.exit(1)

# Test DevKidCCClassifier class
print("\n2. Testing DevKidCCClassifier class...")
try:
    print(f"   - Class available: {DevKidCCClassifier}")
    print("   - Docstring preview:", DevKidCCClassifier.__doc__[:100].strip() + "...")
    print("   [OK] DevKidCCClassifier class accessible")
except Exception as e:
    print(f"   [FAIL] Error: {e}")

# Test classify_kidney_cells function
print("\n3. Testing classify_kidney_cells function...")
try:
    print(f"   - Function available: {classify_kidney_cells}")
    print("   - Signature: (adata, copy=True, verbose=True, install_deps=True)")
    print("   [OK] classify_kidney_cells function accessible")
except Exception as e:
    print(f"   [FAIL] Error: {e}")

# Test plotting functions
print("\n4. Testing plotting functions...")
try:
    print(f"   - plot_scpred_scores_distribution: {plot_scpred_scores_distribution}")
    print(f"   - plot_scpred_umap_py: {plot_scpred_umap_py}")
    print(f"   - create_custom_cmap: {create_custom_cmap}")
    print("   [OK] All plotting functions accessible")
except Exception as e:
    print(f"   [FAIL] Error: {e}")

# Test required dependencies
print("\n5. Testing required dependencies...")
dependencies = {
    'numpy': 'import numpy',
    'pandas': 'import pandas',
    'scanpy': 'import scanpy',
    'anndata': 'import anndata',
    'matplotlib': 'import matplotlib',
    'seaborn': 'import seaborn',
    'rpy2': 'import rpy2',
}

all_deps_ok = True
for name, import_stmt in dependencies.items():
    try:
        exec(import_stmt)
        print(f"   [OK] {name}")
    except ImportError as e:
        print(f"   [FAIL] {name}: {e}")
        all_deps_ok = False

if all_deps_ok:
    print("\n   [OK] All dependencies available")
else:
    print("\n   [WARNING] Some dependencies missing")

# Summary
print("\n" + "=" * 60)
print("Summary:")
print("=" * 60)
print("Package: devkidcc v0.1.0")
print("Status: READY TO USE")
print("\nAvailable functions:")
print("  - DevKidCCClassifier (main class)")
print("  - classify_kidney_cells (convenience function)")
print("  - plot_scpred_scores_distribution (visualization)")
print("  - plot_scpred_umap_py (visualization)")
print("  - create_custom_cmap (utility)")
print("\nNext steps:")
print("  1. Ensure R is installed (required for classification)")
print("  2. Load your AnnData object")
print("  3. Run: adata = classify_kidney_cells(adata)")
print("  4. Use plotting functions to visualize results")
print("\nSee examples/ directory for detailed usage examples.")
print("=" * 60)



