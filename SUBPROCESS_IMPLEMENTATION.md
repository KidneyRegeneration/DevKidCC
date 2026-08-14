# Subprocess-Based DevKidCC Python Wrapper

## Overview

The DevKidCC Python wrapper has been successfully reimplemented using a subprocess-based approach instead of rpy2. This avoids the critical rpy2/reticulate segfault issue that prevented the original implementation from working.

## Architecture

### Previous Approach (rpy2) - FAILED
- Used rpy2 to call R functions directly from Python
- **Critical Issue**: When Seurat loads, it loads the reticulate package, which conflicts with rpy2 in the same process, causing a segmentation fault
- Even with the latest reticulate version (1.44.1), the segfault persisted

### Current Approach (subprocess) - WORKING ✓
- Python saves data to CSV files (count matrix + metadata)
- Calls standalone R script via subprocess
- R script:
  1. Loads CSVs
  2. Creates Seurat object
  3. Runs DKCC classification
  4. Saves results to CSV
- Python reads results CSV and merges back into AnnData

## Implementation Details

### Files Modified/Created

1. **`devkidcc/classifier_subprocess.py`** - New subprocess-based classifier
   - Replaces rpy2 with subprocess calls
   - Uses CSV for data exchange
   - Fully functional!

2. **`devkidcc/run_dkcc.R`** - Standalone R script
   - Can be called independently
   - Takes CSV inputs, outputs CSV
   - Marshalling only: the Seurat v5 compatibility shim, the zero-variance gene
     filter and the KNN rescue all live in `DevKidCC::DKCC()` (>= 0.5.1). This
     script used to carry its own patched copy of that function; it no longer
     does, and refuses to run against a package too old to have `knn.iter`.

3. **`devkidcc/__init__.py`** - Updated to use subprocess classifier

4. **`tests/test_smoke_synthetic.py`** - Test script
   - Verifies end-to-end functionality
   - Uses synthetic data drawn from the reference gene list

### Data Flow

```
Python (AnnData)
    ↓
  Save to CSV (counts.csv, obs.csv)
    ↓
  Call Rscript run_dkcc.R
    ↓
  R: Load CSVs → Create Seurat → DKCC() → Save results.csv
    ↓
  Python: Load results.csv → Merge into AnnData.obs
    ↓
Python (AnnData with classifications)
```

## Usage

```python
import scanpy as sc
from devkidcc import classify_kidney_cells

# Load your kidney scRNA-seq data
adata = sc.read_h5ad("kidney_data.h5ad")

# Classify
adata = classify_kidney_cells(adata)

# View results
print(adata.obs[['LineageID', 'DKCC']].value_counts())
```

## Testing

```bash
pytest                # fast tests: API surface, gene projection, CSV writer
pytest -m slow -s     # end-to-end through R (synthetic data, and regression)
```

The synthetic-data run classifies nearly everything as `unassigned`, which is
expected — random counts are not kidney cells. What it proves is that the
handoff works in both directions.

## Advantages Over rpy2 Approach

1. **No segfault** - Avoids rpy2/reticulate conflict entirely
2. **Simpler** - No complex rpy2 conversions
3. **Isolated** - R runs in separate process
4. **Debuggable** - Can test R script independently

## Disadvantages

1. **File I/O overhead** - CSV writing/reading takes time
2. **Memory** - the CSV is the binding constraint on large inputs; R loads it
   whole before building a Seurat object. Mitigated by projecting onto the
   DevKidCC reference genes and streaming the file out in gene batches rather
   than through a dense pandas DataFrame. The projection carries full-matrix
   library sizes across in the metadata (`dkcc_full_library_size`) so that
   `LogNormalize` divides by the same denominator it would have without it.
3. **No real-time R access** - Can't call arbitrary R functions on the fly

## Future Plans

**Native Python implementation** - Rewrite classification logic in Python using
scikit-learn/scanpy to eliminate the R dependency entirely. This would remove
the subprocess overhead and the CSV round trip, at the cost of no longer being
guaranteed identical to the published R models.

## Requirements

- Python >= 3.8
- R >= 4.0
- R packages: Seurat (v5), scPred, DevKidCC >= 0.5.1 (SeuratDisk no longer needed)
- Python packages: pandas, numpy, scipy, scanpy, anndata

## Status

✅ **Working and tested** - ready for use with real data.
