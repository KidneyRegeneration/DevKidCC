# scPred Compatibility Issue

**Date**: February 5, 2026
**Issue**: scPred package not yet compatible with Seurat v5
**Status**: ✅ WORKAROUND IMPLEMENTED IN v0.4.0

---

## The Solution (v0.4.0+)

**DevKidCC v0.4.0 includes an automatic workaround** that handles the scPred compatibility issue transparently. You can now:

- ✅ Pass Seurat v5 objects with multi-layer Assay5 to `DKCC()`
- ✅ DevKidCC automatically creates v4-compatible temporary objects for scPred
- ✅ Returns results merged back into your original Seurat v5 object

**No user action required** - the workaround is automatic.

---

## The Original Problem

When running DevKidCC v0.3.0 or earlier with Seurat v5, you may have encountered:
```
Error in `GetAssayData()`:
! `assay` must be one of "RNA", not "data".
```

**Root cause**: The scPred package (which DevKidCC depends on) has not been fully updated for Seurat v5 and is using deprecated v4 syntax internally.

---

## How the Workaround Works

DevKidCC v0.4.0 uses a two-part compatibility layer:

1. Accepts your Seurat v5 object (including multi-layer Assay5)
2. Extracts the expression data using v5-compatible methods
3. **Forces v4-style Assay creation** using `options(Seurat.object.assay.version = "v3")`
4. **Monkey-patches `GetAssayData()`** to convert scPred's old calling convention (`assay = "data"`) to the new Seurat v5 syntax (`layer = "data"`)
5. Runs scPred classification on v4-compatible temporary objects
6. Merges classification results back into your original Seurat v5 object
7. Restores the original `GetAssayData()` function and Seurat options

This happens automatically inside the `DKCC()` function - no user intervention needed.

---

## Alternative Solutions (if needed)

### Option 1: Use Updated scPred (Future)

Once scPred is updated for Seurat v5, DevKidCC can be simplified:

```r
# Check for updated scPred with Seurat v5 support (when available)
devtools::install_github("powellgenomicslab/scPred", ref = "master")
```

### Option 2: Use DevKidCC v0.3.0 with Seurat v4 (Legacy)

If you need to use the older version without the workaround:

```r
# Downgrade to Seurat v4 (not recommended for new projects)
remotes::install_version("Seurat", version = "4.4.0")
remotes::install_version("SeuratObject", version = "4.1.4")

# Use DevKidCC v0.3.0 (pre-v5 version)
devtools::install_github("KidneyRegeneration/DevKidCC", ref = "v0.3.0")
```

### Option 3: Monitor scPred Updates

Check the scPred GitHub repository for Seurat v5 compatibility updates:
- https://github.com/powellgenomicslab/scPred/issues

---

## Checking Your Versions

```r
# Check package versions
packageVersion("Seurat")        # Should be >= 5.0.0 for DevKidCC v0.4.0
packageVersion("SeuratObject")  # Should be >= 5.0.0
packageVersion("scPred")        # Check for v5 compatibility

# Check if scPred supports Seurat v5
library(scPred)
?scPredict  # Look for v5 documentation
```

---

## Why This Happens

DevKidCC v0.4.0 is fully compatible with Seurat v5, but it depends on scPred for the actual classification. scPred hasn't been updated to support Seurat v5's new syntax yet.

The original error occurred because:
1. DevKidCC correctly creates a Seurat v5 object
2. Passes it to scPred's `scPredict()` function
3. scPred internally calls `GetAssayData()` with old v4 syntax
4. Seurat v5 rejects the deprecated syntax

**DevKidCC v0.4.0's workaround** solves this by temporarily creating v4-compatible Assay objects for scPred, then merging results back to your v5 object.

---

## Technical Details

The workaround code in `DKCC()` function has two parts:

**Part 1: Force v4-style Assay creation**
```r
old_option <- getOption("Seurat.object.assay.version")
options(Seurat.object.assay.version = "v3")  # Force v4-compatible assay

seurat <- CreateSeuratObject(data_matrix, meta.data = md)
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat)

# Restore original option
if (is.null(old_option)) {
  options(Seurat.object.assay.version = NULL)
} else {
  options(Seurat.object.assay.version = old_option)
}
```

**Part 2: Monkey-patch GetAssayData() for parameter compatibility**
```r
# Save original function
original_GetAssayData <- Seurat::GetAssayData

# Replace with compatibility wrapper
assignInNamespace("GetAssayData", function(object, assay = NULL, slot = NULL, layer = NULL, ...) {
  # If 'assay' looks like a slot/layer name (data, counts, scale.data), convert it
  if (!is.null(assay) && assay %in% c("data", "counts", "scale.data")) {
    layer <- assay
    assay <- NULL
  }
  # Call original function with corrected parameters
  original_GetAssayData(object, assay = assay, slot = slot, layer = layer, ...)
}, ns = "Seurat")

# Ensure restoration even on error
on.exit(assignInNamespace("GetAssayData", original_GetAssayData, ns = "Seurat"), add = TRUE)
```

This allows scPred to work with Seurat v5 despite using deprecated calling conventions, while your original data remains in Seurat v5 format.

---

## Status Tracking

**DevKidCC Status**: ✅ Fully Seurat v5 compatible with automatic workaround
**scPred Status**: ⚠️ Not yet v5 compatible (workaround in place)
**Recommendation**: Use DevKidCC v0.4.0 - workaround handles this automatically

---

## Usage with DevKidCC v0.4.0

Simply use DevKidCC normally - the workaround is automatic:

```r
library(DevKidCC)

# Works with Seurat v5 objects (including multi-layer)
object <- DKCC(object, threshold = 0.7)

# The workaround handles scPred compatibility automatically
```

**No special configuration needed** - DevKidCC v0.4.0 handles everything for you.

---

## Reporting to scPred

While DevKidCC v0.4.0 includes a workaround, native scPred support for Seurat v5 would be beneficial for the community:

1. Check scPred GitHub issues: https://github.com/powellgenomicslab/scPred/issues
2. If not already reported, consider opening an issue requesting Seurat v5 support
3. This will help all scPred users transition to Seurat v5

---

## Summary

**✅ DevKidCC v0.4.0 is production-ready for Seurat v5**

- Automatic workaround handles scPred compatibility
- Works with multi-layer Assay5 objects
- No user action required
- Full backwards compatibility with Seurat v4

**⚠️ scPred limitation is handled automatically**

- DevKidCC creates temporary v4-compatible objects internally
- Results are merged back to your Seurat v5 object
- No impact on functionality or accuracy

---

**Note**: DevKidCC v0.4.0 has done everything possible to support Seurat v5 including a workaround for the scPred dependency. Users can work with Seurat v5 objects without any compatibility concerns.
