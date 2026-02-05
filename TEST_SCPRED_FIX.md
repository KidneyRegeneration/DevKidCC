# Testing the scPred Compatibility Fix

## What Changed

The workaround has been improved to handle scPred's incompatibility with Seurat v5. The new approach uses **two-part compatibility layer**:

1. **v4-style Assay creation** - Creates v4-compatible assays for scPred
2. **GetAssayData monkey-patch** - Temporarily replaces `GetAssayData()` to convert scPred's old calling convention (`assay = "data"`) to the new Seurat v5 syntax (`layer = "data"`)

## How to Test

### Quick Test
```r
# Reload the package
devtools::load_all("DevKidCC")

# Test with your Seurat v5 object
result <- DKCC(liao, threshold = 0.7)

# Check results
table(result$DKCC)
```

### Expected Behavior

**Before the fix:**
```
Error in `GetAssayData()`:
! `assay` must be one of "RNA", not "data".
```

**After the fix:**
```
Detected 3 data layers. Joining layers for processing...
Performing log-normalization
...
[Classification proceeds normally]
```

## What the Workaround Does

### Part 1: v4 Assay Creation
```r
options(Seurat.object.assay.version = "v3")
seurat <- CreateSeuratObject(data_matrix, meta.data = md)
```

### Part 2: GetAssayData Monkey-Patch
```r
original_GetAssayData <- Seurat::GetAssayData

assignInNamespace("GetAssayData", function(object, assay = NULL, slot = NULL, layer = NULL, ...) {
  # If 'assay' looks like a slot/layer name, convert it
  if (!is.null(assay) && assay %in% c("data", "counts", "scale.data")) {
    layer <- assay
    assay <- NULL
  }
  # Call original with corrected parameters
  original_GetAssayData(object, assay = assay, slot = slot, layer = layer, ...)
}, ns = "Seurat")

# Ensure cleanup even on error
on.exit(assignInNamespace("GetAssayData", original_GetAssayData, ns = "Seurat"), add = TRUE)
```

This allows scPred to work with Seurat v5 by:
1. Creating v4-compatible temporary objects
2. Converting scPred's old `GetAssayData(object, assay = "data")` calls to `GetAssayData(object, layer = "data")`
3. Automatically restoring the original `GetAssayData()` function when done

## Verification

If the fix works, you should see:
- ✅ No `GetAssayData()` errors
- ✅ Classification completes successfully
- ✅ Results returned with DKCC labels

## Troubleshooting

If you still see errors:

1. **Check package is reloaded:**
   ```r
   devtools::load_all("DevKidCC")
   ```

2. **Check Seurat version:**
   ```r
   packageVersion("Seurat")  # Should be >= 5.0.0
   ```

3. **Check if GetAssayData is being called differently:**
   ```r
   # Enable debugging
   options(error = rlang::entrace)
   result <- DKCC(liao)
   ```

4. **Check scPred version:**
   ```r
   packageVersion("scPred")
   # Try updating if needed:
   # devtools::install_github("powellgenomicslab/scPred")
   ```

## Why This Approach Works

The monkey-patch intercepts scPred's calls to `GetAssayData()` and translates them:

| scPred calls (old)             | Translated to (Seurat v5) |
|-------------------------------|---------------------------|
| `GetAssayData(obj, assay="data")` | `GetAssayData(obj, layer="data")` |
| `GetAssayData(obj, assay="counts")` | `GetAssayData(obj, layer="counts")` |
| `GetAssayData(obj, assay="scale.data")` | `GetAssayData(obj, layer="scale.data")` |

This allows scPred to work without modification while DevKidCC handles the compatibility internally.

## Safety

The workaround is safe because:
- ✅ Only affects the current DKCC function call
- ✅ Automatically restores original `GetAssayData()` at the end
- ✅ Uses `on.exit()` to ensure restoration even if errors occur
- ✅ Doesn't modify your input data
- ✅ Doesn't persist changes to the global environment

---

**Questions?** See [SCPRED_COMPATIBILITY.md](SCPRED_COMPATIBILITY.md) for more details.
