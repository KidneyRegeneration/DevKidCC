# DevKidCC Version 0.4.0 - Release Notes

**Release Date**: February 5, 2026
**Major Update**: Full Seurat v5 Compatibility with Multi-Layer Support

---

## 🎯 What's New in 0.4.0

This release brings **complete Seurat v5 compatibility** including critical fixes for multi-layer objects, making DevKidCC work seamlessly with modern single-cell analysis workflows.

### 🔥 Critical Features

#### Multi-Layer Object Support
DevKidCC now automatically handles Seurat v5 objects with split layers - a common scenario after:
- Integration workflows (`IntegrateLayers()`)
- Multi-sample normalization
- SCTransform processing

**What happens**: When DevKidCC detects multiple layers, it automatically joins them and notifies you:
```r
object <- DKCC(object)
# Detected 3 data layers. Joining layers for processing...
```

No user action required - it just works!

### ✨ API Improvements

#### Updated Functions
- `DKCC()` - Main classification function now handles multi-layer objects
- `GeneSummary()` - Gene summary function now handles multi-layer objects
- `join_layers_if_needed()` - New helper function in compatibility layer

### 🐛 Bug Fixes

#### Seurat v5 Compatibility
- ✅ Fixed deprecated `slot=` parameter → `layer=` in FetchData calls
- ✅ Updated `@meta.data` → `[[]]` accessor (modern Seurat v5 syntax)
- ✅ Critical fix: DKCC now works with integrated/multi-sample datasets
- ✅ Critical fix: GeneSummary now works with multi-layer objects

### 📦 Dependencies

Updated to require:
- `Seurat (>= 5.0.0)` - Now explicitly required
- `SeuratObject (>= 5.0.0)` - Already required, now enforced

### 📚 Documentation

New documentation files:
- `MULTI_LAYER_FIX.md` - Technical documentation on multi-layer support
- `SEURAT_V5_MIGRATION.md` - Complete migration guide
- Updated `NEWS.md` - Comprehensive changelog

---

## 🔄 Migration from 0.3.0

### Do I need to change my code?
**No!** Version 0.4.0 is fully backwards compatible.

### What if I'm using Seurat v4?
Your code will continue to work exactly as before.

### What if I'm using Seurat v5 with single-layer objects?
Everything works the same, no changes needed.

### What if I have multi-layer Seurat v5 objects?
**This is what 0.4.0 fixes!** Your previously failing code will now work automatically.

---

## 💻 Installation

### From GitHub (Latest)
```r
# Install/update DevKidCC
devtools::install_github("KidneyRegeneration/DevKidCC", ref = "main")
```

### Local Installation
```r
# If you have the source code
devtools::install("path/to/DevKidCC")
```

---

## 🧪 Testing Your Upgrade

### Quick Test
```r
library(DevKidCC)
library(Seurat)

# Load your Seurat v5 object (single or multi-layer)
object <- readRDS("your_kidney_data.rds")

# Run classification
object <- DKCC(object)

# Check results
table(object$LineageID)
table(object$DKCC)
```

### With Multi-Layer Object
```r
# If you have an integrated object
object <- IntegrateLayers(object, method = CCAIntegration)

# This now works in v0.4.0!
object <- DKCC(object)
# You'll see: "Detected X data layers. Joining layers for processing..."
```

---

## 🔍 What Changed Under the Hood

### Code Changes

#### DKCC.R (lines 28-43)
Added automatic layer detection and joining before data extraction:
```r
if (inherits(old.seurat[["RNA"]], "Assay5")) {
  layers <- SeuratObject::Layers(old.seurat, search = "data")
  if (length(layers) > 1) {
    message("Detected ", length(layers), " data layers. Joining layers...")
    old.seurat[["RNA"]] <- JoinLayers(old.seurat[["RNA"]])
  }
}
data_matrix <- get_layer_data(old.seurat, layer = "data", assay = "RNA")
```

#### gene_summary.R (lines 36-42)
Same layer joining logic for GeneSummary function

#### seurat_compat.R (lines 150-171)
New `join_layers_if_needed()` helper function for reusable layer management

---

## 🆚 Comparison with 0.3.0

| Feature | 0.3.0 | 0.4.0 |
|---------|-------|-------|
| Seurat v4 support | ✅ | ✅ |
| Seurat v5 single-layer | ⚠️ Partial | ✅ Full |
| Seurat v5 multi-layer | ❌ Broken | ✅ Fixed |
| Integrated datasets | ❌ Failed | ✅ Works |
| Multi-sample objects | ❌ Failed | ✅ Works |
| SCTransform objects | ⚠️ Sometimes | ✅ Always |

---

## 📖 Additional Resources

- **Technical Details**: See `MULTI_LAYER_FIX.md`
- **Migration Guide**: See `SEURAT_V5_MIGRATION.md`
- **Package Website**: https://kidneyregeneration.github.io/DevKidCC/
- **Publication**: Wilson et al., 2022, Genome Medicine

---

## 🙏 Acknowledgments

This update ensures DevKidCC works with all modern Seurat v5 workflows, making it compatible with the latest single-cell analysis standards.

---

## 📝 Full Changelog

### Added
- Multi-layer object detection and automatic joining
- `join_layers_if_needed()` compatibility helper function
- Comprehensive documentation (MULTI_LAYER_FIX.md, SEURAT_V5_MIGRATION.md)
- User-friendly messages when joining layers

### Changed
- DKCC() now handles multi-layer Assay5 objects
- GeneSummary() now handles multi-layer Assay5 objects
- Updated all Seurat v4 syntax to v5 equivalents
- Explicit Seurat >= 5.0.0 dependency

### Fixed
- Critical bug: DKCC failing with multi-layer objects
- Critical bug: GeneSummary failing with multi-layer objects
- Deprecated `slot=` parameter in FetchData
- Direct `@meta.data` accessor usage

### Removed
- None (fully backwards compatible)

---

**Upgrade recommended for all users working with Seurat v5!**
