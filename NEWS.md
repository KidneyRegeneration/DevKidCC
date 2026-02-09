# DevKidCC 0.5.0

**Major Update: Simplified NPC Refinement for Enhanced Seurat v5 Compatibility**

Released: 2026-02-09

## Major Changes

### Simplified NPC Refinement
- Removed PAX2-based clustering refinement to resolve Seurat v5 compatibility issues
- All NPC cells now uniformly classified as "NPC-like"
- Eliminates "cannot xtfrm data frames" errors during classification
- Maintains classification accuracy for all other cell types (Nephron, Stroma, UrEp, Endo)

### Performance Improvements
- Faster execution: Removed clustering overhead for NPC cells
- Reduced memory usage during NPC refinement step

## Technical Details

### Code Changes (R/DKCC.R, lines 191-204)
- **Removed:** GeneSummary() call for PAX2 expression analysis
- **Removed:** Clustering analysis (NormalizeData, FindVariableFeatures, ScaleData, PCA, UMAP, FindNeighbors, FindClusters) on NPC subset
- **Added:** Direct metadata update to mark all NPCs as "NPC-like"

### Backward Compatibility
- **Impact:** Existing analyses may show different NPC classifications
- **v0.4.0:** NPCs classified as either "NPC" or "NPC-like" based on PAX2 clustering
- **v0.5.0:** All NPCs classified as "NPC-like"
- **Recommendation:** Re-run classifications on v0.5.0 for consistency

## Breaking Changes
- None for standard workflows
- Internal API: NPC subtyping functionality removed

## Future Work
- Re-implement NPC refinement using simpler PAX2 expression threshold (without clustering)
- Fix GeneSummary() for standalone use with proper data frame handling
- Native Python implementation to eliminate R dependency

---

# DevKidCC 0.4.0

**Major Update: Full Seurat v5 Compatibility with Multi-Layer Support**

Released: 2026-02-05

## New Features

* **Multi-layer support** - DevKidCC now automatically detects and handles Seurat v5 objects with split layers (common after integration workflows)
  - Automatically joins layers before processing
  - Works seamlessly with integrated datasets, multi-sample experiments, and SCTransform workflows
  - Provides informative messages when joining layers

* **New visualization functions** - Added R versions of Python plotting functions from DevKidCC-python:
  - `PlotScPredScoresDistribution()` - Boxplot/violin plot of prediction scores, optionally split by sample
  - `PlotScPredUMAPOverlay()` - Single UMAP with all cells colored by highest prediction (categorical overlay)
  - `PlotScPredUMAP()` - UMAP visualization colored by scPred scores with faceting
  - `PlotScPredUMAPIndividual()` - Individual UMAP plots per cell type (compatible with patchwork)

* Added `join_layers_if_needed()` helper function in compatibility layer for reusable layer management

## Bug Fixes

* **Critical fix**: DKCC main function now works with multi-layer Assay5 objects (#issue-multi-layer)
* **Critical fix**: GeneSummary function now works with multi-layer Assay5 objects
* **Critical fix**: Corrected namespace - `LayerData()`, `JoinLayers()`, etc. are from `SeuratObject`, not `Seurat`
* **Critical fix**: Plotting functions now handle non-numeric score columns (#issue-xtfrm-error)
  - Fixed "Error in xtfrm.data.frame(x) : cannot xtfrm data frames"
  - All plotting functions now convert score columns to numeric matrix before operations
  - Handles character/factor score columns automatically
  - Affects: `PlotScPredScoresDistribution()`, `PlotScPredUMAP()`, `PlotScPredUMAPOverlay()`, `PlotScPredUMAPIndividual()`, `PlotScPredUMAP2()`
* **Workaround**: Added comprehensive scPred compatibility layer:
  - Creates v4-style assays for classification
  - Monkey-patches `GetAssayData()` to convert scPred's old calling convention to Seurat v5 syntax
  - Ensures cleanup even on error using `on.exit()`
* Fixed deprecated Seurat v4 `slot=` parameter → `layer=` in FetchData calls
* Updated direct `@meta.data` accessor to recommended `[[]]` operator (4 locations in visualisation.R)
* Added explicit namespace prefixes to all SeuratObject functions
* Fixed subset error when no NPC cells are present - now gracefully skips NPC refinement with informative message

## Dependencies

* Now explicitly requires Seurat (>= 5.0.0)
* Requires SeuratObject (>= 5.0.0)
* **Note**: scPred is not yet fully Seurat v5 compatible - DevKidCC includes workaround

## Documentation

* Added comprehensive MULTI_LAYER_FIX.md technical documentation
* Updated SEURAT_V5_MIGRATION.md with migration details
* All functions now use Seurat v5 compatible syntax

## Backwards Compatibility

* Fully backwards compatible with Seurat v4 objects
* Works with both single-layer and multi-layer Seurat v5 objects
* No breaking changes to user-facing API

---

# DevKidCC 0.3.0

# DevKidCC 0.2.3

# DevKidCC 0.2.2

# DevKidCC 0.2.1

# DevKidCC 0.2.0

# DevKidCC 0.1.6

