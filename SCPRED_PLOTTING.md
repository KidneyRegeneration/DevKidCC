# scPred Score Visualization Functions

DevKidCC v0.4.0 includes R versions of the Python plotting functions for visualizing scPred classification scores.

---

## Functions

### 1. PlotScPredScoresDistribution()

Creates boxplot/violin plots showing the distribution of maximum scPred scores for each predicted cell type.

**Features:**
- Visualize score distributions across cell types
- Optional split by sample or other metadata
- Choose between boxplot, violin plot, or both
- Customizable colors

**Example:**
```r
library(DevKidCC)

# Basic usage
PlotScPredScoresDistribution(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN")
)

# Split by sample
PlotScPredScoresDistribution(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN"),
  split_by = "orig.ident",
  plot_type = "both"
)

# With custom colors
PlotScPredScoresDistribution(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN"),
  colors = c(NPC = "green", EN = "black", DN = "orange")
)
```

---

### 2. PlotScPredUMAPOverlay()

Creates a single UMAP plot where all cells are overlayed and colored by their highest scPred prediction.

**Features:**
- Single unified UMAP showing all cells
- Each cell colored by its predicted cell type
- Automatic label positioning at cluster centroids
- Clean categorical visualization

**Example:**
```r
library(DevKidCC)

# Basic usage
PlotScPredUMAPOverlay(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN"),
  colors = c("green", "black", "orange", "red")
)

# Without labels
PlotScPredUMAPOverlay(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN"),
  colors = c("green", "black", "orange", "red"),
  label = FALSE,
  pt_size = 2,
  alpha = 0.6
)
```

---

### 3. PlotScPredUMAP()

Creates a faceted UMAP visualization where each cell is colored by its maximum scPred score.

**Features:**
- Faceted layout showing each cell type
- Only shows cells where that type has the maximum score
- Custom color gradient per cell type
- Adjustable layout and point parameters

**Example:**
```r
library(DevKidCC)

# Basic usage
PlotScPredUMAP(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN"),
  colors = c("green", "black", "orange", "red")
)

# With custom parameters
PlotScPredUMAP(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN",
                 "scpred_RC", "scpred_UrEp", "scpred_Endo", "scpred_Stroma"),
  colors = c("green", "black", "orange", "red", "blue", "brown", "pink", "yellow"),
  pt_size = 2,
  alpha = 0.6,
  ncol = 4  # 4 columns in facet layout
)
```

---

### 4. PlotScPredUMAPIndividual()

Creates separate UMAP plots for each cell type that can be combined with patchwork.

**Features:**
- Individual plots for fine control
- Each plot has custom color gradient
- Compatible with patchwork for flexible layouts
- Returns list of ggplot objects

**Example:**
```r
library(DevKidCC)
library(patchwork)

# Create individual plots
plots <- PlotScPredUMAPIndividual(
  seurat = organoid,
  score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN"),
  colors = c("green", "black", "orange", "red"),
  pt_size = 1.5,
  alpha = 0.8
)

# Combine with patchwork
wrap_plots(plots, ncol = 2)

# Or custom layout
(plots$NPC | plots$EN) / (plots$DN | plots$PN)

# Access individual plots
plots$NPC + ggtitle("NPC Cells - Custom Title")
```

---

## Complete Workflow Example

```r
library(DevKidCC)
library(Seurat)
library(patchwork)

# Load your Seurat object
seurat <- readRDS("path/to/seurat_object.rds")

# Run DKCC classification
seurat <- DKCC(seurat, threshold = 0.7)

# Define cell types and colors
cell_types <- c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN",
                "scpred_RC", "scpred_UrEp", "scpred_Endo", "scpred_Stroma")
colors <- c("green", "black", "orange", "red", "blue", "brown", "pink", "yellow")

# 1. Score distribution plot
p1 <- PlotScPredScoresDistribution(
  seurat = seurat,
  score_cols = cell_types,
  split_by = "sample",
  plot_type = "both"
)

# 2. Overlay UMAP - single plot with all cell types
p2 <- PlotScPredUMAPOverlay(
  seurat = seurat,
  score_cols = cell_types,
  colors = colors,
  label = TRUE
)

# 3. Faceted UMAP - separate panels for each cell type
p3 <- PlotScPredUMAP(
  seurat = seurat,
  score_cols = cell_types,
  colors = colors,
  ncol = 4
)

# 4. Individual UMAPs - for custom layouts
plots <- PlotScPredUMAPIndividual(
  seurat = seurat,
  score_cols = cell_types,
  colors = colors
)

# Combine individual plots
p4 <- wrap_plots(plots, ncol = 4)

# Save
ggsave("score_distribution.pdf", p1, width = 12, height = 6)
ggsave("umap_overlay.pdf", p2, width = 10, height = 8)
ggsave("umap_faceted.pdf", p3, width = 16, height = 12)
ggsave("umap_individual.pdf", p4, width = 16, height = 12)
```

---

## Notes

**Requirements:**
- ggplot2 (automatically loaded with Seurat)
- dplyr and tidyr (for data manipulation)
- patchwork (optional, for combining individual plots)

**Color Specifications:**
Colors can be specified as:
- Named colors: `"red"`, `"blue"`, `"green"`
- Hex codes: `"#FF0000"`, `"#0000FF"`
- RGB: Using `rgb()` function

**Performance Tips:**
- For large datasets (>50k cells), consider:
  - Reducing `pt_size` (e.g., 0.5)
  - Reducing `alpha` for transparency
  - Using `sample()` to plot a subset of cells for preview

**Customization:**
All functions return ggplot objects that can be further customized:
```r
p <- PlotScPredScoresDistribution(seurat, score_cols)
p +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Custom Title")
```

---

## Comparison with Python Versions

These R functions provide equivalent functionality to the Python functions in DevKidCC-python:

| Python Function | R Function |
|----------------|------------|
| `plot_scpred_scores_distribution()` | `PlotScPredScoresDistribution()` |
| `plot_scpred_umap_py()` | `PlotScPredUMAP()` or `PlotScPredUMAPIndividual()` |
| N/A | `PlotScPredUMAPOverlay()` (R-specific addition) |

**Key Differences:**
- R versions use ggplot2 instead of matplotlib
- R versions are designed for Seurat objects instead of AnnData
- R versions use `patchwork` for combining plots instead of matplotlib subplots
- Color handling is slightly different but functionally equivalent
- R version includes additional `PlotScPredUMAPOverlay()` for categorical overlay visualization

---

For questions or issues, see [NEWS.md](NEWS.md) for changelog.
