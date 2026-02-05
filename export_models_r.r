# export_models.R
# Script to export DevKidCC models from R to Python-compatible format
#
# This script extracts trained scPred models and exports them in a format
# that can be loaded by scikit-learn in Python

library(DevKidCC)
library(scPred)
library(Seurat)
library(tidyverse)

# Create output directory
output_dir <- "python_models"
dir.create(output_dir, showWarnings = FALSE)

cat("Exporting DevKidCC models for Python...\n")

# ============================================================================
# Step 1: Load the trained reference model from DevKidCC package
# ============================================================================

cat("\n1. Loading reference data...\n")

# Load your reference Seurat object with trained models
# This should be the object you used to train DevKidCC
# Adjust the path to your actual reference data location
data("sysdata", package = "DevKidCC")  # or load from file

# If loading from file:
# reference <- readRDS("path/to/your/trained_reference.rds")

# ============================================================================
# Step 2: Extract feature genes for each classification tier
# ============================================================================

cat("\n2. Extracting feature genes...\n")

# Get scPred metadata
scpred_metadata <- reference@misc$scPred

# Extract genes used for lineage classification
lineage_genes <- scpred_metadata$features

# Extract genes for each cell type model
celltype_genes <- list()
for (lineage in names(scpred_metadata$models)) {
  if (!is.null(scpred_metadata$models[[lineage]])) {
    celltype_genes[[lineage]] <- scpred_metadata$models[[lineage]]$features
  }
}

# Create a data frame with all feature genes
max_genes <- max(
  length(lineage_genes),
  sapply(celltype_genes, length)
)

genes_df <- data.frame(
  lineage = c(lineage_genes, rep(NA, max_genes - length(lineage_genes)))
)

# Add columns for each lineage
for (lineage_name in names(celltype_genes)) {
  genes <- celltype_genes[[lineage_name]]
  genes_df[[lineage_name]] <- c(genes, rep(NA, max_genes - length(genes)))
}

# Save feature genes
write.csv(
  genes_df, 
  file.path(output_dir, "feature_genes.csv"), 
  row.names = FALSE
)

cat(sprintf("  ✓ Saved feature genes for %d tiers\n", ncol(genes_df)))

# ============================================================================
# Step 3: Extract model parameters and training data
# ============================================================================

cat("\n3. Extracting model parameters...\n")

# For each model, we need to extract:
# - Training data (X, y)
# - Model type and parameters
# - Class names

export_model_data <- function(model, tier_name, output_dir) {
  
  if (is.null(model)) {
    cat(sprintf("  Warning: No model found for %s\n", tier_name))
    return(NULL)
  }
  
  # Extract training data
  training_data <- list(
    X = model$trainData,
    y = model$trainRes,
    features = model$features,
    classes = model$classes
  )
  
  # Extract model metadata
  model_info <- list(
    method = model$method,
    params = model$params,
    n_features = length(model$features),
    n_classes = length(model$classes)
  )
  
  # Save to RDS for Python to read
  saveRDS(
    list(
      training_data = training_data,
      model_info = model_info,
      model_object = model
    ),
    file.path(output_dir, paste0(tier_name, "_model_data.rds"))
  )
  
  cat(sprintf("  ✓ Exported %s model (%d features, %d classes)\n", 
              tier_name, model_info$n_features, model_info$n_classes))
  
  return(model_info)
}

# Export lineage model
lineage_model_info <- export_model_data(
  scpred_metadata$lineage_model,  # Adjust based on your actual model structure
  "lineage",
  output_dir
)

# Export cell type models
celltype_model_info <- list()
for (lineage in names(scpred_metadata$models)) {
  celltype_model_info[[lineage]] <- export_model_data(
    scpred_metadata$models[[lineage]],
    paste0("celltype_", lineage),
    output_dir
  )
}

# ============================================================================
# Step 4: Export class names and metadata
# ============================================================================

cat("\n4. Exporting class names...\n")

class_names <- list(
  lineage = levels(reference$LineageID),
  celltypes = list()
)

for (lineage in names(scpred_metadata$models)) {
  if (!is.null(scpred_metadata$models[[lineage]])) {
    class_names$celltypes[[lineage]] <- scpred_metadata$models[[lineage]]$classes
  }
}

saveRDS(
  class_names,
  file.path(output_dir, "class_names.rds")
)

cat(sprintf("  ✓ Saved class names for %d lineages\n", length(class_names$celltypes)))

# ============================================================================
# Step 5: Create summary document
# ============================================================================

cat("\n5. Creating export summary...\n")

summary_text <- sprintf("
DevKidCC Model Export Summary
=============================
Date: %s

Lineage Model:
- Features: %d genes
- Classes: %s

Cell Type Models:
%s

Next Steps:
-----------
1. Run the Python conversion script: python convert_models.py
2. This will convert the .rds files to scikit-learn .joblib format
3. The models will be ready to use with the devkidcc Python package

Files exported:
- feature_genes.csv: Gene lists for each tier
- *_model_data.rds: Model parameters and training data
- class_names.rds: Cell type and lineage labels
",
  Sys.time(),
  length(lineage_genes),
  paste(class_names$lineage, collapse = ", "),
  paste(
    sapply(names(celltype_model_info), function(x) {
      sprintf("  %s: %d features, %d classes", 
              x, 
              celltype_model_info[[x]]$n_features,
              celltype_model_info[[x]]$n_classes)
    }),
    collapse = "\n"
  )
)

writeLines(summary_text, file.path(output_dir, "export_summary.txt"))
cat(summary_text)

cat("\n✓ Export complete! Models saved to:", output_dir, "\n")
cat("\nNext: Run 'python convert_models.py' to convert to Python format\n")
