#!/usr/bin/env Rscript
#' Standalone R script to run DevKidCC classification
#' Called by Python wrapper via subprocess
#'
#' Usage: Rscript run_dkcc.R <input.csv> <output.csv> <input.obs.csv>

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript run_dkcc.R <input.csv> <output.csv> <input.obs.csv>")
}

input_csv <- args[1]
output_csv <- args[2]
obs_csv <- args[3]

# Validate input files exist
if (!file.exists(input_csv)) {
  stop(paste("Input file not found:", input_csv))
}
if (!file.exists(obs_csv)) {
  stop(paste("Obs file not found:", obs_csv))
}

cat("============================================================\n")
cat("DevKidCC Classification (R Backend)\n")
cat("============================================================\n\n")

# Load required libraries
cat("Loading required R packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(DevKidCC)
})

cat("  Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("  DevKidCC version:", as.character(packageVersion("DevKidCC")), "\n\n")

# Load count matrix from CSV
cat("Loading expression data from CSV...\n")
counts <- read.csv(input_csv, row.names = 1, check.names = FALSE)
cat("  Matrix shape:", nrow(counts), "genes x", ncol(counts), "cells\n")

# Load obs metadata
cat("Loading metadata...\n")
obs <- read.csv(obs_csv, row.names = 1, check.names = FALSE)
cat("  Metadata:", nrow(obs), "cells x", ncol(obs), "columns\n\n")

# Create Seurat object
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = obs)
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Features:", nrow(seurat_obj), "\n")
cat("  Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n\n")

# Normalize data
cat("Normalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj)
cat("  [OK] Normalization complete\n\n")

# Run DKCC classification
cat("Running DKCC classification...\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()

tryCatch({
  seurat_obj <- DKCC(seurat_obj, threshold = 0.7, max.iter = 1)

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")

  cat("\n[OK] Classification complete (", round(elapsed, 1), " seconds)\n\n")

  # Show results
  if ("LineageID" %in% colnames(seurat_obj[[]])) {
    cat("Lineage assignments:\n")
    print(table(seurat_obj$LineageID, useNA = "ifany"))
    cat("\n")
  }

  if ("DKCC" %in% colnames(seurat_obj[[]])) {
    cat("Cell type assignments (top 10):\n")
    dkcc_counts <- sort(table(seurat_obj$DKCC, useNA = "ifany"), decreasing = TRUE)
    print(head(dkcc_counts, 10))
    cat("\n")
  }

}, error = function(e) {
  cat("\n[ERROR] DKCC classification failed!\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("\nFull traceback:\n")
  traceback()
  stop(paste("DKCC classification failed:", conditionMessage(e)))
})

# Save results metadata to CSV
cat("Saving results...\n")
result_metadata <- seurat_obj[[]]
write.csv(result_metadata, output_csv, row.names = TRUE)
cat("  [OK] Results saved to:", output_csv, "\n\n")

cat("============================================================\n")
cat("[OK] DevKidCC classification complete!\n")
cat("============================================================\n")
