#!/usr/bin/env Rscript
#' Standalone R script to run DevKidCC classification
#' Called by the Python wrapper via subprocess
#'
#' Usage: Rscript run_dkcc.R <input.csv> <output.csv> <input.obs.csv> [threshold] [max_iter] [knn_iter]
#'
#' Everything that used to be patched in here -- the Seurat v5 GetAssayData
#' shim, the zero-variance gene filter, the PCA/UMAP needed by KNN smoothing,
#' and the KNN rescue itself -- now lives in DevKidCC::DKCC(). This script only
#' marshals data across the process boundary.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 6) {
  stop("Usage: Rscript run_dkcc.R <input.csv> <output.csv> <input.obs.csv> [threshold] [max_iter] [knn_iter]")
}

input_csv  <- args[1]
output_csv <- args[2]
obs_csv    <- args[3]

threshold <- if (length(args) >= 4) as.numeric(args[4])  else 0.7
max_iter  <- if (length(args) >= 5) as.integer(args[5])  else 1L
knn_iter  <- if (length(args) >= 6) as.integer(args[6])  else 20L

if (is.na(threshold)) stop("threshold must be numeric, got: ", args[4])
if (is.na(max_iter))  stop("max_iter must be an integer, got: ", args[5])
if (is.na(knn_iter))  stop("knn_iter must be an integer, got: ", args[6])

if (!file.exists(input_csv)) stop(paste("Input file not found:", input_csv))
if (!file.exists(obs_csv))   stop(paste("Obs file not found:", obs_csv))

cat("============================================================\n")
cat("DevKidCC Classification (R Backend)\n")
cat("============================================================\n\n")

cat("Loading required R packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(DevKidCC)
})

cat("  Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("  DevKidCC version:", as.character(packageVersion("DevKidCC")), "\n")

# knn.iter arrived in DevKidCC 0.5.1. Fail loudly rather than silently running
# an older DKCC() that ignores it -- the no-smoothing arm would otherwise be
# reported as if it had run when it had not.
if (!"knn.iter" %in% names(formals(DevKidCC::DKCC))) {
  stop("Installed DevKidCC::DKCC() has no knn.iter parameter. ",
       "Upgrade to >= 0.5.1: remotes::install_github('KidneyRegeneration/DevKidCC')")
}
cat("\n")

cat("Loading expression data from CSV...\n")
counts <- read.csv(input_csv, row.names = 1, check.names = FALSE)
cat("  Matrix shape:", nrow(counts), "genes x", ncol(counts), "cells\n")

cat("Loading metadata...\n")
obs <- read.csv(obs_csv, row.names = 1, check.names = FALSE)
cat("  Metadata:", nrow(obs), "cells x", ncol(obs), "columns\n\n")

cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = obs)
rm(counts)
invisible(gc())
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Features:", nrow(seurat_obj), "\n")
cat("  Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n\n")

cat("Normalizing data...\n")
# The Python wrapper sends only the DevKidCC reference genes to keep the CSV (and
# so R's peak memory) manageable. NormalizeData() would then divide each cell by
# the sum over *those* genes, which is not the library size the models were fit
# against -- it shifts every normalised value and every scPred score with it. The
# wrapper therefore ships the full-matrix totals in this metadata column, and we
# reproduce LogNormalize by hand from them: log1p(count / total * 10000), which is
# exactly what NormalizeData() computes, only with the right denominator.
#
# Absent (a caller driving this script directly with an unprojected matrix), fall
# back to stock NormalizeData.
lib_col <- "dkcc_full_library_size"
if (lib_col %in% colnames(seurat_obj[[]])) {
  totals <- as.numeric(seurat_obj[[lib_col]][, 1])
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Column '", lib_col, "' contains non-positive or non-finite library sizes.")
  }
  cat("  Using full-matrix library sizes from '", lib_col, "'\n", sep = "")
  cts <- SeuratObject::LayerData(seurat_obj, layer = "counts")
  # Column scaling via a diagonal matrix rather than `cts / totals`: sparse
  # arithmetic against a plain vector recycles column-major over every entry,
  # including the structural zeros, which densifies the matrix.
  norm <- cts %*% Matrix::Diagonal(x = 1e4 / totals)
  dimnames(norm) <- dimnames(cts)
  if (methods::is(norm, "sparseMatrix")) {
    norm@x <- log1p(norm@x)          # log1p(0) == 0, so the zeros stay structural
  } else {
    norm <- log1p(norm)
  }
  seurat_obj <- SeuratObject::SetAssayData(seurat_obj, layer = "data", new.data = norm)
  rm(cts, norm)
  invisible(gc())
} else {
  cat("  No '", lib_col, "' column; using NormalizeData() column sums\n", sep = "")
  seurat_obj <- NormalizeData(seurat_obj)
}
cat("  [OK] Normalization complete\n\n")

cat("Running DKCC classification...\n")
cat("  Parameters:\n")
cat("    - threshold:", threshold, "\n")
cat("    - max.iter: ", max_iter, "\n")
cat("    - knn.iter: ", knn_iter,
    if (knn_iter == 0) " (KNN smoothing disabled)" else "", "\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()

tryCatch({
  seurat_obj <- DKCC(seurat_obj, threshold = threshold, max.iter = max_iter,
                     knn.iter = knn_iter)

  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat("\n[OK] Classification complete (", round(elapsed, 1), " seconds)\n\n")

  if ("LineageID" %in% colnames(seurat_obj[[]])) {
    cat("Lineage assignments:\n")
    print(table(seurat_obj$LineageID, useNA = "ifany"))
    cat("\n")
  }

  if ("DKCC" %in% colnames(seurat_obj[[]])) {
    cat("Cell type assignments (top 10):\n")
    print(head(sort(table(seurat_obj$DKCC, useNA = "ifany"), decreasing = TRUE), 10))
    cat("\n")
  }

}, error = function(e) {
  cat("\n[ERROR] DKCC classification failed!\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("\nFull traceback:\n")
  traceback()
  stop(paste("DKCC classification failed:", conditionMessage(e)))
})

cat("Saving results...\n")
write.csv(seurat_obj[[]], output_csv, row.names = TRUE)
cat("  [OK] Results saved to:", output_csv, "\n\n")

cat("============================================================\n")
cat("[OK] DevKidCC classification complete!\n")
cat("============================================================\n")
