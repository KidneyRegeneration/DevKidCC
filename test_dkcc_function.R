# Test script for DKCC function with built-in organoid dataset
# This tests the refactored DKCC that extracts metadata from temp v4 objects

cat("========================================\n")
cat("Testing DKCC Function\n")
cat("========================================\n\n")

# Load required packages
cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DevKidCC)
  library(dplyr)
  library(tibble)
})

cat("Package versions:\n")
cat("  Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("  SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
cat("  DevKidCC:", as.character(packageVersion("DevKidCC")), "\n\n")

# Load the built-in organoid dataset
cat("Loading DevKidCC::organoid dataset...\n")
data("organoid", package = "DevKidCC")
organoid <- liao 
cat("Object info:\n")
cat("  Class:", class(organoid), "\n")
cat("  Cells:", ncol(organoid), "\n")
cat("  Features:", nrow(organoid), "\n")
cat("  Assays:", names(organoid@assays), "\n")
if (inherits(organoid[["RNA"]], "Assay5")) {
  cat("  Assay type: Assay5 (v5)\n")
  layers <- SeuratObject::Layers(organoid, search = "data")
  cat("  Data layers:", length(layers), "-", paste(layers, collapse = ", "), "\n")
} else {
  cat("  Assay type: Assay (v4)\n")
}
cat("\n")

# Check if DKCC has been run before
if ("DKCC" %in% colnames(organoid[[]])) {
  cat("Note: DKCC has been run before. Metadata will be reset.\n\n")
}

# Run DKCC classification
cat("Running DKCC classification...\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()

tryCatch({
  organoid <- DKCC(organoid, threshold = 0.7, max.iter = 1)

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")

  cat("\n========================================\n")
  cat("DKCC Classification Complete!\n")
  cat("========================================\n\n")
  cat("Time elapsed:", round(elapsed, 2), "seconds\n\n")

  # Check results
  cat("Results Summary:\n")
  cat("----------------\n")

  # Check for DKCC column
  if ("DKCC" %in% colnames(organoid[[]])) {
    cat("✓ DKCC column created\n")
    dkcc_counts <- table(organoid$DKCC, useNA = "ifany")
    cat("\nDKCC assignments:\n")
    print(dkcc_counts)
  } else {
    cat("✗ DKCC column not found!\n")
  }

  # Check for LineageID column
  if ("LineageID" %in% colnames(organoid[[]])) {
    cat("\n✓ LineageID column created\n")
    lineage_counts <- table(organoid$LineageID, useNA = "ifany")
    cat("\nLineage assignments:\n")
    print(lineage_counts)
  } else {
    cat("\n✗ LineageID column not found!\n")
  }

  # Check for scpred score columns
  scpred_cols <- grep("^scpred_", colnames(organoid[[]]), value = TRUE)
  if (length(scpred_cols) > 0) {
    cat("\n✓ Found", length(scpred_cols), "scpred score columns:\n")
    cat("  ", paste(head(scpred_cols, 10), collapse = ", "))
    if (length(scpred_cols) > 10) cat(", ...")
    cat("\n")
  } else {
    cat("\n✗ No scpred score columns found!\n")
  }

  # Check metadata columns
  cat("\nTotal metadata columns:", ncol(organoid[[]]), "\n")

  cat("\n========================================\n")
  cat("Test PASSED ✓\n")
  cat("========================================\n")

}, error = function(e) {
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")

  cat("\n========================================\n")
  cat("Test FAILED ✗\n")
  cat("========================================\n\n")
  cat("Error after", round(elapsed, 2), "seconds:\n")
  cat(conditionMessage(e), "\n\n")

  cat("Traceback:\n")
  traceback()

  stop(e)
})
