#!/usr/bin/env Rscript

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
} else if (file.exists("../renv/activate.R")) {
  source("../renv/activate.R")
}

# DevKidCC Seurat v5 Migration Helper Script
# This script helps identify and optionally fix compatibility issues

#' Scan R files for Seurat v4 patterns
#' @param directory Path to package R directory
scan_for_issues <- function(directory = "R") {
  
  cat("Scanning for Seurat v4 compatibility issues...\n\n")
  
  # Patterns to find
  patterns <- list(
    slot_usage = list(
      pattern = 'slot\\s*=\\s*["\']',
      replacement = 'layer =',
      description = "Using deprecated 'slot' parameter"
    ),
    direct_slot_access_counts = list(
      pattern = '@assays\\$[A-Za-z]+@counts',
      replacement = 'LayerData(object, layer = "counts")',
      description = "Direct @ access to counts slot"
    ),
    direct_slot_access_data = list(
      pattern = '@assays\\$[A-Za-z]+@data',
      replacement = 'LayerData(object, layer = "data")',
      description = "Direct @ access to data slot"
    ),
    direct_slot_access_scaled = list(
      pattern = '@assays\\$[A-Za-z]+@scale\\.data',
      replacement = 'LayerData(object, layer = "scale.data")',
      description = "Direct @ access to scale.data slot"
    ),
    var_features_slot = list(
      pattern = '@assays\\$[A-Za-z]+@var\\.features',
      replacement = 'VariableFeatures(object, assay = "RNA")',
      description = "Direct @ access to var.features"
    ),
    meta_data_direct = list(
      pattern = '@meta\\.data\\$',
      replacement = '$',
      description = "Using @meta.data$ instead of $"
    ),
    embeddings_direct = list(
      pattern = '@reductions\\$[a-z]+@cell\\.embeddings',
      replacement = 'Embeddings(object, reduction = "...")',
      description = "Direct access to embeddings"
    )
  )
  
  # Get all R files
  r_files <- list.files(directory, pattern = "\\.R$", 
                        full.names = TRUE, recursive = TRUE)
  
  if (length(r_files) == 0) {
    cat("No R files found in", directory, "\n")
    return(invisible(NULL))
  }
  
  cat("Found", length(r_files), "R files to scan\n\n")
  
  # Results storage
  issues <- data.frame(
    file = character(),
    line = integer(),
    issue_type = character(),
    description = character(),
    code_snippet = character(),
    stringsAsFactors = FALSE
  )
  
  # Scan each file
  for (file in r_files) {
    lines <- readLines(file, warn = FALSE)
    
    for (i in seq_along(lines)) {
      line <- lines[i]
      
      # Check each pattern
      for (issue_name in names(patterns)) {
        pattern_info <- patterns[[issue_name]]
        
        if (grepl(pattern_info$pattern, line)) {
          issues <- rbind(issues, data.frame(
            file = basename(file),
            line = i,
            issue_type = issue_name,
            description = pattern_info$description,
            code_snippet = trimws(line),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Print results
  if (nrow(issues) == 0) {
    cat("No compatibility issues found!\n")
  } else {
    cat("Found", nrow(issues), "potential compatibility issues:\n\n")
    
    # Group by issue type
    issue_summary <- table(issues$issue_type)
    cat("Summary by type:\n")
    for (type in names(issue_summary)) {
      cat(sprintf("  - %s: %d occurrences\n", 
                  patterns[[type]]$description,
                  issue_summary[type]))
    }
    
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("Detailed Issues:\n")
    cat(rep("=", 80), "\n\n", sep = "")
    
    for (i in seq_len(nrow(issues))) {
      cat(sprintf("Issue #%d:\n", i))
      cat(sprintf("  File: %s (line %d)\n", 
                  issues$file[i], issues$line[i]))
      cat(sprintf("  Type: %s\n", issues$description[i]))
      cat(sprintf("  Code: %s\n", issues$code_snippet[i]))
      cat("\n")
    }
  }
  
  return(invisible(issues))
}

#' Generate compatibility layer functions
#' @param output_file Where to save the compatibility layer
generate_compatibility_layer <- function(output_file = "R/seurat_compat.R") {
  
  cat("Generating compatibility layer functions...\n")
  
  compat_code <- '# Seurat v5 Compatibility Layer
# Auto-generated compatibility functions for DevKidCC

#\' Check if assay is Seurat v5 (Assay5 class)
#\' @param object Seurat object
#\' @param assay Assay name (default: NULL = DefaultAssay)
#\' @return Logical
#\' @keywords internal
is_v5_assay <- function(object, assay = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  
  if (!assay %in% names(object@assays)) {
    stop(paste("Assay", assay, "not found"))
  }
  
  inherits(object[[assay]], "Assay5")
}

#\' Get expression data compatible with v4 and v5
#\' @param object Seurat object
#\' @param layer Layer/slot name ("data", "counts", or "scale.data")
#\' @param assay Assay name (default: NULL = DefaultAssay)
#\' @return Expression matrix
#\' @keywords internal
#\' @export
get_layer_data <- function(object, layer = "data", assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  
  # Try new method first (v5)
  if (is_v5_assay(object, assay)) {
    data <- Seurat::LayerData(object, layer = layer, assay = assay)
  } else {
    # Fallback to compatible method
    data <- Seurat::GetAssayData(object, layer = layer, assay = assay)
  }
  
  return(data)
}

#\' Set expression data compatible with v4 and v5
#\' @param object Seurat object
#\' @param layer Layer/slot name
#\' @param value New data matrix
#\' @param assay Assay name
#\' @return Updated Seurat object
#\' @keywords internal
#\' @export
set_layer_data <- function(object, layer = "data", value, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  
  if (is_v5_assay(object, assay)) {
    # v5 method
    Seurat::LayerData(object, layer = layer, assay = assay) <- value
  } else {
    # v4 method
    object <- Seurat::SetAssayData(
      object, 
      layer = layer, 
      new.data = value, 
      assay = assay
    )
  }
  
  return(object)
}

#\' Get variable features
#\' @keywords internal
#\' @export
get_var_features <- function(object, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(object)
  Seurat::VariableFeatures(object, assay = assay)
}

#\' Set variable features
#\' @keywords internal
#\' @export
set_var_features <- function(object, features, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(object)
  Seurat::VariableFeatures(object, assay = assay) <- features
  return(object)
}

#\' Get reduction embeddings
#\' @keywords internal
#\' @export
get_embeddings <- function(object, reduction = "umap") {
  if (!reduction %in% names(object@reductions)) {
    stop(paste("Reduction", reduction, "not found"))
  }
  Seurat::Embeddings(object, reduction = reduction)
}

#\' Check Seurat version and warn if issues
#\' @keywords internal
check_seurat_version <- function() {
  version <- packageVersion("Seurat")
  
  if (version < "5.0.0") {
    warning(
      "DevKidCC is optimized for Seurat v5.x. ",
      "You are using v", version, ". ",
      "Some features may not work correctly."
    )
  }
  
  invisible(version)
}

#\' Update old Seurat object to current version
#\' @param object Seurat object
#\' @return Updated object
#\' @export
update_seurat_object <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  old_version <- if (!is.null(object@version)) {
    as.character(object@version)
  } else {
    "unknown"
  }
  
  message("Original object version: ", old_version)
  
  # Update using Seurat function
  if (exists("UpdateSeuratObject", where = asNamespace("Seurat"))) {
    object <- Seurat::UpdateSeuratObject(object)
    message("Updated to Seurat v", packageVersion("Seurat"))
  } else {
    message("UpdateSeuratObject not available")
  }
  
  return(object)
}
'
  
  # Write to file
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(compat_code, output_file)
  
  cat("Compatibility layer written to:", output_file, "\n")
}

#' Update DESCRIPTION file
#' @param description_file Path to DESCRIPTION
update_description <- function(description_file = "DESCRIPTION") {
  
  if (!file.exists(description_file)) {
    cat("Warning: DESCRIPTION file not found at", description_file, "\n")
    return(invisible(NULL))
  }
  
  cat("Updating DESCRIPTION file...\n")
  
  desc <- readLines(description_file, warn = FALSE)
  
  # Update Seurat version requirement
  desc <- gsub(
    "Seurat\\s*\\([^)]*\\)",
    "Seurat (>= 5.0.0)",
    desc
  )
  
  # Update R version in Depends
  desc <- gsub(
    "R\\s*\\([^)]*\\)",
    "R (>= 4.0.0)",
    desc
  )
  
  # Check if SeuratObject is in Imports
  if (!any(grepl("SeuratObject", desc))) {
    # Find Imports line
    imports_line <- grep("^Imports:", desc)
    if (length(imports_line) > 0) {
      # Add SeuratObject after Imports line
      desc <- c(
        desc[1:imports_line],
        "    SeuratObject (>= 5.0.0),",
        desc[(imports_line + 1):length(desc)]
      )
    }
  }
  
  # Backup original
  backup_file <- paste0(description_file, ".backup")
  file.copy(description_file, backup_file, overwrite = TRUE)
  cat("  Backup created:", backup_file, "\n")
  
  # Write updated
  writeLines(desc, description_file)
  cat("DESCRIPTION updated\n")
}

#' Create test file for v5 compatibility
#' @param test_dir Directory for tests
create_tests <- function(test_dir = "tests/testthat") {
  
  cat("Creating compatibility tests...\n")
  
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create test file content as a character vector to avoid quote issues
  test_lines <- c(
    "# Test Seurat v5 compatibility",
    "",
    "library(testthat)",
    "library(Seurat)",
    "library(DevKidCC)",
    "",
    "# Helper: create test object",
    "create_test_object <- function(n_genes = 100, n_cells = 50) {",
    "  set.seed(42)",
    "  counts <- matrix(",
    "    rpois(n_genes * n_cells, lambda = 5),",
    "    nrow = n_genes,",
    "    ncol = n_cells",
    "  )",
    "  rownames(counts) <- paste0('Gene', seq_len(n_genes))",
    "  colnames(counts) <- paste0('Cell', seq_len(n_cells))",
    "  ",
    "  obj <- CreateSeuratObject(counts = counts, project = 'test')",
    "  obj <- NormalizeData(obj, verbose = FALSE)",
    "  obj <- FindVariableFeatures(obj, verbose = FALSE)",
    "  ",
    "  return(obj)",
    "}",
    "",
    "test_that('Compatibility layer works', {",
    "  obj <- create_test_object()",
    "  ",
    "  # Test get_layer_data",
    "  expect_no_error(data <- get_layer_data(obj, 'data'))",
    "  expect_true(is.matrix(data) || inherits(data, 'Matrix'))",
    "  ",
    "  # Test is_v5_assay",
    "  is_v5 <- is_v5_assay(obj)",
    "  expect_type(is_v5, 'logical')",
    "})",
    "",
    "test_that('Variable features handling works', {",
    "  obj <- create_test_object()",
    "  ",
    "  # Get features",
    "  features <- get_var_features(obj)",
    "  expect_true(length(features) > 0)",
    "  ",
    "  # Set features",
    "  new_features <- paste0('Gene', 1:10)",
    "  expect_no_error(obj <- set_var_features(obj, new_features))",
    "  ",
    "  # Verify",
    "  retrieved <- get_var_features(obj)",
    "  expect_equal(retrieved, new_features)",
    "})"
  )
  
  test_file <- file.path(test_dir, "test-seurat-v5-compat.R")
  writeLines(test_lines, test_file)
  
  cat("Test file created:", test_file, "\n")
}

#' Main migration workflow
#' @param package_dir Package root directory
#' @param generate_compat Generate compatibility layer
#' @param update_desc Update DESCRIPTION
#' @param create_test Create test files
run_migration <- function(package_dir = ".",
                          generate_compat = TRUE,
                          update_desc = TRUE,
                          create_test = TRUE) {
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("DevKidCC Seurat v5 Migration Tool\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # Change to package directory
  old_dir <- getwd()
  on.exit(setwd(old_dir))
  setwd(package_dir)
  
  # 1. Scan for issues
  cat("Step 1: Scanning for compatibility issues\n")
  cat(rep("-", 80), "\n", sep = "")
  issues <- scan_for_issues("R")
  
  # 2. Generate compatibility layer
  if (generate_compat) {
    cat("\n")
    cat("Step 2: Generating compatibility layer\n")
    cat(rep("-", 80), "\n", sep = "")
    generate_compatibility_layer()
  }
  
  # 3. Update DESCRIPTION
  if (update_desc) {
    cat("\n")
    cat("Step 3: Updating DESCRIPTION\n")
    cat(rep("-", 80), "\n", sep = "")
    update_description()
  }
  
  # 4. Create tests
  if (create_test) {
    cat("\n")
    cat("Step 4: Creating tests\n")
    cat(rep("-", 80), "\n", sep = "")
    create_tests()
  }
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("Migration Complete!\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  cat("Next steps:\n")
  cat("1. Review the generated compatibility layer (R/seurat_compat.R)\n")
  cat("2. Update your functions to use the compatibility layer:\n")
  cat("   - Replace GetAssayData() with get_layer_data()\n")
  cat("   - Replace SetAssayData() with set_layer_data()\n")
  cat("   - Replace direct @ access with accessor functions\n")
  cat("3. Test your package: devtools::test()\n")
  cat("4. Build documentation: devtools::document()\n")
  cat("5. Check package: devtools::check()\n")
  
  invisible(issues)
}

# If run as script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  package_dir <- if (length(args) > 0) args[1] else "."
  
  run_migration(package_dir)
}
