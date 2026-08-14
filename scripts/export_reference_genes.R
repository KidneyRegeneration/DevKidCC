#!/usr/bin/env Rscript
#' Export the DevKidCC reference gene list used by the Python wrapper.
#'
#' The wrapper restricts the counts matrix it hands to R to these genes: the
#' full matrix produces a multi-gigabyte CSV that R must load whole, and the
#' resulting Seurat object is what gets OOM-killed on large inputs. Every gene
#' outside this list is discarded by scPred anyway.
#'
#' The list is the union of the feature loadings across all seven scPred models
#' shipped in the package (`R/sysdata.rda`), because classification is
#' hierarchical: model1.all assigns a lineage, then a stage-2/3 model refines it
#' using its own feature set. Taking model1.all alone -- 9,977 genes -- misses 67
#' genes that only the downstream models use.
#'
#' Usage:
#'   Rscript scripts/export_reference_genes.R [output.txt]
#'
#' Default output: devkidcc/data/reference_genes.txt

args <- commandArgs(trailingOnly = TRUE)

script_dir <- tryCatch({
  file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", file[1])))
}, error = function(e) ".")

out_path <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(dirname(script_dir), "devkidcc", "data", "reference_genes.txt")
}

suppressPackageStartupMessages({
  library(DevKidCC)
  library(scPred)
})

MODELS <- c("model1.all", "model2.nephron", "model2.stroma", "model2.urep",
            "model3.dn", "model3.pn", "model3.rc")

ns <- asNamespace("DevKidCC")

genes <- character(0)
for (m in MODELS) {
  if (!exists(m, envir = ns, inherits = FALSE)) {
    stop("Model '", m, "' not found in the DevKidCC namespace. ",
         "Installed version: ", as.character(packageVersion("DevKidCC")))
  }
  feats <- rownames(get(m, envir = ns)@feature_loadings)
  cat(sprintf("  %-16s %5d genes\n", m, length(feats)))
  genes <- union(genes, feats)
}

genes <- sort(genes)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writeLines(genes, out_path)

cat("\nDevKidCC version: ", as.character(packageVersion("DevKidCC")), "\n", sep = "")
cat("Wrote ", length(genes), " genes to ", out_path, "\n", sep = "")
