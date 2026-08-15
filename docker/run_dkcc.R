#!/usr/bin/env Rscript

# Read single-cell files into Seurat, run DevKidCC::DKCC(),
# and save back out as h5ad or RDS.
#
# Supported input:  .h5ad  .h5  .h5seurat  .rds  .RData
# Supported output: .h5ad  .rds
#
# Usage:
#   Rscript /opt/run_dkcc.R \
#       --input sample.h5ad \
#       --output sample_DKCC.h5ad \
#       --format h5ad

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(SeuratDisk)
    library(sceasy)
    library(reticulate)
    library(DevKidCC)
    library(scCustomize)
})

# sceasy reads h5ad through Python's anndata via reticulate. Left to itself,
# reticulate >= 1.41 does not fall back to the system interpreter: finding no
# configured Python it downloads uv, provisions a fresh CPython, installs only
# what the calling package declares, and uses that. sceasy declares no anndata,
# so the conversion lands in an empty environment and fails with
# ModuleNotFoundError -- having never touched this image's environment at all.
#
# The Dockerfile sets RETICULATE_PYTHON, but binding it here too means the script
# is correct however the caller's environment treats that variable. Outside the
# container the path is absent and reticulate behaves as it normally would.
dkcc_python <- Sys.getenv("RETICULATE_PYTHON", "/opt/micromamba/envs/devkid/bin/python")
if (file.exists(dkcc_python)) {
    reticulate::use_python(dkcc_python, required = TRUE)
}

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Input single-cell file"),
    make_option(c("-o", "--output"), type = "character", help = "Output file path"),
    make_option(c("-f", "--format"), type = "character", default = "h5ad",
                help = "Output format: h5ad or rds [default = %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input))  stop("--input is required")
if (is.null(opt$output)) stop("--output is required")

load_rdata_object <- function(path) {
    env  <- new.env()
    objs <- load(path, envir = env)
    for (obj_name in objs) {
        obj <- env[[obj_name]]
        if (inherits(obj, "Seurat")) {
            message("Loaded Seurat object: ", obj_name)
            return(obj)
        }
    }
    stop("No Seurat object found in RData file")
}

read_input <- function(path) {
    ext <- tolower(tools::file_ext(path))
    message("Reading input: ", path)

    if (ext == "h5ad") {
        # SeuratDisk::Convert() can't parse the modern AnnData HDF5 schema
        # (encoding-type/encoding-version attrs from anndata>=0.8 / scanpy>=1.9)
        # and silently produces no output. sceasy reads the file via Python's
        # own anndata module (through reticulate) instead, which handles the
        # modern schema natively.
        obj <- sceasy::convertFormat(path, from = "anndata", to = "seurat", main_layer = "counts")
        # sceasy always builds a legacy (v4) Assay. DKCC() branches on
        # inherits(assay, "Assay5"): when it's not, DKCC() rebuilds the
        # Seurat object from scratch via CreateSeuratObject(), silently
        # dropping any reductions (PCA/UMAP) computed upstream. Converting
        # to Assay5 here makes DKCC() take its layer-joining branch instead,
        # which preserves them.
        obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay5")
    } else if (ext == "h5seurat") {
        obj <- LoadH5Seurat(path)
    } else if (ext == "h5") {
        obj <- CreateSeuratObject(counts = Read10X_h5(path))
    } else if (ext == "rds") {
        obj <- readRDS(path)
    } else if (ext %in% c("rdata", "rdata")) {
        obj <- load_rdata_object(path)
    } else {
        stop("Unsupported file extension: ", ext)
    }

    if (!inherits(obj, "Seurat")) stop("Loaded object is not a Seurat object")
    obj
}

write_output <- function(obj, out_path, out_format) {
    out_format <- tolower(out_format)
    if (out_format == "rds") {
        saveRDS(obj, out_path)
    } else if (out_format == "h5ad") {
        obj@misc  <- list()
        obj@tools <- list()
        as.anndata(x = obj, file_path = dirname(out_path), file_name = basename(out_path))
    } else {
        stop("Output format must be 'h5ad' or 'rds'")
    }
}

# NOTE: this script used to patch DevKidCC::DKCC() here at runtime, rewriting its
# body with deparse()/sub()/assignInNamespace() to guard a missing 'orig.ident'.
# That guard now lives in the package itself (R/DKCC.R), together with a matching
# guard for a PAX2 stripped by the zero-variance filter, so the rewrite is gone.

message("Loading object...")
seu <- read_input(opt$input)

message("Preprocessing...")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)

message("Running DevKidCC...")
seu <- DKCC(seu)

# Check the classification landed before writing. An h5ad that converted cleanly
# but carries no DKCC/LineageID column is the worst failure mode available here:
# exit status 0, an output file of the right size, and labels silently absent.
# Failing here also separates the two causes -- DKCC() not producing the columns
# versus the h5ad writer dropping them -- which is otherwise indistinguishable
# from the outside.
expected_cols <- c("DKCC", "LineageID")
message("  metadata after DKCC(): ", paste(colnames(seu[[]]), collapse = ", "))
missing_cols <- setdiff(expected_cols, colnames(seu[[]]))
if (length(missing_cols) > 0) {
    stop("DKCC() returned no ", paste(missing_cols, collapse = "/"),
         " column; classification did not run as expected.")
}

message("Saving output...")
write_output(seu, opt$output, opt$format)

# And that they survived the conversion, which is a separate library's problem.
# scCustomize::as.anndata builds obs itself, and whatever it does with the Seurat
# metadata it does not reliably carry the classification through -- the columns
# are demonstrably on the object above and absent from the file below.
#
# Rather than depend on that behaviour, write the metadata back over obs from the
# object we already verified. This runs only when something is missing, so a
# conversion that already worked is left untouched.
if (tolower(opt$format) == "h5ad") {
    anndata <- reticulate::import("anndata")
    written <- anndata$read_h5ad(opt$output)
    written_cols <- names(reticulate::py_to_r(written$obs))
    message("  obs written: ", paste(written_cols, collapse = ", "))

    dropped <- setdiff(expected_cols, written_cols)
    if (length(dropped) > 0) {
        message("  h5ad writer dropped ", paste(dropped, collapse = "/"),
                "; restoring obs from the Seurat object")

        md <- seu[[]]
        cell_names <- as.character(reticulate::py_to_r(written$obs_names$tolist()))
        if (!all(cell_names %in% rownames(md))) {
            stop("Cannot restore obs: cell names in the written h5ad do not match ",
                 "the Seurat object.")
        }
        written$obs <- reticulate::r_to_py(md[cell_names, , drop = FALSE])
        written$write_h5ad(opt$output)

        recheck <- names(reticulate::py_to_r(anndata$read_h5ad(opt$output)$obs))
        still_missing <- setdiff(expected_cols, recheck)
        if (length(still_missing) > 0) {
            stop("Classification columns ", paste(still_missing, collapse = "/"),
                 " still absent after restoring obs.")
        }
        message("  [OK] obs restored: ", paste(recheck, collapse = ", "))
    }
}

message("Done.")
