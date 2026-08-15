#!/usr/bin/env Rscript

# Read an R-native single-cell file into Seurat, run DevKidCC::DKCC(),
# and save it back out as RDS.
#
# Supported input:  .h5  .h5seurat  .rds  .RData
# Supported output: .rds
#
# h5ad is deliberately not here. It is AnnData's format, this image contains
# anndata, and reading it from R meant calling back out to Python through
# reticulate and converting in both directions -- which is where every h5ad bug
# in this image has come from. h5ad now goes to /opt/run_dkcc.py; /opt/dkcc
# routes by extension so callers need not choose.
#
# Usage:
#   Rscript /opt/run_dkcc.R \
#       --input sample.rds \
#       --output sample_DKCC.rds

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(DevKidCC)
})

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Input single-cell file"),
    make_option(c("-o", "--output"), type = "character", help = "Output file path"),
    make_option(c("-f", "--format"), type = "character", default = "rds",
                help = "Output format: rds [default = %default]")
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
        # h5ad is AnnData's format and this image reads it in Python, where it
        # is a plain file read rather than a reticulate round trip back out of
        # R. /opt/dkcc routes by extension; this branch exists to say so rather
        # than fail obscurely for anyone calling the R script directly.
        stop("h5ad input is handled by the Python entry point, not this script.\n",
             "  Use:  python /opt/run_dkcc.py --input ", path, " --output <out.h5ad>\n",
             "  Or:   /opt/dkcc --input ", path, " --output <out.h5ad>   (routes by extension)")
    } else if (ext == "h5seurat") {
        # Required only for this branch, so it is loaded here rather than up
        # front: every other format works on an R install without it.
        if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
            stop("Reading .h5seurat needs the SeuratDisk package, which is not installed.")
        }
        obj <- SeuratDisk::LoadH5Seurat(path)
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
        stop("h5ad output is handled by the Python entry point, not this script.\n",
             "  scCustomize::as.anndata built obs itself and silently dropped the\n",
             "  classification columns; writing through anndata avoids the conversion\n",
             "  entirely. Use /opt/run_dkcc.py, or --format rds here.")
    } else {
        stop("Output format must be 'rds' (h5ad goes through /opt/run_dkcc.py)")
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

message("Done.")
