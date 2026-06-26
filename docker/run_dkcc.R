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
    library(DevKidCC)
    library(scCustomize)
})

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
        tmp_h5s <- sub("\\.h5ad$", ".h5seurat", tempfile())
        Convert(path, dest = "h5seurat", overwrite = TRUE, filename = tmp_h5s)
        obj <- LoadH5Seurat(tmp_h5s)
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

message("Saving output...")
write_output(seu, opt$output, opt$format)

message("Done.")
