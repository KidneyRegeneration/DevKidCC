#!/usr/bin/env Rscript

# Build a Seurat .rds from a MatrixMarket trio (matrix, gene names, cell names).
#
#   Rscript build_rds_from_mtx.R sample.mtx sample.genes.txt sample.cells.txt sample.rds
#
# Why a trio rather than a conversion: the container reads h5ad in Python and
# .rds in R, and deliberately does not convert between them. So the example
# inputs are built for each language from a common, language-neutral source
# instead. MatrixMarket keeps the matrix sparse, where a dense CSV of the same
# data would be an order of magnitude larger.

suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: build_rds_from_mtx.R <mtx> <genes.txt> <cells.txt> <out.rds>")
}
mtx_path <- args[1]; genes_path <- args[2]; cells_path <- args[3]; out_path <- args[4]

counts <- readMM(mtx_path)
genes  <- readLines(genes_path)
cells  <- readLines(cells_path)

stopifnot(nrow(counts) == length(genes), ncol(counts) == length(cells))
rownames(counts) <- genes
colnames(counts) <- cells

# dgCMatrix is what Seurat wants; readMM gives a coordinate-format triplet.
counts <- as(counts, "CsparseMatrix")

seu <- CreateSeuratObject(counts = counts, project = "dkcc_example")
saveRDS(seu, out_path)

cat(sprintf("Wrote %s: %d cells x %d genes, counts max %.0f\n",
            out_path, ncol(seu), nrow(seu), max(counts)))
