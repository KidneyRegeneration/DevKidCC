#!/usr/bin/env Rscript

# Dump per-cell DevKidCC labels from a classified .rds to CSV.
#
#   Rscript export_labels.R sample_DKCC.rds sample_labels.csv
#
# Exists so the R and Python routes can be compared cell by cell without either
# language having to read the other's file format.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: export_labels.R <classified.rds> <out.csv>")

seu <- readRDS(args[1])
md  <- seu[[]]
stopifnot(all(c("DKCC", "LineageID") %in% colnames(md)))

out <- data.frame(cell      = rownames(md),
                  LineageID = as.character(md$LineageID),
                  DKCC      = as.character(md$DKCC),
                  row.names = NULL)
write.csv(out, args[2], row.names = FALSE)
cat("Wrote", args[2], "-", nrow(out), "cells\n")
