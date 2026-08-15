#!/usr/bin/env Rscript

# Classify a Seurat object with DevKidCC, from R.
#
# This is the interactive counterpart to /opt/run_dkcc.R: the same call, written
# out step by step so it can be read, adapted, or lifted into a methods section.
# Run it as a script or paste the body into an R session.
#
#   Rscript classify_seurat.R --input organoid.rds --output organoid_DKCC.rds
#
# What it demonstrates
# --------------------
# 1. DevKidCC expects **raw counts** in the counts layer. The models were fit
#    against counts, so pre-normalised input is scored on the wrong scale and
#    fails silently rather than erroring.
# 2. A Seurat **v5** object must be normalised before DKCC(). This is the trap
#    worth knowing: DKCC() reads the `data` layer via Layers(seurat,
#    search = "data") and only rebuilds-and-normalises for legacy v4 assays, so
#    a v5 object handed straight from CreateSeuratObject() is classified on
#    whatever sits in `data` -- nothing, or raw counts. The CLI does this for
#    you; an interactive session must do it itself, as below.
# 3. The labels arrive as two metadata columns:
#      LineageID -- Nephron / Stroma / UrEp / NPC / NPC-like / Endo / unassigned
#      DKCC      -- the finer type within that lineage
#    `unassigned` is a real answer, not a failure: scPred rejects a cell whose
#    probability does not clear the threshold, and organoid data legitimately
#    contains cells outside the reference.

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(DevKidCC)
})

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "Input .rds holding a Seurat object of raw counts"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Where to write the classified .rds"),
    make_option(c("-s", "--summary"), type = "character", default = NULL,
                help = "Where to write a label-count CSV"),
    make_option(c("-t", "--threshold"), type = "double", default = 0.7,
                help = "scPred probability below which a cell is left unassigned [default %default]"),
    make_option(c("-k", "--knn-iter"), type = "integer", default = NULL,
                dest = "knn_iter",
                help = "Rounds of KNN smoothing over unassigned cells; 0 disables it")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("--input is required")

cat("DevKidCC", as.character(packageVersion("DevKidCC")),
    "| Seurat", as.character(packageVersion("Seurat")), "\n")

# ---------------------------------------------------------------------------
# 1. Load and inspect
# ---------------------------------------------------------------------------

seu <- readRDS(opt$input)
if (!inherits(seu, "Seurat")) stop("Input does not hold a Seurat object")

counts <- SeuratObject::LayerData(seu, layer = "counts", assay = "RNA")
cat(sprintf("Input : %d cells x %d genes\n", ncol(seu), nrow(seu)))
cat(sprintf("        counts max %.0f, class %s\n",
            max(counts[, seq_len(min(200, ncol(counts)))]), class(seu[["RNA"]])[1]))

# Log-normalised data tops out near 8-10 and is fractional; counts are
# non-negative integers reaching the hundreds. Refuse the former.
sample_vals <- counts[, seq_len(min(200, ncol(counts)))]
if (max(sample_vals) < 50 && any(sample_vals %% 1 != 0)) {
    stop("This looks log-normalised, not raw counts. DevKidCC normalises ",
         "internally -- pass the raw matrix.")
}

# ---------------------------------------------------------------------------
# 2. Normalise, if the object does not already carry a data layer
# ---------------------------------------------------------------------------
#
# See note 2 in the header. For a v5 assay this is the caller's job.

has_data <- length(SeuratObject::Layers(seu, search = "data")) > 0
if (!has_data) {
    cat("No 'data' layer present -- running NormalizeData() before DKCC()\n")
    seu <- NormalizeData(seu, verbose = FALSE)
} else {
    cat("Existing 'data' layer found; using it\n")
}

# ---------------------------------------------------------------------------
# 3. Classify
# ---------------------------------------------------------------------------

cat("\nClassifying...\n")
dkcc_args <- list(seurat = seu, threshold = opt$threshold)
if (!is.null(opt$knn_iter)) dkcc_args$knn.iter <- opt$knn_iter
seu <- do.call(DKCC, dkcc_args)

md <- seu[[]]
missing <- setdiff(c("DKCC", "LineageID"), colnames(md))
if (length(missing) > 0) {
    stop("Classification returned no ", paste(missing, collapse = "/"), " column")
}

# ---------------------------------------------------------------------------
# 4. Report
# ---------------------------------------------------------------------------

summary_rows <- do.call(rbind, lapply(c("LineageID", "DKCC"), function(column) {
    tab <- sort(table(md[[column]]), decreasing = TRUE)
    tab <- tab[tab > 0]
    data.frame(column  = column,
               label   = names(tab),
               n_cells = as.integer(tab),
               percent = round(100 * as.integer(tab) / nrow(md), 2),
               row.names = NULL)
}))

for (column in c("LineageID", "DKCC")) {
    cat("\n", column, ":\n", sep = "")
    block <- summary_rows[summary_rows$column == column, ]
    for (i in seq_len(nrow(block))) {
        cat(sprintf("  %-16s %6d  %6.2f%%\n",
                    block$label[i], block$n_cells[i], block$percent[i]))
    }
}

assigned <- sum(md$LineageID != "unassigned")
cat(sprintf("\nAssigned: %d/%d (%.1f%%)\n",
            assigned, nrow(md), 100 * assigned / nrow(md)))

if (!is.null(opt$output)) {
    saveRDS(seu, opt$output)
    cat("Wrote", opt$output, "\n")
}
if (!is.null(opt$summary)) {
    write.csv(summary_rows, opt$summary, row.names = FALSE)
    cat("Wrote", opt$summary, "\n")
}
