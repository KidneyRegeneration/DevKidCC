#' DKCC
#'
#' @param seurat seurat object
#' @param threshold minimum value for an identity to be assigned within the model call, default is 0.7
#' @param max.iter Can ask scPred to run this number of integrations, set to 0 be default
#' @param knn.iter Maximum KNN smoothing iterations for unassigned cell rescue. Set to 0 to disable KNN smoothing (useful for benchmarking). Default is 20.
#'
#' @return seurat object with additional metadata columns
#' @export
#'
#' @aliases DKCC
#'
#' @examples
#' organoid <- DKCC(organoid)
#' organoid_no_knn <- DKCC(organoid, knn.iter = 0)

DKCC <- function(seurat, threshold = 0.7, max.iter = 1, knn.iter = 20) {

  if (("dkcc" %in% colnames(seurat[[]])) == FALSE){
    seurat@misc$old.meta <- seurat[[]]
  } else {
    seurat[[]] <- seurat@misc$old.meta
  }
  md <- seurat[[]]

  # Handle Seurat v5 multi-layer objects
  DefaultAssay(seurat) <- "RNA"

  if (inherits(seurat[["RNA"]], "Assay5")) {
    # Join layers if needed
    layers <- SeuratObject::Layers(seurat, search = "data")
    if (length(layers) > 1) {
      message("Joining ", length(layers), " data layers for processing...")
      seurat[["RNA"]] <- SeuratObject::JoinLayers(seurat[["RNA"]])
    }
  } else {
    # For v4 objects, create a new v5-compatible object
    message("Converting v4 Assay to v5 format...")

    # Clean metadata - convert factors to characters to prevent xtfrm errors
    md_clean <- md
    for (col in colnames(md_clean)) {
      if (is.factor(md_clean[[col]])) {
        md_clean[[col]] <- as.character(md_clean[[col]])
      }
    }

    # Create new object with cleaned metadata
    seurat <- Seurat::CreateSeuratObject(
      counts = Seurat::GetAssayData(seurat, assay = "RNA"),
      meta.data = md_clean
    )
    seurat <- Seurat::NormalizeData(seurat)
  }

  # Patch GetAssayData for scPred compatibility
  # scPred uses old v4 parameter syntax, we translate it to v5 on-the-fly
  trace(
    SeuratObject:::GetAssayData.Seurat,
    tracer = quote({
      if (!missing(assay) && is.character(assay) && assay %in% c("data", "counts", "scale.data")) {
        if (missing(layer) || is.null(layer)) {
          layer <- assay
          assay <- NULL
        }
      }
      if (!missing(slot) && is.character(slot) && slot %in% c("data", "counts", "scale.data")) {
        if (missing(layer) || is.null(layer)) {
          layer <- slot
          slot <- NULL
        }
      }
    }),
    print = FALSE,
    where = asNamespace("SeuratObject")
  )

  on.exit({
    try(untrace(SeuratObject:::GetAssayData.Seurat, where = asNamespace("SeuratObject")), silent = TRUE)
  }, add = TRUE)

  # Now run the standard DKCC classification workflow on v4 object
  dkcc <- data.frame()

  # Remove zero-variance genes before classification.
  # scPred standardises its feature genes internally; any gene with zero variance
  # across the cells in this object will produce NA/NaN/Inf and crash the call.
  # This filters the RNA assay in-place so all downstream scPredict calls are safe.
  message("Filtering zero-variance genes...")
  mat <- SeuratObject::LayerData(seurat, assay = "RNA", layer = "data")
  gene_vars <- Matrix::rowMeans(mat * mat) - Matrix::rowMeans(mat)^2
  genes_keep <- rownames(mat)[!is.na(gene_vars) & gene_vars > 0]
  n_removed <- nrow(mat) - length(genes_keep)
  if (n_removed > 0) {
    message("  Removed ", n_removed, " zero-variance gene(s) from ", nrow(mat), " total.")
    seurat <- seurat[genes_keep, ]
  } else {
    message("  No zero-variance genes found.")
  }

  # Step 1: Lineage classification
  message("Running lineage classification...")
  seurat <- suppressWarnings(scPred::scPredict(seurat, reference = model1.all, threshold = threshold, max.iter.harmony = max.iter))
  seurat$scpred_prediction <- gsub("Endothelial", "Endo", seurat$scpred_prediction)
  colnames(seurat[[]]) <- gsub("Endothelial", "Endo", colnames(seurat[[]]))
  seurat$LineageID <- seurat$scpred_prediction
  seurat$LineageID_max <- seurat$scpred_max

  # KNN smoothing requires a low-dimensional embedding (UMAP by default).
  # Compute PCA + UMAP if they are absent and KNN is requested.
  if (knn.iter > 0) {
    if (!"umap" %in% Reductions(seurat)) {
      message("Computing PCA and UMAP for KNN smoothing...")
      n_pcs <- min(30, ncol(seurat) - 1)
      seurat <- suppressWarnings(
        seurat %>%
          FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
          ScaleData(verbose = FALSE) %>%
          RunPCA(npcs = n_pcs, verbose = FALSE) %>%
          RunUMAP(dims = 1:n_pcs, verbose = FALSE)
      )
      message("  UMAP computed (", n_pcs, " PCs)")
    }
    seurat <- fill_unassigned_by_knn_seurat(seurat, "LineageID", threshold = 0.4, k=25, max_iter=knn.iter)
  } else {
    # knn.iter = 0 is the documented way to disable smoothing for benchmarking.
    # It cannot go through fill_unassigned_by_knn_seurat(): that function stops
    # if the umap reduction is absent, which -- since the block above no longer
    # computes one -- is exactly the case knn.iter = 0 produces.
    message("KNN smoothing disabled (knn.iter = 0).")
  }

  dkcc <- seurat[[]] %>% rownames_to_column("cell") %>% filter(LineageID %in% c("unassigned", "NPC", "Endo")) %>% transmute(cell = cell, dkcc = LineageID)

  # Accumulate cell metadata without merging Seurat objects.
  # merge.Seurat in v5 calls merge.Assay5 -> LayerData.Assay5 which crashes when
  # merging subsets that went through different scPredict runs. Since t1 is only
  # ever used for t1[[]] (metadata) at the end, we collect metadata as a plain
  # data.frame and bind_rows instead.
  t1_meta <- seurat[[]] %>% rownames_to_column("cell") %>% filter(LineageID %in% c("unassigned", "NPC", "Endo"))

  # Step 2: Nephron classification
  if (nrow(seurat[[]] %>% filter(LineageID == "Nephron")) > 2) {
    message("Running nephron classification...")
    nephronid <- suppressWarnings(scPred::scPredict(seurat[, seurat$LineageID == "Nephron"], reference = model2.nephron,
                                    threshold = 0.0, max.iter.harmony = max.iter))
    nephronid$NephronID <- nephronid$scpred_prediction

    dkcc <- bind_rows(dkcc,
                      nephronid[[]] %>% rownames_to_column("cell") %>% filter(NephronID %in% c("EN")) %>% transmute(cell = cell, dkcc = NephronID))

    if (nrow(nephronid[[]] %>% filter(NephronID %in% c("EN"))) > 0) {
      t1_meta <- bind_rows(t1_meta,
                           nephronid[[]] %>% rownames_to_column("cell") %>% filter(NephronID %in% c("EN")))
    }

    # Step 2a: Proximal nephron
    if (nrow(nephronid[[]] %>% filter(NephronID == "PN")) > 2) {
      message("Running proximal nephron classification...")
      proximalid <- suppressWarnings(scPred::scPredict(nephronid[, nephronid$NephronID == "PN"],
                                       reference = model3.pn, threshold = 0.0, max.iter.harmony = max.iter))
      proximalid$SegmentID <- proximalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        proximalid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1_meta <- bind_rows(t1_meta, proximalid[[]] %>% rownames_to_column("cell"))
    }

    # Step 2b: Distal nephron
    if (nrow(nephronid[[]] %>% filter(NephronID == "DN")) > 2) {
      message("Running distal nephron classification...")
      distalid <- suppressWarnings(scPred::scPredict(nephronid[, nephronid$NephronID == "DN"],
                                     reference = model3.dn, threshold = 0.0, max.iter.harmony = max.iter))
      distalid$SegmentID <- distalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        distalid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1_meta <- bind_rows(t1_meta, distalid[[]] %>% rownames_to_column("cell"))
    }

    # Step 2c: Renal corpuscle
    if (nrow(nephronid[[]] %>% filter(NephronID == "RC")) > 2) {
      message("Running renal corpuscle classification...")
      rcid <- suppressWarnings(scPred::scPredict(nephronid[, nephronid$NephronID == "RC"], reference = model3.rc,
                                 threshold = 0.0, max.iter.harmony = max.iter))
      rcid$SegmentID <- rcid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        rcid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1_meta <- bind_rows(t1_meta, rcid[[]] %>% rownames_to_column("cell"))
    }
  }

  # Step 3: Stroma classification
  if (nrow(seurat[[]] %>% filter(LineageID == "Stroma")) > 2) {
    message("Running stroma classification...")
    stromaid <- suppressWarnings(scPred::scPredict(seurat[, seurat$LineageID == "Stroma"], reference = model2.stroma,
                                   threshold = 0.0, max.iter.harmony = max.iter))
    stromaid$StromaID <- stromaid$scpred_prediction

    dkcc <- bind_rows(dkcc,
                      stromaid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = StromaID))
    t1_meta <- bind_rows(t1_meta, stromaid[[]] %>% rownames_to_column("cell"))
  }

  # Step 4: Ureteric epithelium classification
  if (nrow(seurat[[]] %>% filter(LineageID == "UrEp")) > 2) {
    message("Running ureteric epithelium classification...")
    urepid <- suppressWarnings(scPred::scPredict(seurat[, seurat$LineageID == "UrEp"], reference = model2.urep,
                                 threshold = 0.0, max.iter.harmony = max.iter))
    urepid$scpred_prediction <- gsub("^Tip", "UTip", urepid$scpred_prediction)
    colnames(urepid[[]]) <- gsub("^Tip", "UTip", colnames(urepid[[]]))
    urepid$UrEpID <- urepid$scpred_prediction

    dkcc <- bind_rows(dkcc,
                      urepid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = UrEpID))
    t1_meta <- bind_rows(t1_meta, urepid[[]] %>% rownames_to_column("cell"))
  }

  # Finalize identity assignments
  if ("NephronID" %in% colnames(t1_meta)) {
    t1_meta$NephronID <- factor(t1_meta$NephronID, levels = c("EN", "DN", "PN", "RC"))
  }

  levels <- c("unassigned", "NPC", "NPC-like", "Endo",
              "EN", "DN", "PN", "RC",
              "SPC", "CS", "MS", "MesS", "Stroma_NC",
              "UTip", "UOS", "UIS",
              "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC")

  # Join DKCC assignments by cell name (order-safe join replaces old direct assignment)
  t1_meta <- left_join(t1_meta, dkcc, by = "cell")
  t1_meta$DKCC <- factor(t1_meta$dkcc, levels = levels[levels %in% c(unique(dkcc$dkcc), "NPC-like")])
  t1_meta$dkcc <- NULL

  t1_meta$LineageID <- factor(t1_meta$LineageID, levels = c("unassigned", "Endo", "Stroma", "NPC", "NPC-like", "Nephron", "UrEp"))

  # Transfer results back to original seurat object
  seurat[[]] <- left_join(md %>% rownames_to_column("cell") %>% select(cell),
                          t1_meta, by = "cell") %>% column_to_rownames("cell")

  # NPC refinement
  npc_count <- sum(seurat$LineageID == "NPC", na.rm = TRUE)

  if (npc_count > 2) {
    message("Refining NPC classification (", npc_count, " NPC cells)...")
    npcs <- seurat[, which(seurat$LineageID == "NPC")]
    DefaultAssay(npcs) <- "RNA"
    npcs <- npcs %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
    npcs$RNA_snn_res.0.5 <- 0

    if (ncol(npcs) > 50) {
      npcs <- npcs %>%
        RunPCA(npcs = 20) %>% RunUMAP(dims = 1:20) %>% FindNeighbors() %>% FindClusters(resolution = 0.5)
    }

    # 'orig.ident' does not survive the metadata rebuild above for objects that
    # arrived as h5ad (the per-lineage scPredict() outputs are recombined via
    # bind_rows()/left_join(), which drops it). Its only use is to populate
    # 'Identity', which nothing downstream reads -- the very next line replaces
    # orig.ident with "all" for GeneSummary()'s grouping -- so guard rather than
    # require it.
    npcs$Identity <- if ("orig.ident" %in% colnames(npcs[[]])) npcs$orig.ident else NA_character_
    npcs$orig.ident <- "all"

    if ("PAX2" %in% rownames(npcs)) {
      markers <- DevKidCC::GeneSummary(npcs, identity = "orig.ident", split.by = "RNA_snn_res.0.5", features = c("PAX2"))
      pax2null <- (markers %>% filter(pct.exp < 33))$Component
    } else {
      # PAX2 was dropped by the zero-variance filter above -- common when a small
      # chunk of NPC cells expresses no PAX2 at all. Every cluster is then
      # PAX2-null by definition, which is what the else branch below concludes
      # anyway; treating them all as NPC-like keeps that behaviour instead of
      # erroring out inside GeneSummary().
      message("PAX2 not present in NPC subset (zero variance in this chunk) - treating all NPC clusters as NPC-like")
      pax2null <- levels(factor(as.character(npcs$RNA_snn_res.0.5)))
    }

    if (length(pax2null) > 0) {
      names <- colnames(npcs[, npcs$RNA_snn_res.0.5 %in% pax2null])
      seurat[[]] <- within(seurat[[]], LineageID[LineageID %in% c('NPC') & rownames(seurat[[]]) %in% names] <- 'NPC-like')
      seurat[[]] <- within(seurat[[]], DKCC[DKCC %in% c('NPC') & rownames(seurat[[]]) %in% names] <- 'NPC-like')
    } else {
      seurat[[]] <- within(seurat[[]], DKCC[DKCC %in% c('NPC')] <- 'NPC-like')
    }

    seurat@misc$NPC.seu <- npcs
  } else if (npc_count > 0) {
    message("Too few NPC cells (", npc_count, ") for refinement, labeling as NPC-like")
    seurat[[]] <- within(seurat[[]], DKCC[DKCC %in% c('NPC')] <- 'NPC-like')
  } else {
    message("No NPC cells found, skipping NPC refinement")
  }

  return(seurat)
}
## KNN smoothing function to classify false negatives. Requires dimensional reduction, UMAP works best in testing.
fill_unassigned_by_knn_seurat <- function(
    obj,
    identity_col,
    reduction = "umap",
    k = 15,
    threshold = 0.7,
    max_iter = 10,
    nn_package = c("RANN", "FNN"),
    verbose = TRUE
) {
  nn_package <- match.arg(nn_package)
  
  # --- checks ---
  if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object.")
  if (!identity_col %in% colnames(obj@meta.data)) {
    stop(paste("Metadata column", identity_col, "not found."))
  }
  if (!reduction %in% Reductions(obj)) {
    stop(paste("Reduction", reduction, "not found in the object."))
  }
  
  # --- extract embeddings + metadata ---
  coords <- Embeddings(obj, reduction)
  labels <- as.character(obj[[identity_col]][,1])
  
  # convert to matrix if needed
  coords <- as.matrix(coords)
  
  # --- compute neighbors (includes self) ---
  if (nn_package == "RANN") {
    nn <- RANN::nn2(coords, coords, k = k)
    idx <- nn$nn.idx
  } else {
    nn <- FNN::get.knnx(coords, coords, k = k)
    idx <- nn$nn.index
  }
  
  if (verbose)
    message(sprintf("Starting KNN fill for 'unassigned' (threshold %.0f%%)", threshold * 100))
  
  # --- iterative fill ---
  current <- labels
  
  for (iter in seq_len(max_iter)) {
    unassigned <- which(current == "unassigned")
    if (length(unassigned) == 0) {
      if (verbose) message("Done: no unassigned left.")
      break
    }
    
    updated <- current
    filled <- 0
    
    for (cell in unassigned) {
      neighbor_labels <- current[idx[cell, ]]
      tab <- table(neighbor_labels)
      tab <- tab[names(tab) != "unassigned"]   # remove unassigned candidates
      if (length(tab) == 0) next
      
      best_label <- names(tab)[which.max(tab)]
      ratio <- max(tab) / k
      
      if (ratio >= threshold) {
        updated[cell] <- best_label
        filled <- filled + 1
      }
    }
    
    current <- updated
    if (verbose)
      message(sprintf("Iteration %d: filled %d cells", iter, filled))
    
    if (filled == 0) {
      if (verbose) message("Reached stability — stopping.")
      break
    }
  }
  
  # --- write the updated labels back into Seurat ---
  obj[[identity_col]] <- current
  
  return(obj)
}


