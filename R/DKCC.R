#' DKCC
#'
#' @param seurat seurat object
#' @param threshold minimum value for an identity to be assigned within the model call, default is 0.7
#' @param max.iter Can ask scPred to run this number of integrations, set to 0 be default
#'
#' @return seurat object with additional metadata columns
#' @export
#'
#' @aliases DKCC
#'
#' @examples
#' organoid <- DKCC(organoid)

DKCC <- function(seurat, threshold = 0.7, max.iter = 1) {

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

  # Step 1: Lineage classification
  message("Running lineage classification...")
  seurat <- scPred::scPredict(seurat, reference = model1.all, threshold = threshold, max.iter.harmony = max.iter)
  seurat$scpred_prediction <- gsub("Endothelial", "Endo", seurat$scpred_prediction)
  colnames(seurat[[]]) <- gsub("Endothelial", "Endo", colnames(seurat[[]]))
  seurat$LineageID <- seurat$scpred_prediction
  seurat$LineageID_max <- seurat$scpred_max
  fill_unassigned_by_knn_seurat(seurat, "LineageID", threshold = 0.5, k=15)

  dkcc <- seurat[[]] %>% rownames_to_column("cell") %>% filter(LineageID %in% c("unassigned", "NPC", "Endo")) %>% transmute(cell = cell, dkcc = LineageID)

  t1 <- seurat[, seurat$LineageID %in% c("unassigned", "NPC", "Endo")]

  # Step 2: Nephron classification
  if (nrow(seurat[[]] %>% filter(LineageID == "Nephron")) > 2) {
    message("Running nephron classification...")
    nephronid <- scPred::scPredict(seurat[, seurat$LineageID == "Nephron"], reference = model2.nephron,
                                    threshold = 0.0, max.iter.harmony = max.iter)
    nephronid$NephronID <- nephronid$scpred_prediction

    dkcc <- bind_rows(dkcc,
                      nephronid[[]] %>% rownames_to_column("cell") %>% filter(NephronID %in% c("EN")) %>% transmute(cell = cell, dkcc = NephronID))

    if (nrow(nephronid[[]] %>% filter(NephronID %in% c("EN"))) > 0) {
      neph <- nephronid[, nephronid$NephronID %in% c("EN")]
      t1 <- merge(t1, neph)
    }

    # Step 2a: Proximal nephron
    if (nrow(nephronid[[]] %>% filter(NephronID == "PN")) > 2) {
      message("Running proximal nephron classification...")
      proximalid <- scPred::scPredict(nephronid[, nephronid$NephronID == "PN"],
                                       reference = model3.pn, threshold = 0.0, max.iter.harmony = max.iter)
      proximalid$SegmentID <- proximalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        proximalid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, proximalid)
    }

    # Step 2b: Distal nephron
    if (nrow(nephronid[[]] %>% filter(NephronID == "DN")) > 2) {
      message("Running distal nephron classification...")
      distalid <- scPred::scPredict(nephronid[, nephronid$NephronID == "DN"],
                                     reference = model3.dn, threshold = 0.0, max.iter.harmony = max.iter)
      distalid$SegmentID <- distalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        distalid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, distalid)
    }

    # Step 2c: Renal corpuscle
    if (nrow(nephronid[[]] %>% filter(NephronID == "RC")) > 2) {
      message("Running renal corpuscle classification...")
      rcid <- scPred::scPredict(nephronid[, nephronid$NephronID == "RC"], reference = model3.rc,
                                 threshold = 0.0, max.iter.harmony = max.iter)
      rcid$SegmentID <- rcid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        rcid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, rcid)
    }
  }

  # Step 3: Stroma classification
  if (nrow(seurat[[]] %>% filter(LineageID == "Stroma")) > 2) {
    message("Running stroma classification...")
    stromaid <- scPred::scPredict(seurat[, seurat$LineageID == "Stroma"], reference = model2.stroma,
                                   threshold = 0.0, max.iter.harmony = max.iter)
    stromaid$StromaID <- stromaid$scpred_prediction

    dkcc <- bind_rows(dkcc,
                      stromaid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = StromaID))
    t1 <- merge(t1, stromaid)
  }

  # Step 4: Ureteric epithelium classification
  if (nrow(seurat[[]] %>% filter(LineageID == "UrEp")) > 2) {
    message("Running ureteric epithelium classification...")
    urepid <- scPred::scPredict(seurat[, seurat$LineageID == "UrEp"], reference = model2.urep,
                                 threshold = 0.0, max.iter.harmony = max.iter)
    urepid$scpred_prediction <- gsub("^Tip", "UTip", urepid$scpred_prediction)
    colnames(urepid[[]]) <- gsub("^Tip", "UTip", colnames(urepid[[]]))
    urepid$UrEpID <- urepid$scpred_prediction

    dkcc <- bind_rows(dkcc,
                      urepid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = UrEpID))
    t1 <- merge(t1, urepid)
  }

  # Finalize identity assignments
  if ("NephronID" %in% colnames(t1[[]])) {
    t1$NephronID <- factor(t1$NephronID, levels = c("EN", "DN", "PN", "RC"))
  }

  levels <- c("unassigned", "NPC", "NPC-like", "Endo",
              "EN", "DN", "PN", "RC",
              "SPC", "CS", "MS", "MesS", "Stroma_NC",
              "UTip", "UOS", "UIS",
              "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC")

  t1$DKCC <- factor(dkcc$dkcc, levels = levels[levels %in% c(unique(dkcc$dkcc), "NPC-like")])
  t1$LineageID <- factor(t1$LineageID, levels = c("unassigned", "Endo", "Stroma", "NPC", "NPC-like", "Nephron", "UrEp"))

  # Transfer results back to original seurat object
  seurat[[]] <- left_join(md %>% rownames_to_column("cell") %>% select(cell),
                          t1[[]] %>% rownames_to_column("cell"), by = "cell") %>% column_to_rownames("cell")

  # NPC refinement - SIMPLIFIED to avoid GeneSummary compatibility issues
  # All NPC cells are labeled as NPC-like
  npc_count <- sum(seurat$LineageID == "NPC", na.rm = TRUE)

  if (npc_count > 0) {
    message("Skipping NPC refinement (GeneSummary disabled), marking ", npc_count, " NPC cells as NPC-like")
    seurat[[]] <- within(seurat[[]], {
      LineageID[LineageID %in% c('NPC')] <- 'NPC-like'
      DKCC[DKCC %in% c('NPC')] <- 'NPC-like'
    })
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
    max_iter = 5,
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


