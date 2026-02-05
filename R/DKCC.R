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
  #if ("cell" %in% colnames(md)) {
  #  seurat[[]] %>% rownames_to_column("celltemp") %>% select(-cell) %>% column_to_rownames("celltemp")
  #}
  old.seurat <- seurat

  # Handle Seurat v5 multi-layer objects
  # Check if we have an Assay5 with multiple layers
  DefaultAssay(old.seurat) <- "RNA"
  if (inherits(old.seurat[["RNA"]], "Assay5")) {
    # Get layer names
    layers <- SeuratObject::Layers(old.seurat, search = "data")
    if (length(layers) > 1) {
      # Multiple layers detected - join them
      message("Detected ", length(layers), " data layers. Joining layers for processing...")
      old.seurat[["RNA"]] <- SeuratObject::JoinLayers(old.seurat[["RNA"]])
    }
  }

  # Now safely extract data
  data_matrix <- get_layer_data(old.seurat, layer = "data", assay = "RNA")

  # Workaround for scPred compatibility: Force v4-style Assay creation
  # scPred hasn't been updated for Seurat v5 yet, so we need a compatibility layer
  old_option <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")  # Force v4-compatible assay

  seurat <- CreateSeuratObject(data_matrix, meta.data = md)
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat)

  # Restore original option
  if (is.null(old_option)) {
    options(Seurat.object.assay.version = NULL)
  } else {
    options(Seurat.object.assay.version = old_option)
  }

  # Additional workaround: Use trace() to intercept GetAssayData calls from scPred
  # scPred calls GetAssayData(object, "data") with old v4 syntax
  # We use trace() to convert the parameters before the actual function is called

  trace(
    SeuratObject:::GetAssayData.Seurat,
    tracer = quote({
      # Convert old calling conventions to Seurat v5 syntax
      # If second arg looks like a layer name, move it to the layer parameter
      if (!missing(assay) && is.character(assay) && assay %in% c("data", "counts", "scale.data")) {
        if (missing(layer) || is.null(layer)) {
          layer <- assay
          assay <- NULL
        }
      }
      # Handle deprecated slot parameter
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

  # Ensure trace is removed even if an error occurs
  on.exit({
    try(untrace(SeuratObject:::GetAssayData.Seurat, where = asNamespace("SeuratObject")), silent = TRUE)
  }, add = TRUE)

  dkcc <- data.frame()

  # Now scPredict should work with the patched GetAssayData
  seurat <- scPred::scPredict(seurat, reference = model1.all, threshold = threshold, max.iter.harmony=max.iter)
  seurat$scpred_prediction <- gsub("Endothelial", "Endo", seurat$scpred_prediction)
  colnames(seurat[[]]) <- gsub("Endothelial", "Endo", colnames(seurat[[]]))
  seurat$LineageID <- seurat$scpred_prediction
  seurat$LineageID_max <- seurat$scpred_max


  dkcc <- seurat[[]] %>% rownames_to_column("cell") %>% filter(LineageID %in% c("unassigned", "NPC", "Endo")) %>% transmute(cell = cell, dkcc = LineageID)

  t1 <- seurat[, seurat$LineageID%in%c("unassigned", "NPC", "Endo")]


  #seurat[[]] <- seurat[[]][, 1:ncol(md)]

  if (nrow(seurat[[]] %>% filter(LineageID=="Nephron")) > 2){
    nephronid <- scPred::scPredict(seurat[, seurat$LineageID=="Nephron"], reference = model2.nephron,
                                   threshold = 0.0, max.iter.harmony=max.iter)
    nephronid$NephronID <- nephronid$scpred_prediction
    #nephronid$NephronID <- gsub("unassigned", "Nephron_NC", nephronid$NephronID)
    dkcc <- bind_rows(dkcc,
                      nephronid[[]] %>% rownames_to_column("cell") %>% filter(NephronID %in% c("EN")) %>% transmute(cell = cell, dkcc = NephronID))

    if (nrow(nephronid[[]] %>% filter(NephronID %in% c("EN"))) > 0){
      neph <- nephronid[, nephronid$NephronID %in% c("EN")]
      t1 <- merge(t1, neph)
    }


    if (nrow(nephronid[[]] %>% filter(NephronID=="PN")) > 2){
      proximalid <- scPred::scPredict(nephronid[, nephronid$NephronID =="PN"],
                                      reference = model3.pn, threshold = 0.0, max.iter.harmony=max.iter)
      proximalid$SegmentID <- proximalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        proximalid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, proximalid)
    }


    if (nrow(nephronid[[]] %>% filter(NephronID=="DN")) > 2){
      distalid <- scPred::scPredict(nephronid[, nephronid$NephronID =="DN"],
                                    reference = model3.dn, threshold = 0.0, max.iter.harmony=max.iter)
      distalid$SegmentID <- distalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        distalid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, distalid)
    }

    if (nrow(nephronid[[]] %>% filter(NephronID=="RC")) > 2){
      rcid <- scPred::scPredict(nephronid[, nephronid$NephronID =="RC"], reference = model3.rc,
                                threshold = 0.0, max.iter.harmony=max.iter)
      rcid$SegmentID <- rcid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        rcid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, rcid)
    }
  }

  if (nrow(seurat[[]] %>% filter(LineageID=="Stroma")) > 2){
    stromaid <- scPred::scPredict(seurat[, seurat$LineageID=="Stroma"], reference = model2.stroma,
                                  threshold = 0.0, max.iter.harmony=max.iter)
    stromaid$StromaID<- stromaid$scpred_prediction

    #stromaid$StromaID <- gsub(pattern = "unassigned", replacement = "Stroma_NC", x = stromaid$StromaID)
    dkcc <- bind_rows(dkcc,
                      stromaid[[]] %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = StromaID))
    t1 <- merge(t1, stromaid)
  }

  if (nrow(seurat[[]] %>% filter(LineageID=="UrEp")) > 2){
    urepid <- scPred::scPredict(seurat[, seurat$LineageID=="UrEp"], reference = model2.urep,
                                threshold = 0.0, max.iter.harmony = max.iter)
    urepid$scpred_prediction <- gsub("^Tip", "UTip", urepid$scpred_prediction)
    colnames(urepid[[]]) <- gsub("^Tip", "UTip", colnames(urepid[[]]))
    urepid$UrEpID <- urepid$scpred_prediction
    #urepid$UrEpID <- gsub(pattern = "unassigned", replacement = "UrEp", x = urepid$UrEpID)
    dkcc <- bind_rows(dkcc,
                      urepid[[]] %>% rownames_to_column("cell") %>%  transmute(cell = cell, dkcc = UrEpID))
    t1 <- merge(t1, urepid)
  }
  # fix up some of the identity names
  if ("NephronID" %in% colnames(t1[[]])){

    t1$NephronID <- factor(t1$NephronID, levels = c("EN", "DN", "PN", "RC"))
  }



  levels <- c("unassigned", "NPC", "NPC-like", "Endo",
              "EN", "DN", "PN", "RC",
              "SPC", "CS", "MS", "MesS", "Stroma_NC",
              "UTip", "UOS", "UIS",
              "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC")

  t1$DKCC <- factor(dkcc$dkcc, levels = levels[levels %in% c(unique(dkcc$dkcc), "NPC-like")])
  t1$LineageID <- factor(t1$LineageID, levels = c("unassigned", "Endo", "Stroma", "NPC", "NPC-like", "Nephron", "UrEp"))



  #


  old.seurat[[]] <- left_join(md %>% rownames_to_column("cell") %>% select(cell),
                                    t1[[]] %>% rownames_to_column("cell"), by = "cell") %>% column_to_rownames("cell")


  # fix the NPC problem - only if there are NPC cells

  # Check if there are any NPC cells before subsetting
  npc_count <- sum(old.seurat$LineageID == "NPC", na.rm = TRUE)

  if (npc_count > 0) {
    npcs <- old.seurat[, old.seurat$LineageID == "NPC"]
    DefaultAssay(npcs) <- "RNA"
    npcs <- npcs %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
    npcs$RNA_snn_res.0.5 <- 0
    if (ncol(npcs)>50) {
      npcs <- npcs %>%
      RunPCA(npcs = 20) %>% RunUMAP(dims = 1:20) %>% FindNeighbors() %>% FindClusters(resolution = 0.5)
    } else {

    }

    npcs$Identity <- npcs$orig.ident
    npcs$orig.ident <- "all"
    markers <- DevKidCC::GeneSummary(npcs, identity = "orig.ident", split.by = "RNA_snn_res.0.5", features = c("PAX2"))
    pax2null <- (markers %>% filter(pct.exp<33))$Component
    if (length(pax2null)>0){
    names <- colnames(npcs[, npcs$RNA_snn_res.0.5 %in% pax2null])
    old.seurat[[]] <- within(old.seurat[[]], LineageID[LineageID %in% c('NPC') & rownames(old.seurat[[]]) %in% names] <- 'NPC-like')
    old.seurat[[]] <- within(old.seurat[[]], DKCC[DKCC %in% c('NPC') & rownames(old.seurat[[]]) %in% names] <- 'NPC-like')
    } else {
      old.seurat[[]] <- within(old.seurat[[]], DKCC[DKCC %in% c('NPC')] <- 'NPC-like')
    }

    old.seurat@misc$NPC.seu <- npcs
  } else {
    message("No NPC cells found, skipping NPC refinement")
  }

  return(old.seurat)
}


#UpdateThreshold <- function(seurat, threshold) {
#  md <- seurat[[]]
#  thresh <- md$LineageID_max < threshold
#  md[thresh,]$LineageID <- "unassigned"
#  md[thresh,]$DKCC <- "unassigned"
#  seurat@md <- md
#  return(seurat)
#}
#
