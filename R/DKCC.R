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


  if (("dkcc" %in% colnames(seurat@meta.data)) == FALSE){
    seurat@misc$old.meta <- seurat@meta.data
  } else {
    seurat@meta.data <- seurat@misc$old.meta
  }
  md <- seurat@meta.data
  #if ("cell" %in% colnames(md)) {
  #  seurat@meta.data %>% rownames_to_column("celltemp") %>% select(-cell) %>% column_to_rownames("celltemp")
  #}
  old.seurat <- seurat
  seurat <- CreateSeuratObject(old.seurat@assays$RNA@counts, meta.data = md)
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat)
  dkcc <- data.frame()

  seurat <- scPred::scPredict(seurat, reference = model1.all, threshold = threshold, max.iter.harmony=max.iter) # changed threshold to 0, all cells assigned
  seurat$scpred_prediction <- gsub("Endothelial", "Endo", seurat$scpred_prediction)
  colnames(seurat@meta.data) <- gsub("Endothelial", "Endo", colnames(seurat@meta.data))
  seurat$LineageID <- seurat$scpred_prediction
  seurat$LineageID_max <- seurat$scpred_max


  dkcc <- seurat@meta.data %>% rownames_to_column("cell") %>% filter(LineageID %in% c("unassigned", "NPC", "Endo")) %>% transmute(cell = cell, dkcc = LineageID)

  t1 <- seurat[, seurat$LineageID%in%c("unassigned", "NPC", "Endo")]


  #seurat@meta.data <- seurat@meta.data[, 1:ncol(md)]

  if (nrow(seurat@meta.data %>% filter(LineageID=="Nephron")) > 2){
    nephronid <- scPred::scPredict(seurat[, seurat$LineageID=="Nephron"], reference = model2.nephron,
                                   threshold = 0.0, max.iter.harmony=max.iter)
    nephronid$NephronID <- nephronid$scpred_prediction
    #nephronid$NephronID <- gsub("unassigned", "Nephron_NC", nephronid$NephronID)
    dkcc <- bind_rows(dkcc,
                      nephronid@meta.data %>% rownames_to_column("cell") %>% filter(NephronID %in% c("EN")) %>% transmute(cell = cell, dkcc = NephronID))

    if (nrow(nephronid@meta.data %>% filter(NephronID %in% c("EN"))) > 0){
      neph <- nephronid[, nephronid$NephronID %in% c("EN")]
      t1 <- merge(t1, neph)
    }


    if (nrow(nephronid@meta.data %>% filter(NephronID=="PN")) > 2){
      proximalid <- scPred::scPredict(nephronid[, nephronid$NephronID =="PN"],
                                      reference = model3.pn, threshold = 0.0, max.iter.harmony=max.iter)
      proximalid$SegmentID <- proximalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        proximalid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, proximalid)
    }


    if (nrow(nephronid@meta.data %>% filter(NephronID=="DN")) > 2){
      distalid <- scPred::scPredict(nephronid[, nephronid$NephronID =="DN"],
                                    reference = model3.dn, threshold = 0.0, max.iter.harmony=max.iter)
      distalid$SegmentID <- distalid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        distalid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, distalid)
    }

    if (nrow(nephronid@meta.data %>% filter(NephronID=="RC")) > 2){
      rcid <- scPred::scPredict(nephronid[, nephronid$NephronID =="RC"], reference = model3.rc,
                                threshold = 0.0, max.iter.harmony=max.iter)
      rcid$SegmentID <- rcid$scpred_prediction

      dkcc <- bind_rows(dkcc,
                        rcid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = SegmentID))
      t1 <- merge(t1, rcid)
    }
  }

  if (nrow(seurat@meta.data %>% filter(LineageID=="Stroma")) > 2){
    stromaid <- scPred::scPredict(seurat[, seurat$LineageID=="Stroma"], reference = model2.stroma,
                                  threshold = 0.0, max.iter.harmony=max.iter)
    stromaid$StromaID<- stromaid$scpred_prediction

    #stromaid$StromaID <- gsub(pattern = "unassigned", replacement = "Stroma_NC", x = stromaid$StromaID)
    dkcc <- bind_rows(dkcc,
                      stromaid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, dkcc = StromaID))
    t1 <- merge(t1, stromaid)
  }

  if (nrow(seurat@meta.data %>% filter(LineageID=="UrEp")) > 2){
    urepid <- scPred::scPredict(seurat[, seurat$LineageID=="UrEp"], reference = model2.urep,
                                threshold = 0.0, max.iter.harmony = max.iter)
    urepid$scpred_prediction <- gsub("^Tip", "UTip", urepid$scpred_prediction)
    colnames(urepid@meta.data) <- gsub("^Tip", "UTip", colnames(urepid@meta.data))
    urepid$UrEpID <- urepid$scpred_prediction
    #urepid$UrEpID <- gsub(pattern = "unassigned", replacement = "UrEp", x = urepid$UrEpID)
    dkcc <- bind_rows(dkcc,
                      urepid@meta.data %>% rownames_to_column("cell") %>%  transmute(cell = cell, dkcc = UrEpID))
    t1 <- merge(t1, urepid)
  }
  # fix up some of the identity names
  if ("NephronID" %in% colnames(t1@meta.data)){

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


  old.seurat@meta.data <- left_join(md %>% rownames_to_column("cell") %>% select(cell),
                                    t1@meta.data %>% rownames_to_column("cell"), by = "cell") %>% column_to_rownames("cell")


  # fix the NPC problem

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
  names <- colnames(npcs[, npcs$RNA_snn_res.0.5 %in% pax2null])
  old.seurat@meta.data <- within(old.seurat@meta.data, LineageID[LineageID %in% c('NPC') & rownames(old.seurat@meta.data) %in% names] <- 'NPC-like')
  old.seurat@meta.data <- within(old.seurat@meta.data, DKCC[DKCC %in% c('NPC') & rownames(old.seurat@meta.data) %in% names] <- 'NPC-like')

  old.seurat@misc$NPC.seu <- npcs


  return(old.seurat)
}


UpdateThreshold <- function(seurat, threshold) {
  md <- seurat@meta.data
  thresh <- md$LineageID_max < threshold
  md[thresh,]$LineageID <- "unassigned"
  md[thresh,]$DKCC <- "unassigned"
  seurat@md <- md
  return(seurat)
}

