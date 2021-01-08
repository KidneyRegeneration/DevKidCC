#' DevKidCC
#'
#' @param seurat seurat object
#' @param threshold minimum value for an identity to be assigned within the model call, default is 0.7
#' @param max.iter Can ask scPred to run this number of integrations, set to 0 be default
#'
#' @return
#' @export
#'
#' @examples
DevKidCC <- function(seurat, threshold = 0.7, max.iter = 0) {
  if (("DKCC" %in% colnames(seurat@meta.data)) == FALSE){
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
  kcc <- data.frame()

  seurat <- scPred::scPredict(seurat, reference = model1.all, threshold = threshold, max.iter.harmony=max.iter)
  seurat$LineageID <- seurat$scpred_prediction

  kcc <- seurat@meta.data %>% rownames_to_column("cell") %>% filter(LineageID %in% c("unassigned", "NPC", "Endothelial")) %>% transmute(cell = cell, kcc = LineageID)
  #kcc$kcc <- gsub("unassigned", "OffTarget", kcc$kcc)
  t1 <- seurat[, seurat$LineageID%in%c("unassigned", "NPC", "Endothelial")]


  #seurat@meta.data <- seurat@meta.data[, 1:ncol(md)]

  if (nrow(seurat@meta.data %>% filter(LineageID=="Nephron")) > 2){
    nephronid <- scPred::scPredict(seurat[, seurat$LineageID=="Nephron"], reference = model2.nephron,
                                   threshold = 0.0, max.iter.harmony=max.iter)
    nephronid$NephronID <- nephronid$scpred_prediction
    nephronid$NephronID <- gsub("unassigned", "Nephron_NC", nephronid$NephronID)
    kcc <- bind_rows(kcc,
                     nephronid@meta.data %>% rownames_to_column("cell") %>% filter(NephronID %in% c("CellCycle", "EarlyNephron", "unassigned")) %>% transmute(cell = cell, kcc = NephronID))

    if (nrow(nephronid@meta.data %>% filter(NephronID %in% c("CellCycle", "EarlyNephron", "unassigned"))) > 0){
      neph <- nephronid[, nephronid$NephronID %in% c("CellCycle", "EarlyNephron", "unassigned")]
      t1 <- merge(t1, neph)
    }


    if (nrow(nephronid@meta.data %>% filter(NephronID=="ProximalNephron")) > 2){
      proximalid <- scPred::scPredict(nephronid[, nephronid$NephronID =="ProximalNephron"],
                                      reference = model3.pn, threshold = 0.0, max.iter.harmony=max.iter)
      proximalid$SegmentID <- proximalid$scpred_prediction

      kcc <- bind_rows(kcc,
                       proximalid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, kcc = SegmentID))
      t1 <- merge(t1, proximalid)
    }


    if (nrow(nephronid@meta.data %>% filter(NephronID=="DistalNephron")) > 2){
      distalid <- scPred::scPredict(nephronid[, nephronid$NephronID =="DistalNephron"],
                                    reference = model3.dn, threshold = 0.0, max.iter.harmony=max.iter)
      distalid$SegmentID <- distalid$scpred_prediction

      kcc <- bind_rows(kcc,
                       distalid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, kcc = SegmentID))
      t1 <- merge(t1, distalid)
    }

    if (nrow(nephronid@meta.data %>% filter(NephronID=="RenalCorpuscle")) > 2){
      rcid <- scPred::scPredict(nephronid[, nephronid$NephronID =="RenalCorpuscle"], reference = model3.rc,
                                threshold = 0.0, max.iter.harmony=max.iter)
      rcid$SegmentID <- rcid$scpred_prediction

      kcc <- bind_rows(kcc,
                       rcid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, kcc = SegmentID))
      t1 <- merge(t1, rcid)
    }
  }

  if (nrow(seurat@meta.data %>% filter(LineageID=="Stroma")) > 2){
    stromaid <- scPred::scPredict(seurat[, seurat$LineageID=="Stroma"], reference = model2.stroma,
                                  threshold = 0.0, max.iter.harmony=max.iter)
    stromaid$StromaID<- stromaid$scpred_prediction

    stromaid$StromaID <- gsub(pattern = "unassigned", replacement = "Stroma_NC", x = stromaid$StromaID)
    kcc <- bind_rows(kcc,
                     stromaid@meta.data %>% rownames_to_column("cell") %>% transmute(cell = cell, kcc = StromaID))
    t1 <- merge(t1, stromaid)
  }

  if (nrow(seurat@meta.data %>% filter(LineageID=="UrEp")) > 2){
    urepid <- scPred::scPredict(seurat[, seurat$LineageID=="UrEp"], reference = model2.urep,
                                threshold = 0.0, max.iter.harmony = max.iter)
    urepid$UrEpID <- urepid$scpred_prediction
    urepid$UrEpID <- gsub(pattern = "unassigned", replacement = "UrEp", x = urepid$UrEpID)
    kcc <- bind_rows(kcc,
                     urepid@meta.data %>% rownames_to_column("cell") %>%  transmute(cell = cell, kcc = UrEpID))
    t1 <- merge(t1, urepid)
  }

  if ("NephronID" %in% colnames(t1@meta.data)){
  t1$NephronID <- gsub("CellCycle", "CC", t1$NephronID)
  t1$NephronID <- gsub("EarlyNephron", "EN", t1$NephronID)
  t1$NephronID <- gsub("DistalNephron", "DN", t1$NephronID)
  t1$NephronID <- gsub("ProximalNephron", "PN", t1$NephronID)
  t1$NephronID <- gsub("RenalCorpuscle", "RC", t1$NephronID)
  t1$NephronID <- factor(t1$NephronID, levels = c("unassigned", "CC", "EN", "DN", "PN", "RC"))
  }

  kcc$kcc <- gsub("CellCycle", "CC", kcc$kcc)
  kcc$kcc <- gsub("EarlyNephron", "EN", kcc$kcc)
  levels <- c("unassigned", "NPC", "Endothelial",
             "CC", "EN", "DN", "PN", "RC",
             "SPC", "Cortex", "Medullary", "Mesangial", "Stroma_NC",
             "Tip", "OuterStalk", "InnerStalk",
             "EDT_EMT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC")

  t1$DKCC <- factor(kcc$kcc, levels = levels[levels %in% unique(kcc$kcc)])
  t1$LineageID <- factor(t1$LineageID, levels = c("unassigned", "Endothelial", "Stroma", "NPC", "Nephron", "UrEp"))



  #


  old.seurat@meta.data <- left_join(md %>% rownames_to_column("cell") %>% select(cell),
                                    t1@meta.data %>% rownames_to_column("cell"), by = "cell") %>% column_to_rownames("cell")

  return(old.seurat)
}










