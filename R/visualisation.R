#' ComponentPlot
#'
#' @param data
#' @param identity
#' @param component
#' @param feature
#' @param show.unassigned
#' @param show.pct
#' @param show.gene.exp
#' @param do.label
#'
#' @return
#' @export
#'
#' @examples
ComponentPlot <- function(data, identity = "orig.ident", component = "DKCC", feature = "MALAT1",
                          show.unassigned = TRUE, show.pct = FALSE, show.gene.exp = FALSE, do.label = T){
  data.plot <- GeneSummary(data, identity = identity, split.by = component, features = feature)

  myColors <- gplots::col2hex(c("grey", "grey",  "royalblue1", "brown", "darkgreen", "green", "skyblue","palevioletred4",
                                "peachpuff2", "goldenrod", "tan2", "wheat3",
                                "lightgreen", "palegreen4", "forestgreen", "goldenrod", "tan3", "lightskyblue3", "cyan", "royalblue3", "grey20",
                                "orchid4", "orchid1", "maroon2", "magenta", "mediumpurple2",
                                "orangered1", "wheat3", "goldenrod4"
  ))
  #Create a custom color scale

  names(myColors) <-  c("OffTarget", "unassigned", "Endo", "Stroma", "NPC-like", "NPC", "Nephron", "UrEp",
                        "EN", "DN", "PN", "RC",
                        "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC",
                        "SPC", "CS", "MS", "MesS", "Stroma_NC",
                        "UTip", "UOS", "UIS"
  )
  data.plot <- data.plot %>% filter(Cells>0)

  data.plot$Component <- factor(data.plot$Component, levels = names(myColors)[names(myColors) %in% unique(as.character(data.plot$Component))])
  if (!is.null(levels(data.plot$Identity))){
    data.plot$Identity <- factor(data.plot$Identity, levels = levels(data@meta.data[, identity]))
  }


  data.plot <- data.plot %>% arrange(desc(Component))
  if (show.unassigned == FALSE) {
    data.plot <- data.plot %>% filter(Component != "unassigned")
  }

  myColors <- myColors[levels(data.plot$Component)]

  fct.order <- levels(data@meta.data[[identity]])
  data.plot <- data.plot %>% mutate(Identity = fct_relevel(Identity, fct.order))

  #ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(Identity==.x) %>% select(Pct))))
  if (show.pct == FALSE){
    p <- ggplot(data.plot, aes(.data$Identity, .data$Cells))

    if (show.gene.exp == T) {
      title <- ggtitle(paste0(feature, " expression across ", component))
      b <- geom_bar(aes(fill = .data$avg.exp.log), stat = "Identity", colour = "black", width = 0.99)
      max <- max(data.plot %>% filter(features.plot == feature) %>% select(avg.exp.log))

      #c <- scale_fill_gradient2(low = "lightgrey", mid = "navy", high = "red", midpoint = max/2, na.value = "black")
      c <- scale_fill_gradient2(low = "lightgrey", mid = "red", high = "green", midpoint = max/2, na.value = "black")
      t <- geom_text(aes(label=ifelse((do.label==T & .data$avg.exp.log >= (max/3)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    } else {
      title <- ggtitle(paste0(component, " cell numbers"))
      b <- geom_bar(aes(fill = .data$Component), stat = "Identity", colour = "black", width = 0.99)
      #c <- scale_fill_manual(name = "Identity", values = (myColors[levels(data.plot$Component)]))
      c <- scale_fill_manual(name = "Identity", values = myColors)
      t <- geom_text(aes(label=ifelse((do.label==T & .data$Cells >= (ymax/10)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    }
    ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(Identity==.x) %>% select(Cells))))
  } else {
    p <- ggplot(data.plot, aes(.data$Identity, .data$Pct))
    if (show.gene.exp == T) {
      title <- ggtitle(paste0(feature, " expression across ", component))
      b <- geom_bar(aes(fill = .data$avg.exp.log), stat = "Identity", colour = "black", width = 0.99)
      max <- max(data.plot %>% filter(features.plot == feature) %>% select(avg.exp.log))

      #c <- scale_fill_gradient2(low = "lightgrey", mid = "navy", high = "red", midpoint = max/2, na.value = "black")
      c <- scale_fill_gradient2(low = "lightgrey", mid = "red", high = "green", midpoint = max/2, na.value = "black")
      t <- geom_text(aes(label=ifelse((do.label==T & .data$avg.exp.log >= (max/3)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    } else {
      title <- ggtitle(paste0(component, " cell proportions"))
      b <- geom_bar(aes(fill = .data$Component), stat = "Identity", colour = "black", width = 0.99)
      #c <- scale_fill_manual(name = "Identity", values = (myColors[levels(data.plot$Component)]))
      c <- scale_fill_manual(name = "Identity", values = myColors)
      t <- geom_text(aes(label=ifelse((do.label==T & .data$Pct >= (ymax/10)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    }
    ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(Identity==.x) %>% select(Pct))))
  }

  p +
    b +
    t +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    theme(legend.title=element_text(size=rel(1.1))) +
    scale_y_continuous(limits = c(0,ymax+1), expand = c(0, 0)) +
    c +
    title
}




#' SankeyPlot
#' @description Plots a Sankey diagram shown how the cells are from a category into the final classification
#'
#' @param data seurat object fot the data to be pulled from
#' @param origin the meta.data column describing the origin identity
#' @param height the height of the plot
#' @param width the width of the plot
#' @param fontsize the fontsize of the plot
#' @param simple defualt is FALSE; TRUE changes to a direct plot going from the origin to final classification without intermediate steps shown
#' @param forced default FALSE IF TRUE changes the outcome to show all cells force called into a cell type instead of being unassigned
#' @param sinksRight default FALSE; graphing parameter from base plotting function that forces final identity to right side of graph
#' @param LinkGroup Changes the way in which linking colours are plotted. Default is "Var3" representing the origin population
#' @param remove.unassigned default FALSE; determines if unassigned values are removed from the plot
#' @param nephron.only default FALSE; TRUE will plot only cell types that are classified into a nephron lineage
#'
#' @importFrom networkD3 sankeyNetwork
#'
#' @return the plot
#' @export
#'
#'
SankeyPlot <- function(data,
                       simple = T,
                       origin = "sample.id",
                       end = "DKCC",
                       height = 900,
                       width = 1200,
                       fontsize = 12,
                       sinksRight = TRUE,
                       LinkGroup = "Var1",
                       remove.unassigned = F,
                       nephron.only = F){

  if (remove.unassigned==T){
    .data$data <- data[, data$DKCC != "unassigned"]

  }

  if (nephron.only == T){
    .data$data[, data$DKCC %!in% c("Stroma", "unassigned", "Non_Renal")]
  }



  if (simple == T){
    links <- as.data.frame(table(data@meta.data[[origin]], data@meta.data[[end]]))
  } else {
    links <- as.data.frame(table(data@meta.data[[origin]], data@meta.data$LineageID)) %>% mutate(Var3 = .data$Var1)
    links <- bind_rows(links, (as.data.frame(table(data$LineageID, data$DKCC, data@meta.data[[origin]]),
                                             stringsAsFactors = F) %>% filter(Var1 != Var2))) %>% filter(.data$Freq!=0)

  }
  links <- filter(links, .data$Freq != 0)

  nodes <- data.frame(name = c(as.character(links$Var1), as.character(links$Var2)),
                      stringsAsFactors = F) %>% unique()

  links$IDsource <- match(links$Var1, nodes$name)-1
  links$IDtarget <- match(links$Var2, nodes$name)-1

  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "Freq", NodeID = "name",
                     sinksRight=sinksRight, fontSize = fontsize,
                     height = height, width = width, LinkGroup = LinkGroup)
  return(p)

}

MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

'%!in%' <- function(x,y) !('%in%'(x,y))

#LoadGE <- function(){
#  https://drive.google.com/file/d/1_6L0EKsYcHq2cS7cRK2NM5G6JCEZClW-/view?usp=sharing
#}

#' DotPlotCompare
#'
#' @param new.object seurat object to compare to reference data
#' @param show this defines the focus of the output; "gene" partitions by gene, anything else partitions by segment
#' @param split.by meta.data column to split the sample by, usually the sample column
#' @param idents Identities to display. Use "all" to show all segments, "nephron" to show all nephron segments, "others" to show all non-nephron segments, or a string of identities to select specific Identities
#' @param dot.min dot size of minimum cell percentage
#' @param dot.scale dot size of maximum cell percentage
#' @param scale.min minimum scale
#' @param scale.max maximum scale
#' @param scale.by "radius" or "size"
#' @param features features to plot
#' @param compare.to "reference" or "organoids" will be displayed on the plot as a comparison
#' @param scaling choose "log", "scale" or "raw" to change how the expression values are displayed
#' @param classification choose "DKCC" or "LineageID"
#' @param col.min minimum expression value to plot
#' @param col.max maximum expression value to plot
#'
#' @return
#' @export
#'


DotPlotCompare <- function(new.object = NULL,
                           show = "gene",
                           split.by = "orig.ident",
                           scaling = "log",
                           size = "pct.exp",
                           idents = "all",
                           classification = "DKCC",
                           dot.min = 0.1,
                           dot.scale = 6,
                           col.min = NA,
                           col.max = NA,
                           scale.min = NA,
                           scale.max = NA,
                           scale.by = "radius",
                           features,
                           compare.to.organoids = FALSE,
                           filter.samples = NULL,
                           columns = NULL) {


  comp.data <- data.frame()
  comp.data <- reference %>% filter(Tier == classification)


  if (compare.to.organoids == TRUE) {
    if ("data/GeneExpressionList.rda" %!in% list.files(path = here::here(), all.files = T, recursive = T)){
      warning("Please refer to https://github.com/KidneyRegeneration/DevKidCC for instructions on downloading the database")
    } else {
    if (!exists("GE")){
      load(file = here::here("data/GeneExpressionList.rda"))
    }
    if (classification == "DKCC"){
      comp.data <- bind_rows(comp.data, GE$dkcc)
    }
    if (classification == "LineageID"){
      comp.data <- bind_rows(comp.data, GE$lineageid)
    }
    if (classification == "NephronID"){
      comp.data <- bind_rows(comp.data, GE$nephronid)
    }
    if (classification %!in% c("DKCC", "LineageID", "NephronID")){
      warning("classification needs to be 'DKCC', 'LineageID' or 'NephronID'")
    }
    }
  }

  if (!is.null(filter.samples)){
    comp.data <- comp.data %>% filter(Identity %in% filter.samples)
  }

  samples <- comp.data$Identity %>% unique()

  #features <- factor(features, levels = features)
  if (classification == "DKCC"){
    identity <- factor(c("Nephron_NC", "Stroma_NC", "NPC-like",
                         unique(comp.data$Component)), levels = c("NPC", "NPC-like", "EN", "EDT_EMT", "DT", "LOH",
                                                           "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC",
                                                           "Tip", "OuterStalk", "InnerStalk", "Stroma_NC",
                                                           "SPC", "Cortex", "Medullary", "Mesangial", "Endothelial", "unassigned"))
  } else if (classification == "LineageID"){
    identity <- factor(c(unique(comp.data$Component), "unassigned", "NPC-like"), levels = c("Endothelial", "Stroma", "NPC", "NPC-like", "Nephron", "UrEp", "unassigned"))
  } else {
    identity <- factor(c(unique(comp.data$Component), "unassigned"), levels = c("EN", "DN", "PN", "RC", "unassigned"))
  }

  if(idents == "all"){
    identity <- identity
  } else if (idents == "nephron"){
    identity <- identity[identity %in% c("NPC", "NPC-like", "EN", "EDT_EMT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod")]
  } else if (idents == "others"){
    identity <- identity[identity %in% c("SPC", "Cortex", "Medullary", "Mesangial", "Endothelial", "unassigned")]
  } else if (idents %in% identity){
    identity <- identity[identity %in% idents]
  } else {
    warning("This call needs to specify identities eg. all, nephron, others, or specific identities")
  }
  if (!is.null(new.object)){
    id <- split.by
    new.data <- GeneSummary(new.object, features = features, split.by = classification, identity = id)
    data.plot <- bind_rows(new.data,
                           comp.data %>% filter(.data$features.plot %in% features))
    if (is.null(levels(new.object@meta.data[[split.by]]))){
      samples <- c(samples, as.character(unique(new.object@meta.data[[split.by]])))
    } else {
      samples <- c(samples, levels(new.object@meta.data[[split.by]]))
    }

  } else {
    data.plot <- comp.data %>% filter(.data$features.plot %in% features)
  }
  data.plot <-  data.plot %>%
    filter(.data$Component %in% identity, !is.na(features.plot)) %>%
    mutate(id = factor(.data$Component, levels = levels(identity),),
           features.plot = factor(.data$features.plot, levels = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  if (scaling == "log"){
    color.by <- "avg.exp.log"
  } else if (scaling == "scale"){
    color.by <- "avg.exp.scaled"
  } else if (scaling == "raw"){
    color.by <- "avg.exp"
  } else (
    warning("Choose 'log', 'scale' or 'raw' please")
  )

  data.plot[,color.by] <- MinMax(data = data.plot[, color.by], min = col.min,
                                 max = col.max)
  max <- max(data.plot[, color.by])
  scale.func <- switch(EXPR = scale.by, size = ggplot2::scale_size,
                       radius = ggplot2::scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

  if (show=="gene"){
    xaxis = "id"
    facet = "features.plot"
  } else {
    xaxis = "features.plot"
    facet = "id"
  }

  if (size == "pct.exp"){
    size.label <- "% Cells w/ GE"
  } else if (size == "Pct") {
    size.label <- "% Cells in Sample"
  } else {
    warning("Size needs to be 'Pct' or 'pct.exp'")
  }

  data.plot[, color.by][data.plot[, color.by]==0] <- NA
  data.plot$Identity <- factor(data.plot$Identity, levels = samples)
  plot <- ggplot(data = data.plot, mapping = aes_string(x = xaxis,
                                                        y = "Identity")) +
    geom_point(mapping = aes_string(size = size,
                                    color = color.by)) +

    scale.func(range = c(1.5, dot.scale),
               limits = c(scale.min, scale.max)) +
    cowplot::theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    #scale_colour_gradient(low = "navy", high = "red") +
    scale_colour_gradient2(low = "lightgrey", mid = "navy", high = "red", midpoint = max/2, na.value = "black") +
    #scale_colour_viridis_c(begin = 0.2, end = 0.8) +
    #guides(size = guide_legend(title = size.label)) +


    guides(size = guide_legend(title = size.label)) +

    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
          panel.grid.major = element_line(colour = "gray"),
          #panel.background = element_rect(colour = "lightgray")
    ) +
    theme(legend.title=element_text(size=rel(0.5))) +

    labs(x = "Features", y = "Identity") +
    facet_wrap(facets = facet, ncol = ifelse(is.null(columns), round(sqrt(length(unique(data.plot[, facet])))), columns))
  return(plot)
}




#' IdentBoxPlot
#'
#' @param data
#' @param identity
#' @param component
#' @param feature
#' @param show.unassigned
#' @param show.pct
#' @param show.gene.exp
#' @param do.label
#'
#' @return
#' @export
#'
#' @examples
IdentBoxPlot <- function(data, group, identity = "orig.ident", component = "DKCC", feature = "MALAT1", column = T,
                          show.unassigned = TRUE, show.pct = FALSE, show.gene.exp = FALSE, do.label = T, scales = "free_y"){
  data.list <- SplitObject(data, split.by = group)
  data.plot <- map2_dfr(data.list, names(data.list), ~GeneSummary(.x, identity = identity, split.by = component, features = feature) %>%
                          mutate(Group = factor(.y, levels = levels(data@meta.data[, group]))))

  myColors <- gplots::col2hex(c("grey", "grey",  "royalblue1", "brown", "darkgreen", "green", "skyblue","palevioletred4",
                                "peachpuff2", "goldenrod", "tan2", "wheat3",
                                "lightgreen", "palegreen4", "forestgreen", "goldenrod", "tan3", "lightskyblue3", "cyan", "royalblue3", "grey20",
                                "orchid4", "orchid1", "maroon2", "magenta", "mediumpurple2",
                                "orangered1", "wheat3", "goldenrod4"
  ))
  #Create a custom color scale

  names(myColors) <-  c("OffTarget", "unassigned", "Endo", "Stroma", "NPC-like", "NPC", "Nephron", "UrEp",
                        "EN", "DN", "PN", "RC",
                        "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC",
                        "SPC", "CS", "MS", "MesS", "Stroma_NC",
                        "UTip", "UOS", "UIS"
  )
  #data.plot <- data.plot %>% filter(Cells>0)

  data.plot$Component <- factor(data.plot$Component, levels = names(myColors)[names(myColors) %in% unique(as.character(data.plot$Component))])
  if (!is.null(levels(data.plot$Identity))){
    data.plot$Identity <- factor(data.plot$Identity, levels = levels(data@meta.data[, identity]))
  }


  data.plot <- data.plot %>% arrange(desc(Component))
  if (show.unassigned == FALSE) {
    data.plot <- data.plot %>% filter(Component != "unassigned")
  }

  myColors <- myColors[levels(data.plot$Component)]

  fct.order <- levels(data@meta.data[[identity]])
  data.plot <- data.plot %>% mutate(Identity = fct_relevel(Identity, fct.order))

  #ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(Identity==.x) %>% select(Pct))))
  #if (show.pct == FALSE){
  #  p <- ggplot(data.plot, aes(.data$Identity, .data$Cells))
  #
  #  if (show.gene.exp == T) {
  #    title <- ggtitle(paste0(feature, " expression across ", component))
  #    b <- geom_bar(aes(fill = .data$avg.exp.log), stat = "Identity", colour = "black", width = 0.99)
  #    max <- max(data.plot %>% filter(features.plot == feature) %>% select(avg.exp.log))
  #
  #    #c <- scale_fill_gradient2(low = "lightgrey", mid = "navy", high = "red", midpoint = max/2, na.value = "black")
  #    c <- scale_fill_gradient2(low = "lightgrey", mid = "red", high = "green", midpoint = max/2, na.value = "black")
  #    t <- geom_text(aes(label=ifelse((do.label==T & .data$avg.exp.log >= (max/3)), levels(.data$Component)[.data$Component], "")), size = 3,
  #                   position=position_stack(vjust=0.5), colour="black")
  #  } else {
  #    title <- ggtitle(paste0(component, " cell numbers"))
  #    b <- geom_bar(aes(fill = .data$Component), stat = "Identity", colour = "black", width = 0.99)
  #    #c <- scale_fill_manual(name = "Identity", values = (myColors[levels(data.plot$Component)]))
  #    c <- scale_fill_manual(name = "Identity", values = myColors)
  #    t <- geom_text(aes(label=ifelse((do.label==T & .data$Cells >= (ymax/10)), levels(.data$Component)[.data$Component], "")), size = 3,
  #                   position=position_stack(vjust=0.5), colour="black")
  #  }
  #  ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(Identity==.x) %>% select(Cells))))
  #} else {
 #   p <- ggplot(data.plot, aes(.data$Component, .data$Pct))
 #
 #     title <- ggtitle(paste0(component, " cell proportions"))
 #     b <- geom_jitter(aes(colour = .data$Group))
 #     #c <- scale_fill_manual(name = "Identity", values = (myColors[levels(data.plot$Component)]))
 #     #c <- scale_fill_manual(name = "Identity", values = myColors)
 #     #t <- geom_text(aes(label=ifelse((do.label==T & .data$Pct >= (ymax/10)), levels(.data$Component)[.data$Component], "")), size = 3,
 #     #               position=position_stack(vjust=0.5), colour="black")
 #   }
    #ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(Identity==.x) %>% select(Pct))))
 # }



  ggplot(data.plot, aes(Group, Pct)) +
    geom_boxplot(aes(colour = Group), outlier.shape = NA, position = position_dodge(width=1), na.rm = F) +
    geom_jitter(aes(colour = Group), position = position_dodge(width = 1)) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    theme(legend.title=element_text(size=rel(1.1))) +
    #scale_y_continuous(limits = c(0,ymax+10), expand = c(0, 0)) +
    if (column==T){
      facet_wrap("Component", ncol = 1, scales = scales)
    } else {
      facet_wrap("Component", nrow = 1, scales = scales)
    }

}


#' myColours
#'
#' @return returns a named vector of colours
#' @export
#'
#' @examples
myColours <- function(){

  myColors <- gplots::col2hex(c("grey", "grey",  "royalblue1", "brown", "darkgreen", "green", "skyblue","palevioletred4",
                                "peachpuff2", "goldenrod", "tan2", "wheat3",
                                "lightgreen", "palegreen4", "forestgreen", "goldenrod", "tan3", "lightskyblue3", "cyan", "royalblue3", "grey20",
                                "orchid4", "orchid1", "maroon2", "magenta", "mediumpurple2",
                                "orangered1", "wheat3", "goldenrod4"
  ))
  #Create a custom color scale

  names(myColors) <-  c("OffTarget", "unassigned", "Endo", "Stroma", "NPC-like", "NPC", "Nephron", "UrEp",
                        "EN", "DN", "PN", "RC",
                        "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC",
                        "SPC", "CS", "MS", "MesS", "Stroma_NC",
                        "UTip", "UOS", "UIS"
  )
  return(myColors)
}


GetSampleIDs <- function(sample = NULL,
                         reference = NULL,
                         type = NULL,
                         base_protocol = NULL,
                         modification = NULL,
                         age = NULL,
                         union = FALSE){

  if (!exists(x = "sample_table")){
    load(file = here::here("data/sample_table.rda"))
  }

  filter.list <- list()

  if (!is.null(sample)){
    filter.list[["f.sample"]] <- sample_table[sample_table$sample %in% sample,]
  }

  if (!is.null(reference)){
    filter.list[["f.reference"]] <- sample_table[sample_table$reference %in% reference,]
  }

  if (!is.null(type)){
    filter.list[["f.type"]] <- sample_table[sample_table$type %in% type,]
  }

  if (!is.null(base_protocol)){
    filter.list[["f.base_protocol"]] <- sample_table[sample_table$base_protocol %in% base_protocol,]
  }

  if (!is.null(modification)){
    filter.list[["f.modification"]] <- sample_table[sample_table$modification %in% modification,]
  }

  if (!is.null(age)){
    filter.list[["f.age"]] <- sample_table[sample_table$age %in% age,]
  }

  if (union){
    sample.ids <- Reduce(union, map(filter.list, ~.x$sample))
  } else {
    sample.ids <- Reduce(intersect, map(filter.list, ~.x$sample))
  }

  return(sample.ids)
}


