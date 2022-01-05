#' ComponentPlot
#'
#' @param data seurat object
#' @param identity column name in metadata for sample identity
#' @param component column name in metadata for annotation/classification identity
#' @param feature gene for expression visualisation, only shown if show.gene.exp is TRUE
#' @param show.unassigned include unassigned cells in visualisation, TRUE by default
#' @param show.pct change plot to show percentage of cells instead of total cell number
#' @param show.gene.exp colour by expression of "feature" in each annotation category
#' @param do.label add label to plot, TRUE by default
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' ComponentPlot(organoid, show.pct = TRUE)
ComponentPlot <- function(data, identity = "orig.ident", component = "DKCC", feature = "MALAT1",
                          show.unassigned = TRUE, show.pct = FALSE, show.gene.exp = FALSE, do.label = T){
  data.plot <- GeneSummary(data, identity = identity, split.by = component, features = feature)

  myColors <- myColours()

  data.plot <- data.plot %>% filter(.data$Cells>0)

  data.plot$Component <- factor(data.plot$Component, levels = names(myColors)[names(myColors) %in% unique(as.character(data.plot$Component))])
  if (!is.null(levels(data.plot$Identity))){
    data.plot$Identity <- factor(data.plot$Identity, levels = levels(data@meta.data[, identity]))
  }


  data.plot <- data.plot %>% arrange(desc(.data$Component))
  if (show.unassigned == FALSE) {
    data.plot <- data.plot %>% filter(.data$Component != "unassigned")
  }

  myColors <- myColors[levels(data.plot$Component)]

  fct.order <- levels(data@meta.data[[identity]])
  data.plot <- data.plot %>% mutate(Identity = forcats::fct_relevel(.data$Identity, fct.order))


  if (show.pct == FALSE){
    p <- ggplot(data.plot, aes(.data$Identity, .data$Cells))

    if (show.gene.exp == T) {
      title <- ggtitle(paste0(feature, " expression across ", component))
      b <- geom_bar(aes(fill = .data$avg.exp.log), stat = "Identity", colour = "black", width = 0.99)
      max <- max(data.plot %>% filter(.data$features.plot == feature) %>% select(.data$avg.exp.log))


      c <- scale_fill_gradient2(low = "lightgrey", mid = "red", high = "green", midpoint = max/2, na.value = "black")
      t <- geom_text(aes(label=ifelse((do.label==T & .data$avg.exp.log >= (max/3)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    } else {
      title <- ggtitle(paste0(component, " cell numbers"))
      b <- geom_bar(aes(fill = .data$Component), stat = "Identity", colour = "black", width = 0.99)

      c <- scale_fill_manual(name = "Identity", values = myColors)
      t <- geom_text(aes(label=ifelse((do.label==T & .data$Cells >= (ymax/10)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    }
    ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(.data$Identity==.x) %>% select(.data$Cells))))
  } else {
    p <- ggplot(data.plot, aes(.data$Identity, .data$Pct))
    if (show.gene.exp == T) {
      title <- ggtitle(paste0(feature, " expression across ", component))
      b <- geom_bar(aes(fill = .data$avg.exp.log), stat = "Identity", colour = "black", width = 0.99)
      max <- max(data.plot %>% filter(.data$features.plot == feature) %>% select(.data$avg.exp.log))


      c <- scale_fill_gradient2(low = "lightgrey", mid = "red", high = "green", midpoint = max/2, na.value = "black")
      t <- geom_text(aes(label=ifelse((do.label==T & .data$avg.exp.log >= (max/3)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    } else {
      title <- ggtitle(paste0(component, " cell proportions"))
      b <- geom_bar(aes(fill = .data$Component), stat = "Identity", colour = "black", width = 0.99)

      c <- scale_fill_manual(name = "Identity", values = myColors)
      t <- geom_text(aes(label=ifelse((do.label==T & .data$Pct >= (ymax/10)), levels(.data$Component)[.data$Component], "")), size = 3,
                     position=position_stack(vjust=0.5), colour="black")
    }
    ymax <- max(map_dbl(unique(data.plot$Identity), ~sum(data.plot %>% filter(.data$Identity==.x) %>% select(Pct))))
  }

  p +
    b +
    t +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1)) +
    theme(legend.title=element_text(size=rel(1.1))) +
    scale_y_continuous(limits = c(0,ymax+1), expand = c(0, 0)) +
    c +
    title
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
#' @param scaling choose "log", "scale" or "raw" to change how the expression values are displayed
#' @param classification choose "DKCC" or "LineageID"
#' @param col.min minimum expression value to plot
#' @param col.max maximum expression value to plot
#' @param size how to scale the size of the dot, "pct.exp" for the percent of cells expressing the gene, "Pct" is percent of cells in identity
#' @param compare.to.organoids if function should plot organoid samples from the database
#' @param filter.samples if compare.to.organoids=T, can provide desired sample names here using the output from "GetSampleIDs".
#' @param columns how many columns to plot. Default is to plot as square as possible.
#'
#' @return ggplot object
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
    comp.data <- comp.data %>% filter(.data$Identity %in% filter.samples)
  }

  samples <- comp.data$Identity %>% unique()

  #features <- factor(features, levels = features)
  if (classification == "DKCC"){
    identity <- factor(c("Nephron_NC", "Stroma_NC", "NPC-like",
                         unique(comp.data$Component)), levels = c("NPC", "NPC-like", "EN", "EDT", "DT", "LOH",
                                                           "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC",
                                                           "UTip", "UOS", "UIS", "Stroma_NC",
                                                           "SPC", "CS", "MS", "MesS", "Endo", "unassigned"))
  } else if (classification == "LineageID"){
    identity <- factor(c(unique(comp.data$Component), "unassigned", "NPC-like"), levels = c("Endo", "Stroma", "NPC", "NPC-like", "Nephron", "UrEp", "unassigned"))
  } else {
    identity <- factor(c(unique(comp.data$Component), "unassigned"), levels = c("EN", "DN", "PN", "RC", "unassigned"))
  }

  if(idents == "all"){
    identity <- identity
  } else if (idents == "nephron"){
    identity <- identity[identity %in% c("NPC", "NPC-like", "EN", "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod")]
  } else if (idents == "ureteric"){
    identity <- identity[identity %in% c("UTip", "UIS", "UOS")]
  } else if (idents == "others"){
    identity <- identity[identity %in% c("SPC", "CS", "MS", "MesS", "Endo", "unassigned")]
  } else if (idents %in% identity){
    identity <- identity[identity %in% idents]
  } else {
    warning("This call needs to specify identities eg. all, nephron, ureteric, others, or specific identities")
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
    filter(.data$Component %in% identity, !is.na(.data$features.plot)) %>%
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

    theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1),
          panel.grid.major = element_line(colour = "gray"),
          #panel.background = element_rect(colour = "lightgray")
    ) +
    theme(legend.title=element_text(size=rel(0.5))) +

    labs(x = "Features", y = "Identity") +
    facet_wrap(facets = facet, ncol = ifelse(is.null(columns), round(sqrt(length(unique(data.plot[, facet])))), columns)) +
    geom_vline(xintercept = c(seq(0.5, (length(features)+0.5), 1)), col="black")
  return(plot)
}




#' IdentMeans
#'
#' @param data seurat object with DevKidCC annotation
#' @param group identity to generate mean scores from - x axis
#' @param identity individual samples
#' @param component class of annotation of interest
#'
#' @return
#' @export
#'
#' @examples
IdentMeans <- function(data, group = "protocol", identity = "orig.ident", component = "DKCC"){

  myColors <- myColours()
  #data.plot <- data.plot %>% filter(Cells>0)

  plot.data <- data@meta.data
  plot.data <- plot.data %>% filter(!is.na(component), !is.na(identity), !is.na(group))

  plot.components <- unique(plot.data[, component])
  plot.idents <- unique(plot.data[, identity])
  plot.groups <- unique(plot.data[, group])


  abc <- plot.data %>% transmute(identity = plot.data[,identity],
                                 component = factor(plot.data[, component], levels = unique(plot.data[, component])),
                                 group = plot.data[, group])

  total.n <- abc %>%
    group_by(identity) %>%
    summarise(total.cells = n())


  identity.n <- map_df(plot.groups, ~as.data.frame(table((abc %>% filter(group == .x))$identity, (abc %>% filter(group == .x))$component)) %>%
                         `colnames<-`(c("identity", "component", "ident.cells")))

  identity.group <- as.data.frame(table(plot.data[,identity], plot.data[,group])) %>% filter(Freq>0) %>% transmute(identity = Var1, group = Var2)

  pct.n <- left_join(total.n, identity.n, by = "identity") %>% mutate(percentage = (ident.cells/total.cells*100) %>% round(digits = 2)) %>%
    left_join(identity.group, by = "identity") %>% filter(!is.na(percentage))

  summary.n <- map_df(plot.groups, ~pct.n %>% group_by(group, component) %>% filter(group == .x) %>%
                        summarise(n = n(),
                                  mean = mean(percentage),
                                  sd = sd(percentage)) %>% mutate(sem = sd/sqrt(n)))


  myColors <- myColors[names(myColors) %in% unique(summary.n$component)]

  fct.order <- names(myColors)
  plot.components <- fct.order
  summary.n <- summary.n %>% mutate(component = forcats::fct_relevel(component, fct.order))
  pct.n <- pct.n %>% mutate(component = forcats::fct_relevel(component, fct.order))


  #t.test((pct.n %>% filter(component=="Stroma", group == "Freedman"))$percentage, (pct.n %>% filter(component=="Stroma", group == "Morizane"))$percentage)

  map(plot.components, ~summary.n %>% filter(component==.x) %>% ggplot() +
        geom_point(aes(group, mean, colour = group, group = group), size = 2, position = position_dodge(width = 0.5)) +
        geom_errorbar(aes(x=group, ymin=mean-sem, ymax=mean+sem, colour=group), width=0.5, alpha=1, size=.75, position = position_dodge(width = 0.5)) +
        geom_point(inherit.aes = F, data = (pct.n %>% filter(component==.x)),
                   aes(group, percentage, colour = group, group = group), size = 1, alpha = 0.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1)) +
        theme(legend.title=element_text(size=rel(1.1)),
              axis.title = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle(label = .x) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,NA))) %>%
    patchwork::wrap_plots(nrow = 1, guides = "collect")

}


#' myColours
#'
#' @return returns a named vector of colours
#' @export
#'
#' @examples
#' myColours()
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
  levels(myColors) <-  c("OffTarget", "unassigned", "Endo", "Stroma", "NPC-like", "NPC", "Nephron", "UrEp",
                        "EN", "DN", "PN", "RC",
                        "EDT", "DT", "LOH", "EPT", "PT", "PEC", "EPod", "Pod", "Nephron_NC",
                        "SPC", "CS", "MS", "MesS", "Stroma_NC",
                        "UTip", "UOS", "UIS"
  )
  return(myColors)
}


#' GetSampleIDs
#'
#' @param sample Not required
#' @param reference Original publication reference eg. Howden_2019
#' @param type can select either "Organoid" or "iUB" if required
#' @param base_protocol Select base protocol/s of interest, including "Takasato", "Morizane", "Freedman", "Kumar", "Low", "Howden-Wilson", "Mae"
#' @param modification TRUE/FALSE for if that base protocol was modified in some way
#' @param age age range of organoids, string of integers eg. 1:18
#' @param union if multiple fields are entered, use union=TRUE to select all samples or union=FALSE to select only the intersect i.e. samples selected from all fields
#'
#' @return string of sample IDs that can be input into DotPlotCompare.
#' @export
#'
#' @examples
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
  sample.ids <- c(as.character(sample.ids), "Reference")

  return(sample.ids)
}


