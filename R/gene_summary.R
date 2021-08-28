#' GeneSummary
#'
#' @param data seurat object
#' @param identity column name in metadata describing required sample identity
#' @param split.by column name in metadat describing required annotation/classification identity
#' @param cells cells to use, defaults to all
#' @param do.norm option to run normalisation on data if not already done, set TRUE if required
#' @param features features to use, defaults to all
#' @param group.by secondary grouping identity, default is NULL
#'
#' @export
#'
#' @examples
#' output <- GeneSummary(organoid, features = c("GAPDH", "SIX2", "JAG1", "EPCAM", "NPHS1"))
GeneSummary <-function(data,
                       identity = "orig.ident",
                       split.by = "DKCC",
                       #idents = NULL,
                       features = rownames(data),
                       cells = colnames(data),
                       do.norm = FALSE,
                       group.by = NULL
){
  # gene proportion information
  Idents(data) <- identity
  df <- as.data.frame(table(data@meta.data[[paste0(identity)]],
                            data@meta.data[[paste0(split.by)]]))
  id <- as.data.frame(table(data[[paste0(identity)]]))
  colnames(id) <- c("Identity", "CellTotal")
  colnames(df) <- c("Identity", "Component", "Cells")
  df <- left_join(df, id, "Identity")
  df$Pct <- round((df$Cells / df$CellTotal * 100), digits = 2)

  # gene expression information

  #cells <- unlist(x = CellsByIdentities(object = data, idents = NULL))
  if (do.norm==T){data <- NormalizeData(data)}
  data.features <- FetchData(object = data, vars = features,
                             cells = cells, slot = "data")

  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = data)[cells, drop = TRUE]
  }  else {
    data[[group.by, drop = TRUE]][cells, drop = TRUE]
  }

  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }

  data.features$id <- as.vector(x = data.features$id) #possibly counter productive to previous call but leave it in anyway

  splits <- data[[split.by, drop = TRUE]][cells, drop = TRUE]
  # rename id to be "sample___identity"
  data.features$id <- paste(data.features$id, splits,
                            sep = "___")

  # this function required for next section
  PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
  }

  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
                     threshold = 0)

    cell.pct <- apply(X = data.use, MARGIN = 2, FUN = length)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp, cell.pct = cell.pct))
  })

  names(x = data.plot) <- unique(x = data.features$id)

  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })

  data.plot <- do.call(what = "rbind", args = data.plot)

  if (length(x = unique(x = data.plot$id)) < 3) {
    data.transform <- ""
    warning("Only less than three identities present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE)
  }

  avg.exp.log <- sapply(X = unique(x = data.plot$features.plot),
                        FUN = function(x) {
                          data.use <- data.plot[data.plot$features.plot ==
                                                  x, "avg.exp"]
                          data.use <- log1p(x = data.use)
                          return(data.use)
                        })
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot ==
                                                     x, "avg.exp"]

                             data.use <- scale(x = data.use)

                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  avg.exp.log <- as.vector(x = t(x = avg.exp.log))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$avg.exp.log <- avg.exp.log
  data.plot$features.plot <- factor(x = data.plot$features.plot,
                                    levels = features)
  data.plot$pct.exp <- data.plot$pct.exp * 100

  if (!is.null(split.by)){
    sample.id <- stringr::str_split(data.plot$id, pattern = "___", simplify = T)
    data.plot$id <- sample.id[,1]
    data.plot$sample <- sample.id[,2]
  }

  #if (is.null(levels(data@meta.data[[split.by]]))){
  #  data.plot$sample <- factor(data.plot$sample, levels = as.character(unique(data@meta.data[[split.by]])))
  #} else {
  #  data.plot$sample <- factor(data.plot$sample, levels = levels(data@meta.data[[split.by]]))
  #}

  data.plot <- left_join(df, data.plot, by = c("Component" = "sample", "Identity" = "id")) %>% select(c("Identity", "Component", "Cells", "Pct", "avg.exp", "pct.exp",
                                                                                                        "features.plot", "avg.exp.scaled", "avg.exp.log"))

  data.plot$Component <- data.plot$Component %>% factor(levels = unique(data.plot$Component))

  return(data.plot)
}



