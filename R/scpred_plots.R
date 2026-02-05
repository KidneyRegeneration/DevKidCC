# scPred Visualization Functions
# R versions of Python plotting functions from DevKidCC-python

#' Plot scPred Scores Distribution
#'
#' Generates a boxplot/violin plot combo for prediction scores, optionally split by a metadata variable
#'
#' @param seurat Seurat object with scPred results
#' @param score_cols Character vector of scPred score column names (e.g., c("scpred_NPC", "scpred_EN"))
#' @param split_by Optional metadata column name to split plots by (e.g., "sample", "orig.ident")
#' @param plot_type Type of plot: "box" (boxplot), "violin" (violin plot), or "both" (default)
#' @param colors Optional named vector of colors for cell types
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' PlotScPredScoresDistribution(
#'   seurat = organoid,
#'   score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN"),
#'   split_by = "sample"
#' )
PlotScPredScoresDistribution <- function(seurat,
                                          score_cols,
                                          split_by = NULL,
                                          plot_type = "both",
                                          colors = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }

  # Extract metadata with scPred scores
  df <- seurat[[]]

  # Check if score columns exist
  missing_cols <- score_cols[!score_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing score columns: ", paste(missing_cols, collapse = ", "))
  }

  # Extract scores and reshape to long format
  score_data <- df[, score_cols, drop = FALSE]

  # Calculate max score and highest type for each cell
  score_data$Max_Value <- apply(score_data[, score_cols], 1, max, na.rm = TRUE)
  score_data$Highest_Type <- score_cols[apply(score_data[, score_cols], 1, which.max)]
  score_data$Highest_Type <- gsub("scpred_", "", score_data$Highest_Type)

  # Add split_by variable if specified
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(df)) {
      stop("Split variable '", split_by, "' not found in metadata")
    }
    score_data[[split_by]] <- df[[split_by]]
  }

  # Order cell types
  ordered_types <- gsub("scpred_", "", score_cols)
  score_data$Highest_Type <- factor(score_data$Highest_Type, levels = ordered_types)

  # Create base plot
  p <- ggplot(score_data, aes(x = Highest_Type, y = Max_Value))

  # Add geoms based on plot_type
  if (plot_type %in% c("both", "box")) {
    if (!is.null(split_by)) {
      p <- p + geom_boxplot(aes(fill = .data[[split_by]]),
                            outlier.shape = NA,
                            position = position_dodge(width = 0.8),
                            alpha = 0.7)
    } else {
      p <- p + geom_boxplot(fill = "white", outlier.shape = NA)
    }
  }

  if (plot_type %in% c("both", "violin")) {
    if (!is.null(split_by)) {
      p <- p + geom_violin(aes(fill = .data[[split_by]]),
                           position = position_dodge(width = 0.8),
                           alpha = 0.5)
    } else {
      p <- p + geom_violin(fill = "lightgrey", alpha = 0.5)
    }
  }

  # Add points
  if (!is.null(split_by)) {
    p <- p + geom_jitter(aes(color = .data[[split_by]]),
                         position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                         size = 0.5, alpha = 0.3)
  } else {
    p <- p + geom_jitter(width = 0.2, size = 0.5, alpha = 0.3, color = "black")
  }

  # Styling
  plot_title <- if (!is.null(split_by)) {
    paste0("Max Score Distribution (Split by ", split_by, ")")
  } else {
    "Max Score Distribution"
  }

  p <- p +
    labs(title = plot_title,
         x = "Predicted Identity",
         y = "Max scPred Score") +
    ylim(-0.05, 1.05) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")

  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors) +
      scale_color_manual(values = colors)
  }

  return(p)
}


#' Plot scPred Scores on UMAP
#'
#' Creates UMAP visualization colored by scPred scores with custom color gradients
#'
#' @param seurat Seurat object with UMAP and scPred results
#' @param score_cols Character vector of scPred score column names
#' @param colors Character vector of colors for each cell type (should match length of score_cols)
#' @param reduction Name of dimensionality reduction to use (default: "umap")
#' @param pt_size Point size for plotting (default: 1.5)
#' @param alpha Point transparency (default: 0.8)
#' @param low_color Color for low scores (default: "lightgrey")
#' @param ncol Number of columns for facet layout (default: 3)
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' PlotScPredUMAP(
#'   seurat = organoid,
#'   score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN"),
#'   colors = c("green", "black", "orange")
#' )
PlotScPredUMAP <- function(seurat,
                           score_cols,
                           colors,
                           reduction = "umap",
                           pt_size = 1.5,
                           alpha = 0.8,
                           low_color = "lightgrey",
                           ncol = 3) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }

  # Check if reduction exists
  if (!reduction %in% names(seurat@reductions)) {
    stop("Reduction '", reduction, "' not found in Seurat object")
  }

  # Check colors match score_cols
  if (length(colors) != length(score_cols)) {
    stop("Length of 'colors' must match length of 'score_cols'")
  }

  # Extract UMAP coordinates
  umap_coords <- get_embeddings(seurat, reduction = reduction)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  # Extract metadata with scPred scores
  df <- seurat[[]]

  # Check if score columns exist
  missing_cols <- score_cols[!score_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing score columns: ", paste(missing_cols, collapse = ", "))
  }

  # Combine UMAP and scores
  plot_data <- cbind(umap_coords, df[, score_cols, drop = FALSE])

  # Calculate max score and highest type
  plot_data$Max_Value <- apply(plot_data[, score_cols], 1, max, na.rm = TRUE)
  plot_data$Highest_Type <- score_cols[apply(plot_data[, score_cols], 1, which.max)]
  plot_data$Highest_Type <- gsub("scpred_", "", plot_data$Highest_Type)

  # Reshape to long format for faceting
  plot_data_long <- plot_data %>%
    tidyr::pivot_longer(
      cols = all_of(score_cols),
      names_to = "CellType",
      values_to = "Score"
    ) %>%
    mutate(CellType = gsub("scpred_", "", CellType))

  # Create color mapping
  color_map <- setNames(colors, gsub("scpred_", "", score_cols))

  # Filter to only show cells where this is the max type
  plot_data_long <- plot_data_long %>%
    filter(CellType == Highest_Type)

  # Create faceted plot
  p <- ggplot(plot_data_long, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
    geom_point(size = pt_size, alpha = alpha) +
    facet_wrap(~ CellType, ncol = ncol) +
    labs(title = "UMAP Projection by Max scPred Score",
         x = "UMAP 1",
         y = "UMAP 2",
         color = "scPred\nScore") +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )

  # Apply custom color gradient for each facet
  # Note: ggplot2 doesn't easily support per-facet color scales,
  # so we'll use a single gradient from low_color to the dominant color
  p <- p + scale_color_gradient(low = low_color, high = "darkblue", limits = c(0, 1))

 return(p)
}


#' Plot Individual scPred Score UMAPs
#'
#' Creates separate UMAP plots for each scPred score with custom color gradients
#' This version creates individual plots that can be combined with patchwork
#'
#' @param seurat Seurat object with UMAP and scPred results
#' @param score_cols Character vector of scPred score column names
#' @param colors Character vector of colors for each cell type
#' @param reduction Name of dimensionality reduction to use (default: "umap")
#' @param pt_size Point size for plotting (default: 1.5)
#' @param alpha Point transparency (default: 0.8)
#' @param low_color Color for low scores (default: "lightgrey")
#'
#' @return List of ggplot objects (one per score column)
#' @export
#' @import ggplot2
#'
#' @examples
#' plots <- PlotScPredUMAPIndividual(
#'   seurat = organoid,
#'   score_cols = c("scpred_NPC", "scpred_EN"),
#'   colors = c("green", "black")
#' )
#' # Combine with patchwork
#' # library(patchwork)
#' # wrap_plots(plots, ncol = 2)
PlotScPredUMAPIndividual <- function(seurat,
                                      score_cols,
                                      colors,
                                      reduction = "umap",
                                      pt_size = 1.5,
                                      alpha = 0.8,
                                      low_color = "lightgrey") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }

  # Check if reduction exists
  if (!reduction %in% names(seurat@reductions)) {
    stop("Reduction '", reduction, "' not found in Seurat object")
  }

  # Extract UMAP coordinates
  umap_coords <- get_embeddings(seurat, reduction = reduction)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  # Extract metadata with scPred scores
  df <- seurat[[]]

  # Combine UMAP and scores
  plot_data <- cbind(umap_coords, df[, score_cols, drop = FALSE])

  # Create a plot for each score column
  plot_list <- list()

  for (i in seq_along(score_cols)) {
    score_col <- score_cols[i]
    cell_type <- gsub("scpred_", "", score_col)
    high_color <- colors[i]

    p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[score_col]])) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_gradient(low = low_color, high = high_color, limits = c(0, 1)) +
      labs(title = cell_type,
           x = "UMAP 1",
           y = "UMAP 2",
           color = "Score") +
      theme_classic() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "right"
      )

    plot_list[[cell_type]] <- p
  }

  return(plot_list)
}


#' Plot scPred Classification on UMAP (Overlay)
#'
#' Creates a single UMAP plot where cells are colored by their highest scPred prediction.
#' Each cell type is assigned a color based on which score column has the maximum value.
#'
#' @param seurat Seurat object with UMAP and scPred results
#' @param score_cols Character vector of scPred score column names
#' @param colors Character vector of colors for each cell type (should match length of score_cols)
#' @param reduction Name of dimensionality reduction to use (default: "umap")
#' @param pt_size Point size for plotting (default: 1.5)
#' @param alpha Point transparency (default: 0.8)
#' @param label Add cell type labels to plot (default: TRUE)
#' @param label_size Size of labels (default: 5)
#' @param legend_position Position of legend (default: "right")
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' PlotScPredUMAPOverlay(
#'   seurat = organoid,
#'   score_cols = c("scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN"),
#'   colors = c("green", "black", "orange", "red")
#' )
PlotScPredUMAPOverlay <- function(seurat,
                                   score_cols,
                                   colors,
                                   reduction = "umap",
                                   pt_size = 1.5,
                                   alpha = 0.8,
                                   label = TRUE,
                                   label_size = 5,
                                   legend_position = "right") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }

  # Check if reduction exists
  if (!reduction %in% names(seurat@reductions)) {
    stop("Reduction '", reduction, "' not found in Seurat object")
  }

  # Check colors match score_cols
  if (length(colors) != length(score_cols)) {
    stop("Length of 'colors' must match length of 'score_cols'")
  }

  # Extract UMAP coordinates
  umap_coords <- get_embeddings(seurat, reduction = reduction)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  # Extract metadata with scPred scores
  df <- seurat[[]]

  # Check if score columns exist
  missing_cols <- score_cols[!score_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing score columns: ", paste(missing_cols, collapse = ", "))
  }

  # Combine UMAP and scores
  plot_data <- cbind(umap_coords, df[, score_cols, drop = FALSE])

  # Calculate max score and highest type for each cell
  plot_data$Max_Score <- apply(plot_data[, score_cols], 1, max, na.rm = TRUE)
  plot_data$Predicted_Type <- score_cols[apply(plot_data[, score_cols], 1, which.max)]
  plot_data$Predicted_Type <- gsub("scpred_", "", plot_data$Predicted_Type)

  # Create color mapping
  cell_types <- gsub("scpred_", "", score_cols)
  color_map <- setNames(colors, cell_types)

  # Order cell types for consistent legend
  plot_data$Predicted_Type <- factor(plot_data$Predicted_Type, levels = cell_types)

  # Calculate centroids for labels
  if (label) {
    centroids <- plot_data %>%
      group_by(Predicted_Type) %>%
      summarise(
        UMAP_1 = median(UMAP_1, na.rm = TRUE),
        UMAP_2 = median(UMAP_2, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Predicted_Type)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = color_map, name = "Predicted\nCell Type") +
    labs(title = "UMAP by Highest scPred Prediction",
         x = "UMAP 1",
         y = "UMAP 2") +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = legend_position
    )

  # Add labels if requested
  if (label) {
    p <- p + geom_text(
      data = centroids,
      aes(label = Predicted_Type),
      color = "black",
      size = label_size,
      fontface = "bold",
      show.legend = FALSE
    )
  }

  return(p)
}


PlotScPredUMAP2 <- function(seurat,
                           score_cols,
                           colors,
                           reduction = "umap",
                           pt_size = 1.5,
                           alpha = 0.8,
                           low_color = "lightgrey") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required")
  if (!requireNamespace("ggnewscale", quietly = TRUE)) stop("Package 'ggnewscale' required")

  # 1. Coordinate Extraction
  umap_coords <- as.data.frame(Seurat::Embeddings(seurat, reduction = reduction))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  # 2. Metadata Extraction
  df <- seurat[[]]
  plot_data <- cbind(umap_coords, df[, score_cols, drop = FALSE])
  
  # Determine which cell type has the max score for each cell
  # We use the exact column names from score_cols to avoid matching errors
  plot_data$Highest_Type <- score_cols[max.col(plot_data[, score_cols], ties.method = "first")]

  # Initialize plot with a background layer of all cells in light grey
  # This ensures the plot isn't "empty" and provides context
  p <- ggplot() +
    geom_point(data = plot_data, aes(x = UMAP_1, y = UMAP_2), 
               color = "snow2", size = pt_size, alpha = 0.2) +
    theme_classic() +
    labs(title = "UMAP Projection by Max scPred Score",
         x = "UMAP 1", y = "UMAP 2") +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  # 3. Layer Iteration
  for (i in seq_along(score_cols)) {
    col_name <- score_cols[i]
    display_name <- gsub("scpred_", "", col_name)
    target_color <- colors[i]
    
    # Filter: Only cells where this specific col_name is the highest
    layer_data <- plot_data[plot_data$Highest_Type == col_name, ]
    
    # Debug check: print if a layer is empty
    if(nrow(layer_data) == 0) next 

    p <- p +
      ggnewscale::new_scale_color() + 
      geom_point(data = layer_data, 
                 aes(x = UMAP_1, y = UMAP_2, color = .data[[col_name]]),
                 size = pt_size, 
                 alpha = alpha) +
      scale_color_gradient(low = low_color, 
                           high = target_color, 
                           name = display_name,
                           limits = c(0, 1))
  }

  return(p)
}
