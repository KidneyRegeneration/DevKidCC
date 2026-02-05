# Seurat v5 Compatibility Layer
# Auto-generated compatibility functions for DevKidCC

#' Check if assay is Seurat v5 (Assay5 class)
#' @param object Seurat object
#' @param assay Assay name (default: NULL = DefaultAssay)
#' @return Logical
#' @keywords internal
is_v5_assay <- function(object, assay = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  
  if (!assay %in% names(object@assays)) {
    stop(paste("Assay", assay, "not found"))
  }
  
  inherits(object[[assay]], "Assay5")
}

#' Get expression data compatible with v4 and v5
#' @param object Seurat object
#' @param layer Layer/slot name ("data", "counts", or "scale.data")
#' @param assay Assay name (default: NULL = DefaultAssay)
#' @return Expression matrix
#' @keywords internal
#' @export
get_layer_data <- function(object, layer = "data", assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  
  # Try new method first (v5)
  if (is_v5_assay(object, assay)) {
    data <- SeuratObject::LayerData(object, layer = layer, assay = assay)
  } else {
    # Fallback to compatible method
    data <- Seurat::GetAssayData(object, layer = layer, assay = assay)
  }
  
  return(data)
}

#' Set expression data compatible with v4 and v5
#' @param object Seurat object
#' @param layer Layer/slot name
#' @param value New data matrix
#' @param assay Assay name
#' @return Updated Seurat object
#' @keywords internal
#' @export
set_layer_data <- function(object, layer = "data", value, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  
  if (is_v5_assay(object, assay)) {
    # v5 method
    SeuratObject::LayerData(object, layer = layer, assay = assay) <- value
  } else {
    # v4 method
    object <- Seurat::SetAssayData(
      object,
      layer = layer,
      new.data = value,
      assay = assay
    )
  }
  
  return(object)
}

#' Get variable features
#' @keywords internal
#' @export
get_var_features <- function(object, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(object)
  SeuratObject::VariableFeatures(object, assay = assay)
}

#' Set variable features
#' @keywords internal
#' @export
set_var_features <- function(object, features, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(object)
  SeuratObject::VariableFeatures(object, assay = assay) <- features
  return(object)
}

#' Get reduction embeddings
#' @keywords internal
#' @export
get_embeddings <- function(object, reduction = "umap") {
  if (!reduction %in% names(object@reductions)) {
    stop(paste("Reduction", reduction, "not found"))
  }
  SeuratObject::Embeddings(object, reduction = reduction)
}

#' Check Seurat version and warn if issues
#' @keywords internal
check_seurat_version <- function() {
  version <- packageVersion("Seurat")
  
  if (version < "5.0.0") {
    warning(
      "DevKidCC is optimized for Seurat v5.x. ",
      "You are using v", version, ". ",
      "Some features may not work correctly."
    )
  }
  
  invisible(version)
}

#' Update old Seurat object to current version
#' @param object Seurat object
#' @return Updated object
#' @export
update_seurat_object <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  old_version <- if (!is.null(object@version)) {
    as.character(object@version)
  } else {
    "unknown"
  }

  message("Original object version: ", old_version)

  # Update using Seurat function
  if (exists("UpdateSeuratObject", where = asNamespace("Seurat"))) {
    object <- Seurat::UpdateSeuratObject(object)
    message("Updated to Seurat v", packageVersion("Seurat"))
  } else {
    message("UpdateSeuratObject not available")
  }

  return(object)
}

#' Join layers in Seurat v5 Assay5 objects if needed
#' @param object Seurat object
#' @param assay Assay name (default: NULL = DefaultAssay)
#' @param layers Type of layers to join: "data", "counts", or "scale.data"
#' @return Updated Seurat object with joined layers
#' @export
join_layers_if_needed <- function(object, assay = NULL, layers = "data") {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }

  # Only process if it's an Assay5 object
  if (!inherits(object[[assay]], "Assay5")) {
    return(object)
  }

  # Check for multiple layers
  layer_names <- SeuratObject::Layers(object, search = layers)

  if (length(layer_names) > 1) {
    message("Joining ", length(layer_names), " '", layers, "' layers in assay '", assay, "'")
    object[[assay]] <- SeuratObject::JoinLayers(object[[assay]])
  }

  return(object)
}

