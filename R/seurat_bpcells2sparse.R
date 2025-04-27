#' Convert BPCells Matrices to Sparse Matrices in Seurat Object
#'
#' Converts any BPCells matrices in a Seurat object to dgCMatrix format, which is
#' more compatible with standard Seurat operations. The conversion is performed
#' in-place for specified assay layers.
#'
#' @param so A Seurat object or path to a Seurat object file
#' @param assay Name of the assay to process (default: "RNA")
#' @param layers Character vector of layer names to convert (default: c("counts", "data", "scale.data"))
#'
#' @return Path to saved Seurat object (.qs format) with converted matrices
#'
#' @details
#' The function:
#' 1. Loads the Seurat object if path provided
#' 2. Processes each specified layer in the assay
#' 3. Converts BPCells matrices to dgCMatrix format
#' 4. Saves the modified object to "{analysis_cache}/bpcells2sparse_out/{project_name}.qs"
#'
#' Memory usage is monitored and reported during conversion.
#'
#' @import Seurat SeuratObject
#' @importFrom methods as
#' @importFrom pryr mem_used
#' @importFrom glue glue
#' @export
seurat_bpcells2sparse <- function(
    so,
    assay = "RNA",
    layers = c("counts", "data", "scale.data")
) {
    message("This is seurat_bpcells2sparse()...\n")
    message("Initial memory: ", pryr::mem_used() / 1e9, " GB")
    so <- load_seurat(so)
    project_name <- Seurat::Project(so)
    message("Project name: ", project_name)

    # Get the assay object directly
    assay_obj <- so[[assay]]
    message(
        "Available layers: ",
        paste(SeuratObject::Layers(assay_obj), collapse = ", ")
    )

    # Helper function to convert if needed and update the Seurat object
    convert_if_bpcells <- function(seurat_obj, assay_name, layer_name) {
        layer_data <- Seurat::GetAssayData(
            seurat_obj[[assay_name]],
            layer = layer_name
        )

        if (is.null(layer_data)) {
            return(seurat_obj)
        }

        if (
            inherits(layer_data, "BPCells") ||
                inherits(layer_data, "RenameDims")
        ) {
            message(sprintf("Converting %s matrix to dgCMatrix...", layer_name))
            converted_data <- methods::as(layer_data, "dgCMatrix")
            seurat_obj[[assay_name]] <- Seurat::SetAssayData(
                seurat_obj[[assay_name]],
                layer = layer_name,
                new.data = converted_data
            )
            gc()
            rm(layer_data, converted_data)
        } else {
            message(sprintf(
                "Layer %s is already in correct format, skipping...",
                layer_name
            ))
        }

        return(seurat_obj)
    }

    # Process specified layers
    for (layer_name in layers) {
        if (layer_name %in% SeuratObject::Layers(assay_obj)) {
            message(sprintf("\nChecking %s...", layer_name))
            so <- convert_if_bpcells(so, assay, layer_name)
        } else {
            message(sprintf("\nLayer '%s' not found, skipping...", layer_name))
        }
    }

    message("\nFinal memory: ", pryr::mem_used() / 1e9, " GB")
    #--------------------------------------------------------------------------------
    # Save results
    object_path <- (glue::glue("{analysis_cache}/bpcells2sparse_out/{project_name}.qs"))
    message("Writing to disk at: ", object_path)
    save_results(so, object_path)
    return(object_path)
}
