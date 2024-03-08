#' Check If Seurat Object Has Been Normalized
#'
#' This function tests whether a given Seurat object has undergone normalization.
#' Normalization is assumed to have occurred if the `data` slot of the Seurat
#' object contains values. The function checks for the existence and non-emptiness
#' of the `data` slot as an indicator of normalization.
#'
#' @param seuratObj A Seurat object.
#'
#' @return Logical indicating whether the Seurat object has been normalized (`TRUE`)
#' or not (`FALSE`). 
#'
#' @examples
#' # Assuming `seuratObj` is your Seurat object
#' if (isNormalizedSeurat(seuratObj)) {
#'   print("The Seurat object has been normalized.")
#' } else {
#'   print("The Seurat object has not been normalized.")
#' }
#'
#' @export
isNormalized <- function(seurat_object) {
  activeAssay <- Seurat::DefaultAssay(seurat_object)
  
  !is.null(seurat_object[[activeAssay]]$data) && length(seurat_object[[activeAssay]]$data) > 0
}

#' Load a Seurat object from a file or variable
#'
#' This function loads a Seurat object from a file (.qs or .rds) or a variable. If the input is not a Seurat object, an error message is displayed.
#'
#' @param seurat_object A Seurat object or the path to a file containing a Seurat object (.qs or .rds).
#'
#' @return A Seurat object.
#'
#' @importFrom qs qread
#' @importFrom glue glue
#' @importFrom Seurat Seurat
#'
#' @examples
#' # Load a Seurat object from a .qs file
#' seurat_obj <- load_seurat("path/to/seurat_object.qs")
#'
#' # Load a Seurat object from a variable
#' seurat_obj <- create_seurat_object()
#' seurat_obj <- load_seurat(seurat_obj)
#'
#' @export
load_seurat <- function(seurat_object) {
  input_name <- deparse(substitute(seurat_object))
  if (!inherits(seurat_object, "Seurat")) {
    library("Seurat") |> suppressPackageStartupMessages()
    if (grepl("\\.qs$", seurat_object)) {
      seurat_object <- qs::qread(seurat_object)
    } else if (grepl("\\.rds$", seurat_object)) {
      seurat_object <- readRDS(seurat_object)
    } else {
      stop("Invalid Seurat object file format. Please provide the path to a .qs or .rds file.")
    }
    #Test if the object is a Seurat object. If not, show the name of variable forvided to the function and the message "...is not a Seurat object."
    if (!inherits(seurat_object, "Seurat")) {
    stop(glue::glue("{input_name} is not a Seurat object."))
    }   
  }
  return(seurat_object)
}
