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

