#' Perform Log-Normalization and Scaling on a Seurat Object
#'
#' Applies standard log-normalization, identifies variable features, and performs
#' scaling with optional regression of specified variables. The process is performed
#' separately for each group if split_by is specified.
#'
#' @param seurat_object A Seurat object or path to a Seurat object file
#' @param assay Name of assay to process (default: "RNA")
#' @param vars_to_regress Variables to regress out during scaling (default: "batch")
#' @param split_by Column name to split data by during scaling (default: "orig.ident")
#'
#' @return A Seurat object with normalized and scaled data. Project name is updated
#'         to include "_lognorm" suffix.
#'
#' @details
#' The function performs three main steps:
#' 1. Log-normalization using NormalizeData()
#' 2. Variable feature selection using FindVariableFeatures()
#' 3. Data scaling with regression using ScaleData()
#'
#' @import Seurat
#' @importFrom future plan
#' @importFrom glue glue
#' @author Denis O'Meally
#' @export
seurat_lognorm <- function(seurat_object, assay = "RNA", vars_to_regress = "batch", split_by = "orig.ident") {

    future::plan("multisession")
    options(future.globals.maxSize = hprcc::slurm_allocation()$Memory_GB * 1024^3)

    message("This is seurat_lognorm()...\n")
    seurat_object <- load_seurat(seurat_object)

    project_name <- Seurat::Project(seurat_object)

    new_project_name <- glue::glue("{project_name}_lognorm")
   
    Seurat::DefaultAssay(seurat_object) <- assay

    seurat_object <- seurat_object |>
        Seurat::NormalizeData(normalization.method = "LogNormalize", verbose = TRUE) |>
        Seurat::FindVariableFeatures() |>
        Seurat::ScaleData(vars.to.regress = vars_to_regress, split.by = split_by, verbose = TRUE)

    Seurat::Project(seurat_object) <- new_project_name
    #--------------------------------------------------------------------------------
    # Save results
    # output_path <- (glue::glue("{analysis_cache}/lognorm_out/{new_project_name}.qs"))
    # save_results(seurat_object, output_path)
    return(seurat_object)
}
