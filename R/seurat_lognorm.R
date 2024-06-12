#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param seurat_objects_merged
#' @return
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

