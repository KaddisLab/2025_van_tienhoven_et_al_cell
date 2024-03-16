#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param cellbender_h5
#' @param cellranger_run_folder
#' @return
#' @author Denis O'Meally
#' @export
make_seurat_cellbender <- function(cellbender_h5, cellranger_run_folder) {
    
    sample_id <- basename(cellranger_run_folder)

    seurat_object <- scCustomize::Create_CellBender_Merged_Seurat(
        scCustomize::Read_CellBender_h5_Mat(cellbender_h5),
        Seurat::Read10X_h5(glue::glue("{cellranger_run_folder}/outs/filtered_feature_bc_matrix.h5"))
        ) |> suppressWarnings()

    # Set orig.ident
    seurat_object@meta.data$orig.ident <- sample_id

    # Add cellbender cols    
    seurat_object <- scCustomize::Add_CellBender_Diff(  
        seurat_object = seurat_object,
        raw_assay_name = "RAW",
        cell_bender_assay_name = "RNA") 
    
    # Add QC metrics
    seurat_object <- seurat_add_cell_metrics(seurat_object)

    # Save the Seurat object
    seurat_object_path <- glue::glue("{analysis_cache}/cellbender_seurat_objects/{sample_id}.qs")
    dir.create(dirname(seurat_object_path), recursive = TRUE, showWarnings = FALSE)

    qs::qsave(seurat_object, file = seurat_object_path)

    return(seurat_object_path)
}
