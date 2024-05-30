#' Create a Seurat object from CellBender and CellRanger data
#'
#' This function reads HDF5 formatted output from CellBender and a CellRanger run folder
#' to create a merged Seurat object. The function also adds QC metrics and renames cells
#' using sample identifiers, and filters out cells with fewer than 200 features.
#'
#' @param cellbender_h5 The file path to the CellBender output in HDF5 format.
#' @param cellranger_run_folder The directory of the CellRanger run.
#' @param sample_metadata A data frame containing sample metadata. Must include the columns `sample_id` (library_id) and `sample_name` (orig.ident).
#' @param add_cols Additional columns to add to the Seurat object metadata from the sample metadata, passed to select(). Default is `NULL`.
#'
#' @return The file path of the saved Seurat object in .qs format.
#'
#' @examples
#' make_seurat_cellbender("path/to/cellbender_output.h5", "path/to/cellranger_folder")
#'
#' @importFrom Seurat Read10X_h5 RenameCells
#' @importFrom glue glue
#' @importFrom qs qsave
#' @import scCustomize Create_CellBender_Merged_Seurat Read_CellBender_h5_Mat Add_CellBender_Diff
#' @export
make_seurat_cellbender <- function(cellbender_h5, cellranger_run_folder, sample_metadata, add_cols = NULL) {
    sample_id <- basename(cellranger_run_folder)
    sample_name <- sample_metadata$sample_name[which(sample_metadata$sample_id == sample_id)]

    seurat_object <- scCustomize::Create_CellBender_Merged_Seurat(
        min_features = 200,
        scCustomize::Read_CellBender_h5_Mat(cellbender_h5),
        Seurat::Read10X_h5(glue::glue("{cellranger_run_folder}/outs/filtered_feature_bc_matrix.h5"))
    ) |> suppressWarnings()

    # Set orig.ident
    seurat_object@meta.data$orig.ident <- sample_name
    Seurat::Project(seurat_object) <- sample_id

    # Add cellbender cols
    seurat_object <- scCustomize::Add_CellBender_Diff(
        seurat_object = seurat_object,
        raw_assay_name = "RAW",
        cell_bender_assay_name = "RNA"
    )

    # Add QC metrics
    seurat_object <- seurat_add_cell_metrics(seurat_object)

    # Rename cells
    seurat_object <- Seurat::RenameCells(seurat_object, add.cell.id = sample_id)

    meta_data <- dplyr::select(sample_metadata, c(sample_id, sample_name, !!!add_cols))

    # Add sample metadata
    seurat_object <- scCustomize::Add_Sample_Meta(seurat_object = seurat_object, meta_data = meta_data, "orig.ident", "sample_name")


    # Save the Seurat object
    seurat_object_path <- glue::glue("{analysis_cache}/cellbender_seurat_objects/{sample_id}.qs")
    dir.create(dirname(seurat_object_path), recursive = TRUE, showWarnings = FALSE)

    qs::qsave(seurat_object, file = seurat_object_path)

    return(seurat_object_path)
}
