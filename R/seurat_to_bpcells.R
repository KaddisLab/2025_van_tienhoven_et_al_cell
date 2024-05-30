#' Convert Seurat Object to BPCells Format
#'
#' This function takes a Seurat object (or its path) and converts the RNA counts to BPCells format. If the BPCells
#' directory does not exist for the given sample, it creates one and writes the matrix. It then
#' saves the Seurat object in a quicksave (qs) file format and returns the path to this file.
#'
#' @param seurat_object A Seurat object or a path to a Seurat object file that should be converted
#' to the BPCells format. If a path is provided, the function attempts to load the Seurat object
#' from this path.
#'
#' @return A string containing the path to the saved qs file of the Seurat object after conversion.
#'
#' @examples
#' if_interactive({
#'     seurat_obj <- Seurat::CreateSeuratObject(counts = your_counts_data)
#'     # Assume 'analysis_cache' is a predefined variable or replace with actual path
#'     seurat_obj_path <- seurat_to_bpcells(seurat_obj)
#'     print(seurat_obj_path)
#' })
#'
#' @importFrom glue glue
#' @importFrom qs qsave
#' @importFrom BPCells write_matrix_dir open_matrix_dir
#' @export
#'
#' @details
#' The function first checks if the Seurat object is already loaded or needs to be loaded from a file.
#' It then identifies the sample ID from the Seurat object's metadata. It uses this ID to construct
#' a path where the BPCells formatted data will be stored. If necessary, it converts and saves the
#' RNA counts data to BPCells format in a directory specific to the sample ID. Finally, it saves
#' the Seurat object, now containing a reference to the BPCells data, to a `.qs` file and returns the path.
seurat_to_bpcells <- function(seurat_object) {

    seurat_object <- load_seurat(seurat_object)
    sample_id <- Seurat::Project(seurat_object)
    metadata <- seurat_object[[]] |> 
        dplyr::select(
            dplyr::any_of(
            c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_RAW", "nFeature_RAW",
                "nFeature_Diff", "nCount_Diff", "percent_mt", "percent_rb", "percent_hb",
                "percent_pl", "percent_xist", "percent_chrY")))

    bp_dir <- glue::glue("{analysis_cache}/bpcells_out/{sample_id}")
    message("Seurat object for ", sample_id, " is ", format(object.size(seurat_object), units = "Mb", digits = 2), " in memory")

    # Convert counts to BPCells / load the BPCells folder
    if (!dir.exists(bp_dir)) {
        BPCells::write_matrix_dir(mat = seurat_object[["RNA"]]$counts, dir = bp_dir, compress = TRUE)
    }
    seurat_object <- Seurat::CreateSeuratObject(counts = BPCells::open_matrix_dir(dir = bp_dir), meta.data = metadata, project = sample_id)
    message("BPCells object for ", sample_id, " is ", format(object.size(seurat_object), units = "Gb", digits = 2), " in memory")

    seurat_object_path <- glue::glue("{analysis_cache}/bpcells_out/{sample_id}.qs")

    qs::qsave(seurat_object, file = seurat_object_path)

    return(seurat_object_path)
}
