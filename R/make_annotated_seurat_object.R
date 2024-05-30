#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param seurat_object
#' @param cell_metadata
#' @param sample_metadata
#' @return
#' @author Denis O'Meally
#' @export
make_annotated_seurat_object <- function(seurat_object,
                                         cell_metadata,
                                         sample_metadata) {

    # if (is.null(seurat_object)) seurat_object <- tar_read(integrated_seurat_sketch_750)
    # if (is.null(cell_metadata)) cell_metadata <- tar_read(aggregated_cell_annot_csv)
    # if (is.null(sample_metadata)) sample_metadata <- tar_read(pancdb_metadata_gt)

#    seurat_object <- targets::tar_read(integrated_seurat_sketch_750)
#    cell_metadata <- targets::tar_read(aggregated_cell_annot_csv)
#    sample_metadata <- targets::tar_read(pancdb_metadata_gt)

    seurat_object <- load_seurat(seurat_object)
    cell_metadata <- cell_metadata |>
        dplyr::select(c(cell, cell_cycle = "Phase", tosti_cell_type = "pruned.labels",
        cell_type, hpap_celltype, hpap_clusters = "cluster_labels", seurat_clusters = "manual_clusters",
        xbp1u_psi = xbp1_psi, gene_ratio_INS, gene_ratio_XBP1, gene_ratio_GAPDH)) |>
        tibble::column_to_rownames(var = "cell")

    sample_metadata <- sample_metadata |>
        select(c(
            sample_name = donor_id, sample_sex, sample_age, sample_ethnicity, diabetes_status,
            batch, protected, rs3842752_consensus, rs3842753_consensus, rs689_consensus,
            sample_xbp1u_psi = xbp1u_psi))

    seurat_object <- Seurat::AddMetaData(
        object = seurat_object,
        metadata = cell_metadata
    )

    # bug in scCustomize::Add_Sample_Meta, requires "meta_data" as df name???
    meta_data <- sample_metadata
    
    seurat_object <- scCustomize::Add_Sample_Meta(
        seurat_object = seurat_object,
        meta_data = meta_data,
#        meta_data = sample_metadata,
        join_by_seurat = "orig.ident",
        join_by_meta = "sample_name",
        overwrite = TRUE,
        na_ok = TRUE
    )

    seurat_object_path <- glue::glue("{analysis_cache}/data/seurat_object_annotated.qs")
    qs::qsave(seurat_object, seurat_object_path)
    message("Saved seurat object")

    return(seurat_object_path)
}
