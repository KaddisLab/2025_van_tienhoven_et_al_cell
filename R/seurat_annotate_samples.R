#' Add Sample-level Metadata to a Seurat Object
#'
#' Annotates a Seurat object with sample-level metadata including clinical,
#' technical, and genetic information. Joins metadata using sample identifiers
#' from the orig.ident field.
#'
#' @param seurat_object A Seurat object or path to a Seurat object file
#' @param sample_metadata Data frame containing sample metadata with required columns:
#'        donor_id, sample_sex, sample_age, sample_ethnicity, diabetes_status,
#'        number_abs, tissue_source, reagent_kit, technology, protected,
#'        rs3842752_consensus, rs3842753_consensus, rs689_consensus, xbp1u_psi
#'
#' @return A Seurat object with updated metadata containing sample information
#'
#' @details
#' Adds the following metadata columns to the Seurat object:
#' - sample_name (from donor_id)
#' - Clinical: sample_sex, sample_age, sample_ethnicity, diabetes_status
#' - Technical: tissue_source, reagent_kit, technology
#' - Genetic: protected status, SNP consensus genotypes
#' - Molecular: XBP1 PSI values
#'
#' @import Seurat
#' @importFrom scCustomize Add_Sample_Meta
#' @importFrom dplyr select
#' @author Denis O'Meally
#' @export
seurat_annotate_samples <- function(seurat_object, sample_metadata) {
    # if (is.null(seurat_object)) seurat_object <- tar_read(integrated_seurat_sketch_750)
    # if (is.null(cell_metadata)) cell_metadata <- tar_read(aggregated_cell_annot_csv)
    # if (is.null(sample_metadata)) sample_metadata <- tar_read(pancdb_metadata_gt)

    #    seurat_object <- targets::tar_read(integrated_seurat_sketch_750)
    #    cell_metadata <- targets::tar_read(aggregated_cell_annot_csv)
    #    sample_metadata <- targets::tar_read(pancdb_metadata_agg)

    seurat_object <- load_seurat(seurat_object)
    assay <- Seurat::DefaultAssay(seurat_object)

    sample_metadata <- sample_metadata |>
        select(c(
            sample_name = donor_id, sample_sex, sample_age, sample_ethnicity, diabetes_status, number_abs,
            tissue_source, reagent_kit, technology, protected, rs3842752_consensus, rs3842753_consensus, rs689_consensus,
            sample_xbp1u_psi = xbp1u_psi
        ))

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

    return(seurat_object)
}
