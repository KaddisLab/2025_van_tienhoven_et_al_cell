#' Perform CHOIR Analysis on a Seurat Object
#'
#' This function integrates CHOIR clustering analysis with a Seurat object
#' and optionally performs batch correction using Harmony if batch labels are provided.
#' It adds sample metadata, performs the CHOIR analysis, and saves the results.
#'
#' @param seurat_object Path to an RDS file or a Seurat object in the environment.
#' @param assay The name of the assay to use for the CHOIR analysis.
#' @param batch Optional; a vector of batch labels corresponding to the `batch` column in the Seurat object metadata.
#' @param meta_data A data frame containing sample metadata to add to the Seurat object. Must include 'orig.ident' and 'sample_name' columns.
#'
#' @return A string containing the path to the saved CHOIR CSV results.
#'
#' @examples
#' # Assuming that 'seurat' is your Seurat object, 'meta' is the metadata DataFrame
#' # containing 'orig.ident' and 'sample_name', and 'batch_labels' is your batch information:
#' choir_results_path <- seurat_CHOIR(
#'     seurat_object = seurat,
#'     assay = "RNA",
#'     batch = batch_labels,
#'     meta_data = meta
#' )
#'
#' @export
#'
#' @importFrom CHOIR CHOIR
#' @importFrom Seurat Project
#' @importFrom glue glue
#' @importFrom qs qsave
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select contains
seurat_CHOIR <- function(seurat_object, assay, batch = NULL, meta_data = NULL) {
    # https://www.choirclustering.com/articles/CHOIR.html
    require(CHOIR)

    if (nzchar(Sys.getenv("SLURM_JOB_ID"))) {
        hprcc::init_multisession()
        cpus <- hprcc::slurm_allocation()$CPUs
    } else {
        cpus <- future::availableCores()
    }

    message("CPUs allocated: ", cpus)

    seurat_object <- load_seurat(seurat_object)

    if (!is.null(batch)) {
        if (is.null(meta_data)) stop("Batch correction requires sample metadata.")
        batch_method <- "Harmony"
        seurat_object <- scCustomize::Add_Sample_Meta(seurat_object = seurat_object, meta_data = meta_data, "orig.ident", "sample_name")
        message("Sample metadata added to Seurat object; batch correction will be performed.")
    } else {
        batch_method <- "none"
        message("No batch correction will be performed.")
    }

    message("Starting CHOIR analysis...")
        seurat_object <- CHOIR::CHOIR(
            object = seurat_object,
            use_assay = assay,
            batch_correction_method = batch_method,
            batch_labels = batch,
            n_cores = cpus, 
            random_seed = 42,
            verbose = TRUE
        )
        message("CHOIR analysis completed.")

        sample_id <- Seurat::Project(seurat_object)
        choir_csv_path <- glue::glue("{analysis_cache}/CHOIR_out/CHOIR_{assay}_{sample_id}.csv")
        choir_seurat_path <- glue::glue("{analysis_cache}/CHOIR_out/CHOIR_{assay}_{sample_id}.qs")
        choir_csv <- seurat_object[[]] |>
            tibble::rownames_to_column(var = "cell") |>
            select("cell", contains("CHOIR"))

        dir.create(dirname(choir_csv_path), showWarnings = FALSE, recursive = TRUE)
        write.csv(choir_csv, choir_csv_path, row.names = FALSE, quote = FALSE)
        qs::qsave(seurat_object, file = choir_seurat_path)
        message("Results saved.")
    
    #save profile reportt to file
    return(choir_csv_path)
}


# seurat_CHOIR <- function(input_seurat_object, assay, batch = NULL, sample_metadata) {
#     require(CHOIR)

#     options(parallelly.availableCores.methods = "Slurm")
#     hprcc::init_multisession()

#     cpus <- hprcc::slurm_allocation()$CPUs
#     message("CPUs allocated: ", cpus)

#     if (!is.null(batch)) {
#         batch_method <- "Harmony"
#     } else {
#         batch_method <- NULL
#     }

#     message("Loading Seurat object...")
#     loaded_seurat_object <- load_seurat(input_seurat_object)
#     if (is.null(loaded_seurat_object)) {
#         stop("Failed to load Seurat object.")
#     }
#     message("Seurat object loaded successfully.")

#     message("Filtering and preparing sample metadata...")
#     # Check for required columns in sample_metadata
#     required_columns <- c("reagent_kit", "donor_id", "sample_sex", "sample_age", "sample_ethnicity", "ab_positive", "diabetes_status")
#     if (!all(required_columns %in% names(sample_metadata))) {
#         stop("Sample metadata is missing one or more required columns.")
#     }
#     meta_data <- sample_metadata |>
#         dplyr::filter(dplyr::str_detect(reagent_kit, "10X")) |>
#         dplyr::select(c(sample_id = "donor_id", "sample_sex", "sample_age", "sample_ethnicity", "ab_positive", "diabetes_status", "reagent_kit")) |>
#         dplyr::mutate(batch = as.integer(as.factor(reagent_kit))) |>
#         dplyr::select(-reagent_kit)
#     message("Sample metadata prepared.")

#     message("Adding sample metadata to Seurat object...")
#     seurat_with_metadata <- scCustomize::Add_Sample_Meta(seurat_object = loaded_seurat_object, meta_data = meta_data, "orig.ident", "sample_id")
#     if (is.null(seurat_with_metadata)) {
#         stop("Failed to add metadata to Seurat object.")
#     }
#     message("Sample metadata added successfully.")

#     message("Starting CHOIR analysis...")
#     seurat_after_CHOIR <- CHOIR::CHOIR(
#         object = seurat_with_metadata,
#         use_assay = assay,
#         batch_correction_method = batch_method,
#         batch_labels = batch,
#         n_cores = cpus,
#         random_seed = 42,
#         verbose = TRUE
#     )
#     if (is.null(seurat_after_CHOIR)) {
#         stop("CHOIR analysis failed.")
#     }
#     message("CHOIR analysis completed successfully.")

#     message("Projecting data for output...")
#     sample_id <- Seurat::Project(seurat_after_CHOIR)
#     choir_csv_path <- glue::glue("{analysis_cache}/CHOIR_out/CHOIR_{assay}_{sample_id}.csv")
#     choir_seurat_path <- glue::glue("{analysis_cache}/CHOIR_out/CHOIR_{assay}_{sample_id}.qs")

#     choir_csv <- seurat_after_CHOIR[[]] |>
#         tibble::rownames_to_column(var = "cell") |>
#         dplyr::select("cell", dplyr::contains("CHOIR"))

#     message("Saving results...")
#     dir.create(dirname(choir_csv_path), showWarnings = FALSE, recursive = TRUE)
#     write.csv(choir_csv, choir_csv_path, row.names = FALSE, quote = FALSE)
#     qs::qsave(seurat_after_CHOIR, file = choir_seurat_path)
#     message("Results saved at: ", choir_csv_path)

#     # Manual garbage collection to free up memory
#     gc()

#     return(choir_csv_path)
# }
