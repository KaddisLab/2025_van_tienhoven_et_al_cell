#' Calculate Gene Signature Scores Using UCell
#'
#' Calculates gene signature scores for cells in a Seurat object using the UCell algorithm.
#' Saves the results to a CSV file and returns the scores data frame.
#'
#' @param seurat_object A Seurat object or path to a Seurat object file
#' @param features List of gene signatures, where each element is a character vector
#'        of gene symbols comprising that signature
#'
#' @return A data frame containing:
#'   - cell: Cell barcode identifier
#'   - [signature]_UCell: UCell score for each signature
#'
#' @details
#' Uses the UCell algorithm to calculate signature scores, which are saved to
#' "{analysis_cache}/sig_scores_out/{project_name}_sig_scores.csv"
#'
#' @import UCell
#' @importFrom BiocParallel MulticoreParam
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select contains
#' @export
seurat_sig_scores <- function(seurat_object, features) {
    message("This is seurat_sig_scores()")
    seurat_object <- load_seurat(seurat_object)

    num_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = "1"))

    # Get the project name
    project_name <- Seurat::Project(seurat_object)

    obj.score <- UCell::AddModuleScore_UCell(seurat_object,
        features = features,
        assay = "RNA",
        BPPARAM = BiocParallel::MulticoreParam(workers = num_cores),
        ncores = cpus
    )

    scores <- obj.score[[]] |>
        rownames_to_column(var = "cell") |>
        select(cell, contains("_UCell"))
    # Save the scores
    score_path <- glue::glue("{analysis_cache}/sig_scores_out/{project_name}_sig_scores.csv")
    save_results(scores, score_path)

    return(scores)
}
