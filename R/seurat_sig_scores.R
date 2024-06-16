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


