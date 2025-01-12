#' Cluster cells using Seurat package
#'
#' This function clusters cells using the Seurat package. It normalizes the data,
#' finds variable features, scales the data, runs PCA, finds neighbours, then
#' clusters the cells across the range of resolutions set by 'res'. The
#' column 'seurat_clusters' is set according to the first resolution in the
#' series with a stable ARI value, as indicated in the ARI plot (if 
#' 'plot_path' is not FALSE). A metadata table with cell names and clusters
#' is saved in the 'analysis_cache' directory and the path returned.
#' 
#' @param seurat_object A Seurat object (or its path) to be clustered.
#' @param res A vector of resolution values to use for clustering.
#' @param dims The number of dimensions to use for clustering.
#' @param regress_out Vector of column names for variables to regress out.
#' @param do.plot Logical. Save the ARI plot. Default is TRUE.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return The path to a metadata table with cell names and clusters in csv format.
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters
#' @importFrom dplyr select matches
#' @importFrom tibble rownames_to_column as_tibble
#'
#' @examples
#' if (interactive()) {
#' data("pbmc_small")
#' seurat_object <- CreateSeuratObject(counts = pbmc_small)
#' seurat_cluster(seurat_object)
#' }
#' 
#' @export
seurat_cluster_ari <- function(
    seurat_object,
    assay = "RNA",
    res = seq(0.1, 1.5, by = 0.05),
    reduction = "pca",
    dims = 30,
    vars_to_regress = c("percent_mt", "percent_rb"),
    do.plot = TRUE, ...) {
    set.seed(42)
    message("This is seurat_cluster_ari()")

    options(future.globals.maxSize = 100 * 1024^3)
    future::plan("multisession", workers = 6)

    seurat_object <- load_seurat(seurat_object)
    Seurat::DefaultAssay(seurat_object) <- assay

    # Run the standard workflow for visualization and clustering -------------------
    # https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    seurat_object <- seurat_object |>
        Seurat::SCTransform(assay = assay, vars.to.regress = vars_to_regress) |>
        # Seurat::NormalizeData() |>
        # Seurat::FindVariableFeatures() |>
        # Seurat::ScaleData(vars.to.regress = vars_to_regress) |>
        Seurat::RunPCA() |>
        Seurat::FindNeighbors(reduction = reduction, dims = 1:dims) |>
        Seurat::FindClusters(resolution = res, method = "igraph")

    clusters_csv <- seurat_object[[]] |>
        tibble::rownames_to_column(var = "cell") |>
        select("cell", "seurat_clusters", contains("SCT"), -c("nCount_SCT", "nFeature_SCT")) |>
        tibble::as_tibble()
    
    resolution_cols <- grep("SCT_snn_res", names(clusters_csv), value = TRUE)
    resolutions <- gsub("SCT_snn_res\\.", "", resolution_cols) |> as.numeric()

    ari_values <- numeric(length(resolutions) - 1)

    for (i in seq_along(resolutions)[-length(resolutions)]) {
        col1 <- resolution_cols[i]
        col2 <- resolution_cols[i + 1]
        ari_values[i] <- mclust::adjustedRandIndex(as.integer(clusters_csv[[col1]]), as.integer(clusters_csv[[col2]]))
    }
    # Identify the first resolution with a significant change
    stable_res_index <- which(diff(ari_values) > 0.05)[1] # Adjust the condition based on your criteria
    if (is.na(stable_res_index) || length(stable_res_index) == 0) {
        stable_res_index <- which.max(ari_values) # Choose max ARI if no significant change
    }
    stable_res <- resolutions[stable_res_index]

    sample_id <- seurat_object@project.name
    # if TRUE, plot ARI values
    if ((do.plot)) {
        plot_path <- glue::glue("{analysis_cache}/clusters_out/{sample_id}_{assay}_ari_res{stable_res}.png")
        dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)
        ari_data <- data.frame(Resolution = resolutions[-length(resolutions)], ARI = ari_values)
        ggplot(ari_data, aes(x = Resolution, y = ARI)) +
        geom_line() +
        geom_point() +
        geom_point(data = subset(ari_data, Resolution == stable_res), aes(x = Resolution, y = ARI), color = "red", size = 3) +
        geom_text(data = subset(ari_data, Resolution == stable_res), aes(x = Resolution, y = ARI, label = sprintf("Res: %.2f", stable_res)), vjust = -1.5, color = "red") +
        labs(title = glue::glue("Adjusted Rand Index Across Resolutions {sub('_ari.png', '', plot_path|> basename())}"), x = "Resolution", y = "ARI")
        ggsave(plot_path, width = 8, height = 6)
    }

    stable_resolution <- glue::glue("SCT_snn_res.{stable_res}")

    clusters_csv_path = glue::glue("{analysis_cache}/clusters_out/{sample_id}_{assay}_clusters.csv")
    dir.create(dirname(clusters_csv_path), showWarnings = FALSE, recursive = TRUE)

    clusters_csv$seurat_clusters <- clusters_csv[[stable_resolution]]
    clusters_csv |> write.csv(clusters_csv_path, row.names = FALSE, quote = FALSE)

    return(clusters_csv_path)
}
