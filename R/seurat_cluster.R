#' Cluster cells using Seurat package
#'
#' This function clusters cells using the Seurat package. It normalizes the data, finds variable features, scales the data, runs PCA, finds neighbors, and finally clusters the cells. It returns a metadata table with the cell names and their corresponding clusters.
#'
#' @param seurat_object A Seurat object (or its path) containing the single-cell RNA-seq data.
#' @param res A vector of resolution values to use for clustering.
#' @param dims The number of dimensions to use for clustering.
#'
#' @return A metadata table with the cell names and their corresponding clusters.
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters
#' @importFrom dplyr select matches
#' @importFrom tibble rownames_to_column as_tibble
#'
#' @examples
#' data("pbmc_small")
#' seurat_object <- CreateSeuratObject(counts = pbmc_small)
#' seurat_cluster(seurat_object)
#'
#' @export
seurat_cluster_ari <- function(seurat_object, res = seq(0.1, 1.5, by = 0.05), dims = 30, regress_out = NULL, remove_outliers = TRUE, ...) {
    set.seed(42)

    seurat_object <- load_seurat(seurat_object)

    if (isTRUE(remove_outliers)) {
        if ("is_outlier" %in% colnames(seurat_object[[]])) {
            seurat_object <- seurat_object[, !seurat_object[["is_outlier"]]]
        } else {
            message("No \"is_outlier\" column found in seurat_object, run seurat_flag_outliers_mad() first. Skipping outlier removal.")
        }
    }

    # Run the standard workflow for visualization and clustering -------------------
    # https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    seurat_object <- seurat_object |>
        Seurat::SCTransform(vars.to.regress = regress_out) |>
        Seurat::RunPCA() |>
        Seurat::FindNeighbors(reduction = "pca", dims = 1:dims) |>
        Seurat::FindClusters(resolution = res, method = "igraph")

    sample_id <- seurat_object@meta.data$orig.ident[1]
    clusters_csv_path <- glue::glue("{analysis_cache}/clusters_out/{sample_id}_clusters.csv")
    dir.create(dirname(clusters_csv_path), showWarnings = FALSE, recursive = TRUE)

    clusters_csv <- seurat_object[[]] |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "cell") |>
        select("cell", "seurat_clusters", contains("SCT"), -c("nCount_SCT", "nFeature_SCT")) |>
        unique() |>
        tibble::as_tibble()
    
    resolution_cols <- grep(res_prefix, names(clusters_csv), value = TRUE)
    resolutions <- gsub(paste0("^", res_prefix, "\\."), "", resolution_cols) |> as.numeric()

    ari_values <- numeric(length(resolutions) - 1)

    for (i in seq_along(resolutions)[-length(resolutions)]) {
        col1 <- resolution_cols[i]
        col2 <- resolution_cols[i + 1]

        # Ensure both columns exist to calculate ARI
        ari_values[i] <- mclust::adjustedRandIndex(as.integer(clusters_csv[[col1]]), as.integer(clusters_csv[[col2]]))
    }
    # Identify the first resolution with a significant change
    stable_res_index <- which(diff(ari_values) > 0.05)[1] # Adjust the condition based on your criteria
    if (is.na(stable_res_index) || length(stable_res_index) == 0) {
        stable_res_index <- which.max(ari_values) # Choose max ARI if no significant change
    }
    stable_res <- resolutions[stable_res_index]

    plot_path <- "{analysis_cache}/clusters_out/{sample_id}_ari_res{stable_res}.png"

    # if plot path is not null, plot ARI values
    if (!is.null(plot_path)) {
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

    clusters_csv$seurat_clusters <- clusters_csv[[stable_resolution]]

    clusters_csv |> write.csv(clusters_csv_path, row.names = FALSE, quote = FALSE)

    return(clusters_csv_path)
}
