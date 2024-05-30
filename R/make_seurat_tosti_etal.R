#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @return
#' @author Denis O'Meally
#' @export
make_seurat_tosti_etal <- function(integrate_across_donors = FALSE) {
    require(Seurat)
    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    base_url <- "http://singlecell.charite.de/cellbrowser/pancreas/Adult_Pancreas/"
    files_to_download <- c("exprMatrix.tsv.gz", "meta.tsv", "Seurat_umap.coords.tsv.gz")

    download_tosti_etal_files <- lapply(files_to_download, function(file) {
        dest_file_path <- file.path(analysis_cache, "data/tosti_etal", file)
        dir.create(dirname(dest_file_path), showWarnings = FALSE, recursive = TRUE)

        # Check if exists before downloading
        if (!file.exists(dest_file_path)) {
            url <- paste0(base_url, file)
            download.file(url, destfile = dest_file_path, mode = "wb")
        }
        return(dest_file_path) # Return the destination path for each file
    }) |> unlist()

    # https://cellbrowser.readthedocs.io/en/master/load.html
    meta <- read.table(download_tosti_etal_files[2], header = T, sep = "\t", as.is = T, row.names = 1)
    umap <- data.table::fread(download_tosti_etal_files[3])
    umap_coords <- as.matrix(umap[, -1])
    rownames(umap_coords) <- umap$V1
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

    # Read expression matrix and convert to sparse matrix
    exprMatrix <- data.table::fread(download_tosti_etal_files[1], header = TRUE, sep = "\t", data.table = FALSE)
    rownames(exprMatrix) <- gsub(".+[|]", "", exprMatrix[, 1])
    exprMatrix <- exprMatrix[, -1]
    sparse_mat <- as(Matrix::Matrix(as.matrix(exprMatrix), sparse = TRUE), "dgCMatrix")
    rm(exprMatrix)
    # Create Seurat object
    seurat_object <- Seurat::CreateSeuratObject(counts = sparse_mat, project = "tosti_etal", meta.data = meta)
    # Add the UMAP coordinates to the Seurat object
    seurat_object[["umap_tosti_etal"]] <- Seurat::CreateDimReducObject(embeddings = umap_coords, key = "UMAP_")
    # make a UMAP plot
    plot_umap_tosti_etal <- Seurat::DimPlot(seurat_object, group.by = "Cluster", reduction = "umap_tosti_etal", cols = cell_type_palette, label = TRUE, repel = TRUE, label.size = 4, shuffle = TRUE) & NoLegend() & labs(title = "UMAP Tosti et al.") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    ggsave(glue::glue("{analysis_cache}/data/tosti_etal/tosti_etal_UMAP_coords.png", plot = plot_umap_tosti_etal, width = 10, height = 10, dpi = 300))

    # Run PCA
    seurat_object <- Seurat::SCTransform(seurat_object) |> RunPCA()

    if (integrate_across_donors) {
        # Integrate across donors
        seurat_object[["SCT"]] <- split(seurat_object[["SCT"]], f = seurat_object$patient_ID)
        seurat_object <- IntegrateLayers(
            object = seurat_object, method = HarmonyIntegration,
            orig.reduction = "pca", new.reduction = "harmony",
            verbose = TRUE
        )
        seurat_object <- JoinLayers(seurat_object)

        # Add a UMAP model to the Seurat object
        seurat_object <- seurat_object |>
            Seurat::FindNeighbors(dims = 1:30) |>
            Seurat::FindClusters(cluster.name = "harmony_clusters") |>
            Seurat::RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", return.model = TRUE, repulsion.strength = 5)
        plot_umap_ref <- DimPlot(seurat_object, group.by = "Cluster", cols = cell_type_palette, label = TRUE, label.size = 6, repel = TRUE, reduction = "umap_harmony", shuffle = TRUE) & NoLegend() & labs(title = "UMAP reference") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
        plot_umap_ref_path <- glue::glue("{analysis_cache}/data/tosti_etal/tosti_etal_UMAP_harmony.png")
        seurat_object_path <- (glue::glue("{analysis_cache}/data/tosti_etal/tosti_etal_harmony.qs"))
    } else {
        # Add a UMAP model to the Seurat object
        seurat_object <- seurat_object |>
            Seurat::FindNeighbors(dims = 1:30) |>
            Seurat::FindClusters() |>
            Seurat::RunUMAP(dims = 1:30, reduction = "pca", return.model = TRUE, reduction.name = "umap_sct")
        plot_umap_ref <- DimPlot(seurat_object, group.by = "Cluster", cols = cell_type_palette, label = TRUE, label.size = 6, repel = TRUE, reduction = "umap_sct", shuffle = TRUE) & NoLegend() & labs(title = "UMAP reference") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
        plot_umap_ref_path <- glue::glue("{analysis_cache}/data/tosti_etal/tosti_etal_UMAP_sct.png")
        seurat_object_path <- (glue::glue("{analysis_cache}/data/tosti_etal/tosti_etal_sct.qs"))
    }


    ridge_plot <- RidgePlot(seurat_object, features = cell_type_markers, group.by = "Cluster", ncol = 3, log = TRUE, y.max = 200) + ggplot2::theme(aspect.ratio = 1)
    ggsave(plot = ridge_plot, filename = glue::glue("{analysis_cache}/data/tosti_etal/tosti_etal_RidgePlot.png"), width = 20, height = 20, dpi = 300)
    ggsave(plot = plot_umap_ref, filename = plot_umap_ref_path, width = 10, height = 10, dpi = 300)


    qs::qsave(seurat_object, seurat_object_path)

    return(seurat_object_path)
}
