#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @return
#' @author Denis O'Meally
#' @export
make_seurat_tosti_etal <- function() {

    base_url <- "http://singlecell.charite.de/cellbrowser/pancreas/Adult_Pancreas/"
    files_to_download <- c("exprMatrix.tsv.gz", "meta.tsv", "Seurat_umap.coords.tsv.gz")

    download_tosti_etal_files <- lapply(files_to_download, function(file) {
    dest_file_path <- file.path(analysis_cache, "data/tosti_etal", file)
    dir.create(dirname(dest_file_path), showWarnings = FALSE, recursive = TRUE)
    
    # Check if this specific file does not exist before downloading
    if (!file.exists(dest_file_path)) {
        url <- paste0(base_url, file)
        download.file(url, destfile = dest_file_path, mode = "wb")
    }
    return(dest_file_path) # Return the destination path for each file
    })

    # Convert the list of paths to a character vector
    download_tosti_etal_files <- unlist(download_tosti_etal_files)

    #https://cellbrowser.readthedocs.io/en/master/load.html
    meta <- read.table(download_tosti_etal_files[2], header=T, sep="\t", as.is=T, row.names=1)
    umap <- data.table::fread(download_tosti_etal_files[3])
    umap_coords <- as.matrix(umap[, -1])
    rownames(umap_coords) <- umap$V1
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

    # Read expression matrix and convert to sparse matrix
    exprMatrix <- data.table::fread(download_tosti_etal_files[1], header = TRUE, sep = "\t", data.table = FALSE)
    rownames(exprMatrix) <- gsub(".+[|]", "", exprMatrix[, 1])
    exprMatrix <- exprMatrix[, -1]
    sparse_mat <- as(Matrix::Matrix(as.matrix(exprMatrix), sparse = TRUE), "dgCMatrix")

    # Create Seurat object
    so <- Seurat::CreateSeuratObject(counts = sparse_mat, project = "tosti_etal", meta.data = meta) 
    # Normalize the data
    so <- Seurat::NormalizeData(so)
    # Add the UMAP coordinates to the Seurat object
    so[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_coords, key = "UMAP_")

    # make a UMAP plot
    plot <- Seurat::DimPlot(so, group.by = "Cluster", reduction = "umap")
    ggsave(glue::glue("{analysis_cache}/data/seurat_tosti_etal.png", plot = umap_plot, width = 12, height = 8, dpi = 300))

    so_path <- (glue::glue("{analysis_cache}/data/seurat_tosti_etal.qs"))

    qs::qsave(so, so_path)
    
    return(so_path)
}
