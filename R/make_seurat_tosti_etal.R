#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param download_tosti_etal_files
#' @return
#' @author Denis O'Meally
#' @export
make_seurat_tosti_etal <- function(download_tosti_etal_files) {

  #https://cellbrowser.readthedocs.io/en/master/load.html
    mat <- data.table::fread(download_tosti_etal_files[1])
    meta <- read.table(download_tosti_etal_files[2], header=T, sep="\t", as.is=T, row.names=1)
    umap <- data.table::fread(download_tosti_etal_files[3])
    umap_coords <- as.matrix(umap[, -1])
    rownames(umap_coords) <- umap$V1
    genes = mat[,1][[1]]
    genes = gsub(".+[|]", "", genes)
    mat = data.frame(mat[,-1], row.names=genes)
    so <- CreateSeuratObject(counts = mat, project = "adultPancreas", meta.data=meta)
    so[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", name = "umap")

    so_path <- (paste0(analysis_cache, "/data/seurat_tosti_etal_", Sys.Date(), ".qs"))
    qs::qsave(so, so_path)
    return(so_path)

}
