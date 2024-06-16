#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param seurat_object_lognorm_annotated
#' @return
#' @author Denis O'Meally
#' @export
make_beta_cell_object <- function(seurat_object) {
    require(tidyseurat)

    total_endocrine_cells_per_donor <- seurat_object |>
        dplyr::filter(cell_type %in% c("Alpha", "Beta", "Alpha+Beta", "Delta", "Gamma", "Epsilon", "Gamma+Epsilon")) |>
        dplyr::count(orig.ident) |>
        dplyr::rename(donor = orig.ident, total_count = n)


    seurat_object <- seurat_object |>
        dplyr::filter(diabetes_status == "NODM") |>
        dplyr::filter(cell_type %in% c("Beta", "Alpha+Beta")) |>
        Seurat::SCTransform(vst.flavor = "v2", vars.to.regress = c("percent_mt", "percent_rb"), verbose = FALSE) |>
        Seurat::RunPCA(verbose = FALSE) |>
        harmony::RunHarmony(group.by.vars = c("tissue_source", "reagent_kit", "orig.ident"), verbose = FALSE) |>
        Seurat::FindNeighbors(reduction = "harmony", verbose = FALSE) |>
        Seurat::FindClusters(resolution = 0.1, method = "igraph", verbose = FALSE) |>
        Seurat::RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE)

    beta_cell_list <- list(
        seurat_object = seurat_object,
        total_endocrine_cells_per_donor = total_endocrine_cells_per_donor
    )
    return(beta_cell_list)
}
