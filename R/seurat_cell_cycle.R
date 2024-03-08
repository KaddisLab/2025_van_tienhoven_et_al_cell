#' Annotate human or mouse cell cycle
#'
#' Uses the function `Seurat::CellCycleScoring()` to annotate the cells with their cell cycle `Phase`,
#' `S.score` and `G2M.score`. The Seurat function uses an internal list of human gene symbols from
#' Tirosh et al 2016. This function uses the `gprofiler2::gorth()` function to convert the human
#' gene symbols to mouse gene symbols when `mmu=TRUE`.
#'
#' @title Annotate human or mouse cell cycle
#' @description This function annotates the cell cycle phases of cells in a Seurat object.
#' @param seurat_object a Seurat object or a path to a seurat object saved as "qs" format
#' @param mmu Logical indicating whether to convert human gene symbols to mouse gene symbols.
#' @import Seurat
#' @import gprofiler2
#' @return A tibble with columns `cell`, `S.Score`, `G2M.Score`, and `Phase`.
#' @references
#' - Tirosh, I., et al. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 352(6282), 189-196.
#' - Seurat package documentation: https://rdrr.io/cran/Seurat/man/CellCycleScoring.html
#' - GitHub issue discussing the use of `gprofiler2::gorth()`: https://github.com/satijalab/seurat/issues/2493
#' @examples
#' # Annotate cell cycle for a Seurat object
#' annotated_object <- seurat_annotate_cell_cycle(seurat_object)
#' @export
seurat_cell_cycle <- function(seurat_object, mmu = FALSE) {

    # Check if the input is a Seurat object, else read it from qs file
    if (!inherits(seurat_object, "Seurat")) {
        seurat_object <- qs::qread(seurat_object)
    }
    if (!isNormalized(seurat_object)) {
        seurat_object <- Seurat::NormalizeData(seurat_object)
    }

    if (isTRUE(mmu)) {
        message("Converting human gene symbols to mouse gene symbols")
        s_features = gprofiler2::gorth(Seurat::cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
        g2m_features = gprofiler2::gorth(Seurat::cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    } else {
        s_features = Seurat::cc.genes.updated.2019$s.genes
        g2m_features = Seurat::cc.genes.updated.2019$g2m.genes
    }

    seurat_object <- Seurat::CellCycleScoring(
        object = seurat_object,
        g2m.features = g2m_features,
        s.features = s_features
    )
    cell_cycle <- seurat_object[[]] |> as.data.frame() |> tibble::rownames_to_column(var = "cell") |> select("cell", "S.Score", "G2M.Score", "Phase") |> unique() |> tibble::as_tibble()

    sample_id <- seurat_object@meta.data$orig.ident[1]
    cell_cycle_path <- glue::glue("{analysis_cache}/cell_cycle_out/{sample_id}_cell_cycle.tsv")
    dir.create(dirname(cell_cycle_path), showWarnings = FALSE, recursive = TRUE)
    write.table(cell_cycle, cell_cycle_path, sep = "\t", quote = FALSE, row.names = FALSE)
    return(cell_cycle_path)
}
