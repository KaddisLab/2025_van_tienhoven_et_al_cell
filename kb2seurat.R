library(targets) # Hypothetical package for handling Kallisto BUS files, replace with actual if available

tar_load(kb_count_tcc_alltech)

bus_file <- kb_count_tcc_alltech[1]

bus_file <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/kb_out/HPAP-113/output.bus"

# Function to load spliced and unspliced counts from a Kallisto BUS file
library(Seurat)
library(dplyr)
library(Matrix)
library(BUSpaRse)
library(Seurat)
library(dplyr)
library(Matrix)
library(BUSpaRse)

library(Seurat)
library(Matrix)
library(BUSpaRse)
library(SeuratWrappers)
library(velocyto.R)
library(DropletUtils)


###########################

# Adapted read_count_output function
read_count_output <- function(dir, name, tcc = FALSE) {
    dir <- normalizePath(dir, mustWork = TRUE)
    m <- Matrix::readMM(paste0(dir, "/", name, ".mtx"))
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    ge <- if (tcc) ".ec.txt" else ".genes.txt"
    con_genes <- file(paste0(dir, "/cells_x_genes", ge))
    con_bcs <- file(paste0(dir, "/cells_x_genes.barcodes.txt"))
    genes <- readLines(con_genes)
    barcodes <- readLines(con_bcs)
    colnames(m) <- barcodes
    rownames(m) <- genes
    close(con_genes)
    close(con_bcs)
    return(m)
}
# Adapted read_velocity_output function
read_velocity_output <- function(spliced_dir, unspliced_dir, spliced_name, unspliced_name) {
    spliced <- read_count_output(spliced_dir, spliced_name, FALSE)
    unspliced <- read_count_output(unspliced_dir, unspliced_name, FALSE)
    list(spliced = spliced, unspliced = unspliced)
}

# Main function to load spliced and unspliced counts into a Seurat object
import_kallisto_bustools <- function(bus_file) {
    # Define base path and necessary subdirectories
    base_path <- dirname(bus_file)
    spliced_dir <- file.path(base_path, "counts_unfiltered")
    unspliced_dir <- file.path(base_path, "counts_unfiltered")
    spliced_name <- "cells_x_genes.mature"
    unspliced_name <- "cells_x_genes.nascent"

    # Read velocity output using adapted BUSpaRse functions
    velocity_data <- read_velocity_output(
        spliced_dir = spliced_dir,
        unspliced_dir = unspliced_dir,
        spliced_name = spliced_name,
        unspliced_name = unspliced_name
    )

    spliced_counts <- velocity_data$spliced
    unspliced_counts <- velocity_data$unspliced

    # Create Seurat object with spliced counts
    seurat_object <- Seurat::CreateSeuratObject(counts = spliced_counts, assay = "spliced")

    # Add unspliced counts as a second assay
    seurat_object[["unspliced"]] <- Seurat::CreateAssayObject(counts = unspliced_counts)

    return(seurat_object)
}

# x <- import_kallisto_bustools(bus_file)
# qs::qsave(x, file = "HPAP113_velocity_seurat_object.qs")
obj<-qs::qread("HPAP113_velocity_seurat_object.qs")
spliced <- GetAssayData(obj, layer = "counts", assay = "spliced")
unspliced <- GetAssayData(obj, layer = "counts", assay = "unspliced")
sum(unspliced@x) / (sum(unspliced@x) + sum(spliced@x))
dim(spliced)
dim(unspliced)
bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)
tot_count <- Matrix::colSums(spliced)
summary(tot_count)
#' Knee plot for filtering empty droplets
#'
#' Visualizes the inflection point to filter empty droplets. This function plots
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions.
#'
#' @param bc_ranks A named list of output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
#' @importFrom tibble tibble
#' @importFrom purrr map map_dbl
#' @importFrom dplyr distinct
#' @importFrom ggplot2 geom_line geom_hline geom_vline scale_x_log10 scale_y_log10
#' @importFrom tidyr unnest
#' @export
knee_plot <- function(bc_ranks) {
    # purrr pluck shorthand doesn't work on S4Vector DataFrame
    knee_plt <- tibble(
        rank = map(bc_ranks, ~ .x[["rank"]]),
        total = map(bc_ranks, ~ .x[["total"]]),
        dataset = names(bc_ranks)
    ) %>%
        unnest(cols = c(rank, total)) %>%
        distinct() %>%
        dplyr::filter(total > 0)
    annot <- tibble(
        inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
        rank_cutoff = map_dbl(
            bc_ranks,
            ~ max(.x$rank[.x$total >
                metadata(.x)[["inflection"]]])
        ),
        dataset = names(bc_ranks)
    )
    p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
        geom_line() +
        geom_hline(aes(yintercept = inflection, color = dataset),
            data = annot, linetype = 2
        ) +
        geom_vline(aes(xintercept = rank_cutoff, color = dataset),
            data = annot, linetype = 2
        ) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x = "Rank", y = "Total UMIs")
    return(p)
}
knee_plot(list(spliced = bc_rank, unspliced = bc_uns)) +
    coord_flip()

bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]
# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 0]
spliced_filtered <- spliced[genes_use, bcs_use]
unspliced_filtered <- unspliced[genes_use, bcs_use]
dim(spliced_filtered)
dim(unspliced_filtered)
rownames(spliced_filtered) <- str_remove(rownames(spliced_filtered), "\\.\\d+")
rownames(unspliced_filtered) <- str_remove(rownames(unspliced_filtered), "\\.\\d+")

seu <- CreateSeuratObject(spliced_filtered, assay = "sf") %>%
    SCTransform(assay = "sf", new.assay.name = "spliced")
seu <- RunPCA(seu, assay = "spliced", verbose = FALSE)
seu[["uf"]] <- CreateAssayObject(unspliced_filtered)
seu <- SCTransform(seu, assay = "uf", new.assay.name = "unspliced")
seu <- SeuratWrappers::RunVelocity(seu, ncores = 6, reduction = "pca", verbose = FALSE)

DefaultAssay(seu) <- "spliced"
seu<-seu|>RunPCA()|>FindNeighbors()|>FindClusters(res = 0.2)
seu$seurat_clusters|>table()

cols_use <- c("nCount_sf", "nFeature_sf", "nCount_uf", "nFeature_uf")
VlnPlot(seu, cols_use, pt.size = 0.1, ncol = 1, group.by = "seurat_clusters")
ggsave("plot.png")

# Helper functions for ggpairs
log10_diagonal <- function(data, mapping, ...) {
  GGally::ggally_densityDiag(data, mapping, ...) + scale_x_log10()
}
log10_points <- function(data, mapping, ...) {
  GGally::ggally_points(data, mapping, ...) + scale_x_log10() + scale_y_log10()
}
GGally::ggpairs(seu@meta.data, columns = cols_use,
        upper = list(continuous = "cor"),
        diag = list(continuous = log10_diagonal),
        lower = list(continuous = GGally::wrap(log10_points, alpha = 0.1, size=0.3)),
        progress = FALSE)
ggsave("plot.png")

DefaultAssay(seu) <- "spliced"
seu <- RunPCA(seu, verbose = FALSE, npcs = 70)
ElbowPlot(seu, ndims = 70) + theme_classic()
ggsave("plot.png")

seu <- RunUMAP(seu, dims = 1:50)
DimPlot(seu, label = TRUE) + NoLegend()
ggsave("plot.png")


#INS
VlnPlot(seu, "ENSG00000254647", pt.size = 0.1, ncol = 1, group.by = "seurat_clusters")
ggsave("plot.png")

# XBP1
VlnPlot(seu, "ENSG00000100219", pt.size = 0.1, ncol = 1, group.by = "seurat_clusters")
ggsave("plot.png")
