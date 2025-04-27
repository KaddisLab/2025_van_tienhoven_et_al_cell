#' Normalize INS Expression by Housekeeping Genes
#'
#' Calculates a housekeeping gene-normalized expression value for INS (insulin) 
#' in each cell. The normalization is performed by dividing INS counts by the sum
#' of housekeeping gene counts, scaling to 10,000, and log-transforming.
#'
#' @param seurat_object A Seurat object or path to a Seurat object file. Must
#'        contain INS and housekeeping genes in the RNA assay.
#'
#' @return A tibble containing:
#'   - cell: Cell barcode identifier
#'   - INS_hk: Housekeeping-normalized INS expression value
#'
#' @details
#' The function:
#' 1. Verifies presence of housekeeping genes
#' 2. Extracts INS and housekeeping gene counts
#' 3. Calculates normalized expression: log1p((INS / sum(housekeeping)) * 10000)
#' 4. Saves results to "{analysis_cache}/INS_hknorm_out/{project_name}_INS_hknorm_per_cell.csv"
#'
#' If INS gene is not found, returns 0 for all cells.
#'
#' @import Seurat
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @export
seurat_INS_hknorm <- function(seurat_object) {
    message("This is seurat_INS_hknorm()")
    seurat_object <- load_seurat(seurat_object)

    # Get the project name
    project_name <- Seurat::Project(seurat_object)

    message("Normalising INS expression by housekeeping genes for project: ", project_name)

    missing_genes <- setdiff(housekeeping_genes, rownames(seurat_object[["RNA"]]$counts))
    if (length(missing_genes) > 0) {
        stop("The following housekeeping genes are not present in the Seurat object: ", paste(missing_genes, collapse = ", "))
    }

    # Initialize INS_counts
    INS_counts <- NULL

    # Try to extract counts for INS
    tryCatch(
        {
            INS_counts <- seurat_object[["RNA"]]$counts["INS", ]
        },
        error = function(e) {
            # If INS is not found, print a message
            message("INS gene is not present in the Seurat object. Setting INS_hk to 0 for all cells.")
        }
    )

    if (is.null(INS_counts)) {
        # If INS gene is not found, set INS_hk to 0 for all cells
        cell_names <- colnames(seurat_object[["RNA"]]$counts)
        INS_hk_tibble <- tibble(cell = cell_names, INS_hk = 0)
    } else {
        # Extract counts for housekeeping genes
        housekeeping_counts <- seurat_object[["RNA"]]$counts[housekeeping_genes, ]

        # Sum of housekeeping gene counts
        sum_housekeeping_counts <- colSums(housekeeping_counts)

        # Normalize INS counts by the sum of housekeeping gene counts,
        # scaled to a total of 10,000 and natural log transformed
        INS_normalized <- log1p((INS_counts / sum_housekeeping_counts) * 10000)

        # Replace infinite values with 0, i.e., some housekeeping genes have 0
        # counts so we can't normalize INS expression in those cells
        INS_normalized[is.infinite(INS_normalized)] <- 0

        message("INS expression normalised by housekeeping genes. Mean INS expression: ", mean(INS_normalized))

        INS_hk_tibble <- tibble(cell = names(INS_normalized), INS_hk = as.numeric(INS_normalized))
    }

    INS_hk_tibble_path <- glue::glue("{analysis_cache}/INS_hknorm_out/{project_name}_INS_hknorm_per_cell.csv")

    save_results(INS_hk_tibble, INS_hk_tibble_path)

    return(INS_hk_tibble)
}
