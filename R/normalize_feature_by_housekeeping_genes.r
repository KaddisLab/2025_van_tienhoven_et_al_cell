normalize_feature_by_housekeeping_genes <- function(input_data, feature, housekeeping_genes) {
    # Check if input_data is a Seurat object
    if (inherits(input_data, "Seurat")) {
        counts <- GetAssayData(input_data, slot = "counts")
    } else if (inherits(input_data, "dgCMatrix")) {
        counts <- input_data
    } else {
        stop("input_data must be a Seurat object or a dgCMatrix.")
    }

    # Check if target_feature and housekeeping_genes are in the dataset
    if (!(feature %in% rownames(counts))) {
        stop(paste("The target feature", feature, "is not found in the dataset."))
    }

    valid_housekeeping_genes <- housekeeping_genes[housekeeping_genes %in% rownames(counts)]
    if (length(valid_housekeeping_genes) == 0) {
        stop("None of the housekeeping genes are found in the dataset.")
    }

    # Calculate housekeeping genes' sum expression per cell
    housekeeping_sum <- Matrix::colSums(counts[valid_housekeeping_genes, , drop = FALSE])

    # Calculate target gene expression
    target_expression <- counts[target_feature, ]

    # Normalize target gene expression by housekeeping genes
    target_normalized <- target_expression / housekeeping_sum * 10000

    # Log-transform the normalized expression
    target_normalized_log <- log1p(target_normalized)

    # Return the normalized expression vector
    return(target_normalized_log)
}

# # Example usage with a Seurat object
# # Assuming 'seurat_obj' is your Seurat object
# housekeeping_genes <- c("ACTB", "GAPDH", "PGK1", "PPIA", "RPLP0", "B2M", "SDHA", "TFRC", "GUSB", "HMBS", "HPRT1", "TBP")
# normalized_expression <- normalize_feature_by_housekeeping_genes(seurat_obj, "INS", housekeeping_genes)

# # Example usage with a dgCMatrix
# # Assuming 'counts_matrix' is your dgCMatrix
# normalized_expression <- normalize_feature_by_housekeeping_genes(counts_matrix, "INS", housekeeping_genes)
