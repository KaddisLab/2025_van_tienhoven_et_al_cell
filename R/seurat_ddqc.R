# Data driven QC (ddqc) for single cell RNA-seq data
# Adapted from https://raw.githubusercontent.com/ayshwaryas/ddqc_R/master/R/ddqc.R
# Citation https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02820-w

.ddqcBoxplot <- function(df.qc, metric.name, h.line.x = 0, do.log = FALSE) {
    require(Seurat)
    require(ggplot2)
    plt.data <- data.frame(metric = df.qc[[metric.name]], clusters = df.qc$cluster_labels)
    colnames(plt.data) <- c("metric", "clusters")

    plt.data$failed_qc <- !df.qc[[paste0(metric.name, ".passed.qc")]]

    if (do.log) {
        plt.data$metric <- log2(plt.data$metric)
        axis.labels <- labs(y = paste0("log2(", metric.name, ")"))
    } else {
        axis.labels <- labs(y = metric.name)
    }
    plt.data$clusters <- with(plt.data, reorder(clusters, -metric, mean))

    horizontal_line <- NULL
    if (h.line.x > 0) {
        horizontal_line <- geom_hline(yintercept = h.line.x, color = "red", linewidth = 0.5)
    }

    boxplot <- ggplot(plt.data, aes(x = clusters, y = metric)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color = failed_qc), alpha = 0.5, size = 0.5, width = 0.2) +
        axis.labels +
        horizontal_line +
        ggtitle(metric.name) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_color_manual(values = c("black", "red"))
}

.clusterData <- function(data, norm.factor = 10000, n.pcs = 50, k.param = 20, res = 1, random.seed = 42) {
    require(Seurat)
    set.seed(random.seed)
    data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize", scale.factor = norm.factor, verbose = FALSE)
    data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    all.genes <- rownames(x = data)
    data <- Seurat::ScaleData(data, features = all.genes, verbose = FALSE)
    data <- Seurat::RunPCA(data, npcs = n.pcs, features = VariableFeatures(data), verbose = FALSE)
    data <- Seurat::FindNeighbors(data, dims = 1:n.pcs, k.param = k.param, verbose = FALSE)
    data <- Seurat::FindClusters(data, resolution = res, verbose = FALSE)
    return(data)
}

.metricFilter <- function(
    data, df.qc, param = 2, metric.name, do.upper.co = FALSE, do.lower.co = FALSE,
    lower.bound = 10E10, upper.bound = -10E10) {
    passed.qc <- vector(mode = "logical", length = length(colnames(data)))
    names(passed.qc) <- colnames(data)

    df.qc[[metric.name]] <- data[[metric.name]][[metric.name]]
    df.qc[[paste0(metric.name, ".upper.co")]] <- NaN
    df.qc[[paste0(metric.name, ".lower.co")]] <- NaN
    df.qc[[paste0(metric.name, ".passed.qc")]] <- FALSE

    for (cl in levels(data$seurat_clusters)) {
        idx <- data$seurat_clusters == cl
        values <- data[, idx][[metric.name]][[metric.name]]

        median.v <- median(values)
        mad.v <- mad(values)
        lower.co <- min(median.v - param * mad.v, lower.bound)
        upper.co <- max(median.v + param * mad.v, upper.bound)

        qc.pass.cl <- vector(mode = "logical", length = length(values))
        qc.pass.cl <- TRUE

        if (do.lower.co) {
            qc.pass.cl <- qc.pass.cl & (values >= lower.co)
            df.qc[idx, ][[paste0(metric.name, ".lower.co")]] <- lower.co
        }
        if (do.upper.co) {
            qc.pass.cl <- qc.pass.cl & (values <= upper.co)
            df.qc[idx, ][[paste0(metric.name, ".upper.co")]] <- upper.co
        }

        df.qc[idx, ][[paste0(metric.name, ".passed.qc")]] <- qc.pass.cl
        passed.qc[idx] <- qc.pass.cl
    }

    return(df.qc)
}

#' Filter the Seurat object
#'
#' This function filters Seurat object based on ddqc and scDoublet results
#'
#' @param data Seurat object or path to a Seurat object file
#' @param ddqc_out Path to the CSV file with ddqc statistics
#' @param scDoublet_out Path to the CSV file with scDoublet results (optional)
#'
#' @return Path to the filtered Seurat object saved as a Quick Save file
#' @export
seurat_filter_qc <- function(data, ddqc_out, scDoublet_out = NULL) {
    data <- load_seurat(data)
    df.qc <- read.csv(ddqc_out, stringsAsFactors = FALSE)
    if (!is.null(scDoublet_out)) {
        scDoublet_out <- read.csv(scDoublet_out, stringsAsFactors = FALSE)
        # Ensure scDoublet_out and df.qc have the same cell order as data
        scDoublet_out <- scDoublet_out[match(colnames(data), scDoublet_out$cell), ]
        df.qc <- df.qc[match(colnames(data), df.qc$cell), ]

        keep_cells <- (scDoublet_out[["class"]] != "doublet") & df.qc$passed.qc
    } else {
        # Ensure df.qc has the same cell order as data
        df.qc <- df.qc[match(colnames(data), df.qc$cell), ]

        keep_cells <- df.qc$passed.qc
    }

    data <- subset(data, cells = colnames(data)[keep_cells])

    sample_id <- data$orig.ident[1]

    data_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_ddqc.qs")
    dir.create(dirname(data_path), showWarnings = FALSE, recursive = TRUE)
    qs::qsave(data, file = data_path)

    return(data_path)
}

#' Calculate which cells are passing ddqc
#'
#' This function takes a Seurat object after InitialQC, and then performs ddqc on it
#' Returns a path to the CSV file with ddqc statistics
#'
#' @param data Seurat object, or a path to a Seurat object file
#' @param n.pcs number of principal components for clustering. 50 by default
#' @param k.param k for FindNeighbors. 20 by default
#' @param res clustering resolution. 1 by default
#' @param threshold MAD multiplier for ddqc. 2 by default
#' @param do.plots whether to generate plots. TRUE by default
#' @param do.counts whether to consider nCount_RNA for ddqc. TRUE by default
#' @param do.genes whether to consider nFeature_RNA for ddqc. TRUE by default
#' @param do.mito whether to consider percent_mt for ddqc. TRUE by default
#' @param do.ribo whether to consider percent_rb for ddqc. TRUE by default
#' @param n.genes.lower.bound bound for lower nFeature_RNA cluster-level threshold. 200 by default
#' @param percent.mito.upper.bound bound for upper percent_mt cluster-level threshold. 15 by default
#' @param random.state random seed for clustering results reproducibility. 42 by default
#'
#' @return Path to the CSV file with ddqc statistics
#' @export
seurat_ddqc_metrics <- function(data, n.pcs = 50, k.param = 20, res = 1, threshold = 2, do.plots = TRUE, do.counts = TRUE, do.genes = TRUE, do.mito = TRUE, do.ribo = TRUE,
                         n.genes.lower.bound = 200, percent.mito.upper.bound = 15, random.state = 42) {
    data <- .clusterData(load_seurat(data), res = res, n.pcs = n.pcs, k.param = k.param, random.seed = random.state)
    sample_id <- data[[]]$orig.ident[1]
    df.qc <- data.frame("cluster_labels" = data$seurat_clusters, row.names = colnames(data))
    passed.qc <- vector(mode = "logical", length = length(data$seurat_clusters))
    passed.qc <- TRUE

    plots <- list() # Initialize an empty list to store plots

    if (do.counts) {
        df.qc <- .metricFilter(data, df.qc, threshold, "nCount_RNA", do.lower.co = TRUE)
        if (do.plots) {
            plots[["nCount_RNA"]] <- .ddqcBoxplot(df.qc, "nCount_RNA")
        }
        passed.qc <- passed.qc & df.qc$nCount_RNA.passed.qc
    }
    if (do.genes) {
        df.qc <- .metricFilter(data, df.qc, threshold, "nFeature_RNA",
            do.lower.co = TRUE,
            lower.bound = n.genes.lower.bound
        )
        if (do.plots) {
            plots[["nFeature_RNA"]] <- .ddqcBoxplot(df.qc, "nFeature_RNA", log2(n.genes.lower.bound), TRUE)
        }
        passed.qc <- passed.qc & df.qc$nFeature_RNA.passed.qc
    }
    if (do.mito) {
        df.qc <- .metricFilter(data, df.qc, threshold, "percent_mt",
            do.upper.co = TRUE,
            upper.bound = percent.mito.upper.bound
        )
        if (do.plots) {
            plots[["percent_mt"]] <- .ddqcBoxplot(df.qc, "percent_mt", percent.mito.upper.bound, FALSE)
        }
        passed.qc <- passed.qc & df.qc$percent_mt.passed.qc
    }
    if (do.ribo) {
        df.qc <- .metricFilter(data, df.qc, threshold, "percent_rb", do.upper.co = TRUE)
        if (do.plots) {
            plots[["percent_rb"]] <- .ddqcBoxplot(df.qc, "percent_rb")
        }
        passed.qc <- passed.qc & df.qc$percent_rb.passed.qc
    }
    df.qc[["passed.qc"]] <- passed.qc

    if (do.plots) {
        total_cells <- length(passed.qc)
        kept_cells <- sum(passed.qc)
        dropped_cells <- total_cells - kept_cells
        percent_dropped <- round((dropped_cells / total_cells) * 100, 1)
        combined_plot <- patchwork::wrap_plots(plots, ncol = 2) &
            patchwork::plot_annotation(
                title = paste0(sample_id, " ddqc plot"),
                subtitle = glue::glue("Filtering {threshold} Median Absolute Deviations (MADs) per cluster; {kept_cells} cells kept vs <span style='color:red;'>{dropped_cells} dropped</span> (will remove <span style='color:red;'>{percent_dropped}%</span> of {total_cells} total)"),
                theme = theme(
                    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                    plot.subtitle = ggtext::element_markdown(hjust = 0, size = 10)
                )
            )
        combined_plot_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_ddqc.png")
        dir.create(dirname(combined_plot_path), showWarnings = FALSE, recursive = TRUE)
        ggsave(filename = combined_plot_path, plot = combined_plot, width = 10, height = 10)
    }
    ddqc_table_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_ddqc.csv")
    write.csv(df.qc |> tibble::rownames_to_column(var = "cell"), ddqc_table_path, row.names = FALSE, quote = FALSE)

    return(ddqc_table_path)
}
