# Data driven QC (ddqc) for single cell RNA-seq data
# Adapted from https://raw.githubusercontent.com/ayshwaryas/ddqc_R/master/R/ddqc.R
# Citation https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02820-w


.ddqcBoxplot <- function(df.qc, metric.name, h.line.x = 0, do.log = FALSE, scDblFinder_out) {
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  
  plt.data <- data.frame(metric = df.qc[[metric.name]], clusters = df.qc$cluster_labels, row.names = df.qc$cell)
  
  plt.data$failed_qc <- !df.qc[[paste0(metric.name, ".passed.qc")]]
  
  if (!is.null(scDblFinder_out)) {
    plt.data$is_doublet <- scDblFinder_out[match(rownames(plt.data), scDblFinder_out$cell), "class"] == "doublet"
  } else {
    plt.data$is_doublet <- FALSE
  }
  
  plt.data$point_color <- ifelse(plt.data$failed_qc & !plt.data$is_doublet, "red",
                                ifelse(plt.data$failed_qc & plt.data$is_doublet, "darkred",
                                        ifelse(!plt.data$failed_qc & plt.data$is_doublet, "blue", "black")))
  
  plt.data$point_shape <- ifelse(plt.data$is_doublet, 17, 16)
  
  if (do.log) {
    plt.data$metric <- log1p(plt.data$metric)
    axis.labels <- labs(y = paste0("log(", metric.name, ")"))
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
    geom_point(aes(color = point_color, shape = point_shape), alpha = 0.5, size = 1.5, position = position_jitter(width = 0.2)) +
    axis.labels +
    horizontal_line +
    ggtitle(metric.name) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_identity() +
    scale_shape_identity() +
    theme_classic()
  
  return(boxplot)
}


.clusterData <- function(seurat_object, norm.factor = 10000, n.pcs = 50, k.param = 20, res = 1) {
    require(Seurat)
    set.seed(42)
    seurat_object <- Seurat::NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = norm.factor, verbose = FALSE)
    seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    all.genes <- rownames(x = seurat_object)
    seurat_object <- Seurat::ScaleData(seurat_object, features = all.genes, verbose = FALSE)
    seurat_object <- Seurat::RunPCA(seurat_object, npcs = n.pcs, features = VariableFeatures(seurat_object), verbose = FALSE)
    seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:n.pcs, k.param = k.param, verbose = FALSE)
    seurat_object <- Seurat::FindClusters(seurat_object, resolution = res, verbose = FALSE)
    seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:n.pcs)
    return(seurat_object)
}

.metricFilter <- function(seurat_object, df.qc, param = 2, metric.name, do.upper.co = FALSE, do.lower.co = FALSE,
                          lower.bound = -10E10, upper.bound = 10E10, use_mad_filtering = TRUE) {
    passed.qc <- rep(TRUE, length(colnames(seurat_object)))
    names(passed.qc) <- colnames(seurat_object)

    df.qc[[metric.name]] <- seurat_object[[metric.name]][[metric.name]]
    df.qc[[paste0(metric.name, ".upper.co")]] <- NA_real_
    df.qc[[paste0(metric.name, ".lower.co")]] <- NA_real_
    df.qc[[paste0(metric.name, ".passed.qc")]] <- FALSE

    for (cl in levels(seurat_object$seurat_clusters)) {
        idx <- seurat_object$seurat_clusters == cl
        values <- seurat_object[, idx][[metric.name]][[metric.name]]

        median.v <- median(values, na.rm = TRUE)
        
        if (use_mad_filtering) {
            mad.v <- mad(values, na.rm = TRUE)
            lower.co <- max(median.v - param * mad.v, lower.bound)
            upper.co <- min(median.v + param * mad.v, upper.bound)
            } else {
            lower.co <- lower.bound
            upper.co <- upper.bound
            }
        qc.pass.cl <- rep(TRUE, length(values))

        if (do.lower.co) {
            qc.pass.cl <- qc.pass.cl & (values >= lower.co)
            # Handle cases where values might be NA
            qc.pass.cl[is.na(values)] <- FALSE
            df.qc[idx, paste0(metric.name, ".lower.co")] <- lower.co
        }
        if (do.upper.co) {
            qc.pass.cl <- qc.pass.cl & (values <= upper.co)
            # Handle cases where values might be NA
            qc.pass.cl[is.na(values)] <- FALSE
            df.qc[idx, paste0(metric.name, ".upper.co")] <- upper.co
        }

        # Ensure qc.pass.cl does not have NA values
        qc.pass.cl[is.na(qc.pass.cl)] <- FALSE

        df.qc[idx, paste0(metric.name, ".passed.qc")] <- qc.pass.cl
        passed.qc[idx] <- qc.pass.cl
    }

    return(df.qc)
}


#' ddQC and scDblFinder filtering for Seurat objects
#'
#' This function takes a Seurat object after running seurat_cell_metrics(), runs ddQC and
#' returns a path to the CSV file with ddQC statistics, or a path to the filtered Seurat object if "do.filter" is TRUE
#'
#' @param seurat_object Seurat object, or a path to a Seurat object file
#' @param n.pcs number of principal components for clustering. 50 by default
#' @param k.param k for FindNeighbors. 20 by default
#' @param res clustering resolution. 1 by default
#' @param threshold MAD multiplier for ddqc. 3.5 by default
#' @param do.plots whether to generate plots. TRUE by default
#' @param do.counts whether to consider nCount_RNA for ddqc. TRUE by default
#' @param do.genes whether to consider nFeature_RNA for ddqc. TRUE by default
#' @param do.mito whether to consider percent_mt for ddqc. TRUE by default
#' @param do.ribo whether to consider percent_rb for ddqc. TRUE by default
#' @param n.reads.lower.bound bound for lower nCount_RNA cluster-level threshold. 500 by default
#' @param n.genes.lower.bound bound for lower nFeature_RNA cluster-level threshold. 200 by default
#' @param percent.mito.upper.bound bound for upper percent_mt cluster-level threshold. 20 by default
#' @param percent.rb.lower.bound bound for lower percent_rb cluster-level threshold. 0 by default
#' @param do.filter Logical. Filter MAD outliers and doublets if scDblFinder_out is provided. TRUE by default
#' @param scDblFinder_out Path to the CSV file with scDoublet results.
#' @param use_mad_counts Logical. Whether to use MAD filtering for nCount_RNA. TRUE by default
#' @param use_mad_genes Logical. Whether to use MAD filtering for nFeature_RNA. TRUE by default
#' @param use_mad_mito Logical. Whether to use MAD filtering for percent_mt. TRUE by default
#' @param use_mad_ribo Logical. Whether to use MAD filtering for percent_rb. TRUE by default
#'
#' @return Path to the filtered seurat_object at "{analysis_cache}/ddqc_out/{sample_id}_ddqc.qs", or the similarly named CSV file with ddqc statistics if do.filter is FALSE
#' @export
seurat_ddqc <- function(seurat_object, scDblFinder_out, n.pcs = 50, k.param = 20, res = 1, threshold = 3.5, 
                        do.plots = TRUE, do.counts = TRUE, do.genes = TRUE, do.mito = TRUE, do.ribo = TRUE,
                        n.reads.lower.bound = 500, n.genes.lower.bound = 200, percent.mito.upper.bound = 20, 
                        percent.rb.lower.bound = 0, do.filter = TRUE, 
                        use_mad_counts = TRUE, use_mad_genes = TRUE, use_mad_mito = TRUE, use_mad_ribo = TRUE) {
    message("Loading Seurat object...")
    seurat_object <- load_seurat(seurat_object)
    sample_id <- Seurat::Project(seurat_object)
    if (!is.data.frame(scDblFinder_out)) {scDblFinder_out <- read.csv(scDblFinder_out, stringsAsFactors = FALSE)}
    scDblFinder_out <- scDblFinder_out[match(colnames(seurat_object), scDblFinder_out$cell), ]
    message("Loaded scDblFinder output for ", sample_id)
    message("Running ddQC on ", sample_id,"...")
    seurat_object <- .clusterData(seurat_object, res = res, n.pcs = n.pcs, k.param = k.param)
    message("Clustered the Seurat object")
    df.qc <- data.frame("cluster_labels" = seurat_object$seurat_clusters, row.names = colnames(seurat_object)) |> mutate(cell = colnames(seurat_object), .before = 1)
    passed.qc <- rep(TRUE, nrow(seurat_object[[]]))

    plots <- list() # Initialize an empty list to store plots

    if (do.counts) {
        message("Filtering cells based on nCount_RNA")
        df.qc <- .metricFilter(seurat_object, df.qc, threshold, "nCount_RNA",
            do.lower.co = TRUE,
            lower.bound = n.reads.lower.bound,
            use_mad_filtering = use_mad_counts
        )
        if (do.plots) {
            plots[["nCount_RNA"]] <- .ddqcBoxplot(df.qc, "nCount_RNA", n.reads.lower.bound, scDblFinder_out = scDblFinder_out)
        }
        passed.qc <- passed.qc & df.qc$nCount_RNA.passed.qc
    }

    if (do.genes) {
        message("Filtering cells based on nFeature_RNA")
        df.qc <- .metricFilter(seurat_object, df.qc, threshold, "nFeature_RNA",
            do.lower.co = TRUE,
            lower.bound = n.genes.lower.bound,
            use_mad_filtering = use_mad_genes
        )
        if (do.plots) {
            plots[["nFeature_RNA"]] <- .ddqcBoxplot(df.qc, "nFeature_RNA", log1p(n.genes.lower.bound), TRUE, scDblFinder_out = scDblFinder_out)
        }
        passed.qc <- passed.qc & df.qc$nFeature_RNA.passed.qc
    }

    if (do.mito) {
        message("Filtering cells based on percent_mt")
        df.qc <- .metricFilter(seurat_object, df.qc, threshold, "percent_mt",
            do.upper.co = TRUE,
            upper.bound = percent.mito.upper.bound,
            use_mad_filtering = use_mad_mito
        )
        if (do.plots) {
            plots[["percent_mt"]] <- .ddqcBoxplot(df.qc, "percent_mt", percent.mito.upper.bound, FALSE, scDblFinder_out = scDblFinder_out)
        }
        passed.qc <- passed.qc & df.qc$percent_mt.passed.qc
    }

    if (do.ribo) {
        message("Filtering cells based on percent_rb")
        df.qc <- .metricFilter(seurat_object, df.qc, threshold, "percent_rb",
            do.lower.co = TRUE,
            lower.bound = percent.rb.lower.bound,
            use_mad_filtering = use_mad_ribo
        )
        if (do.plots) {
            plots[["percent_rb"]] <- .ddqcBoxplot(df.qc, "percent_rb", percent.rb.lower.bound, scDblFinder_out = scDblFinder_out)
        }
        passed.qc <- passed.qc & df.qc$percent_rb.passed.qc
    }
    df.qc[["passed.qc"]] <- passed.qc

    if (do.plots) {
        message("Generating combined plot")

        # Data merging and basic calculations
        df <- dplyr::inner_join(df.qc, scDblFinder_out, by = "cell")
        total_cells <- nrow(df)
        ddqc_kept_cells <- sum(df$passed.qc)
        ddqc_dropped_cells <- total_cells - ddqc_kept_cells
        ddqc_percent_dropped <- round((ddqc_dropped_cells / total_cells) * 100, 1)
        n_doublets <- sum(df$class == "doublet")
        percent_doublets <- round((n_doublets / total_cells) * 100, 1)
        total_remaining <- total_cells - ddqc_dropped_cells - n_doublets
        percent_remaining <- round((total_remaining / total_cells) * 100, 1)

        # Identify and categorize cells
        qc_fail <- df %>%
            dplyr::filter(!passed.qc) %>%
            dplyr::pull(cell)
        doublets <- df %>%
            dplyr::filter(class == "doublet") %>%
            dplyr::pull(cell)
        both <- intersect(qc_fail, doublets)
        qc_fail <- setdiff(qc_fail, both)
        doublets <- setdiff(doublets, both)

        # Assign status based on conditions
        df$status <- with(df, ifelse(cell %in% qc_fail, "QC Fail",
            ifelse(cell %in% doublets, "Doublets",
                ifelse(cell %in% both, "Both", "Passed QC")
            )
        ))

        cells_to_highlight <- list("QC Fail" = qc_fail, "Doublets" = doublets, "Both" = both)

        # Cell highlight plot configuration
        # Assuming the `plots[["umap"]]` assignment starts here
        plots[["umap"]] <- scCustomize::Cell_Highlight_Plot(
            seurat_object = seurat_object,
            label = TRUE,
            label.size = 6,
            repel = TRUE,
            cells_highlight = cells_to_highlight,
            highlight_color = c("QC Fail" = "red", "Doublets" = "blue", "Both" = "purple"),
            background_color = "gray"
        ) +
            ggplot2::labs(title = "UMAP", subtitle = "Gene expression space") +
            ggplot2::theme_classic() +
            NoLegend() +
            ggplot2::theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "plain"))


        # Scatter plot configuration
        plots[["reads_vs_genes"]] <- ggplot2::ggplot(df, aes(x = nFeature_RNA, y = nCount_RNA, color = status)) +
            ggplot2::geom_point(alpha = ifelse(df$passed.qc & df$class != "doublet", 0.2, 1)) +
            ggplot2::labs(title = "nCount_RNA vs nFeature_RNA", x = "log(nFeature_RNA)", y = "log(nCount_RNA)", color = "Quality") +
            ggplot2::scale_color_manual(values = c("Passed QC" = "black", "QC Fail" = "red", "Doublets" = "blue", "Both" = "purple")) +
            ggplot2::scale_x_log10() +
            ggplot2::scale_y_log10() +
            ggplot2::theme_classic() +
            NoLegend()

        # Combining plots
        combined_plot <- patchwork::wrap_plots(plots, ncol = 2) +
            patchwork::plot_annotation(
                title = glue::glue("{sample_id} ddQC plot"),
                subtitle = glue::glue("Filtering {total_cells} cells; <span style='color:red;'>{ddqc_dropped_cells} ({ddqc_percent_dropped}%)</span> cells over {threshold} MADs and <span style='color:blue;'>{n_doublets} ({percent_doublets}%)</span> doublets removed. {total_remaining} ({percent_remaining}%) cells remaining."),
                theme = ggplot2::theme(
                    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                    plot.subtitle = ggtext::element_markdown(hjust = 0, size = 10)
                )
            )

        # Save plot to file
        combined_plot_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_ddqc.png")
        dir.create(dirname(combined_plot_path), showWarnings = FALSE, recursive = TRUE)
        ggsave(filename = combined_plot_path, plot = combined_plot, width = 10, height = 12)
    }
    ddqc_table_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_ddqc.csv")
    write.csv(df.qc, ddqc_table_path, row.names = FALSE, quote = FALSE)
    if(do.filter) {
        message("Filtering Seurat object")
        if (!is.null(scDblFinder_out)) {
            keep_cells <- (scDblFinder_out[["class"]] != "doublet") & df.qc$passed.qc
        } else {
            keep_cells <- df.qc$passed.qc
        }
        seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[keep_cells])

        seurat_object_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_ddqc.qs")
        dir.create(dirname(seurat_object_path), showWarnings = FALSE, recursive = TRUE)
        qs::qsave(seurat_object, file = seurat_object_path)
        message("Filtered Seurat object saved to", seurat_object_path)
        return(seurat_object_path)
    }
    message("ddQC statistics saved to", ddqc_table_path)
    return(ddqc_table_path)
}
