## Pinched from https://raw.githubusercontent.com/ayshwaryas/ddqc_R/master/R/ddqc.R
suppressPackageStartupMessages({
    require(Seurat)
    require(ggplot2)
})


.clusterData <- function(data, norm.factor=10000, n.pcs=50, k.param=20, res=1, random.seed=42) {
  set.seed(random.seed)
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = norm.factor, verbose=FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000,  verbose=FALSE)
  all.genes <- rownames(x = data)
  data <- ScaleData(data, features = all.genes, verbose=FALSE)
  data <- RunPCA(data, npcs=n.pcs, features = VariableFeatures(data), verbose=FALSE)
  data <- FindNeighbors(data, dims = 1:n.pcs, k.param = k.param, verbose=FALSE)
  data <- FindClusters(data, resolution = res, verbose=FALSE)
  return(data)
}


.metricFilter <- function(data, df.qc, param=2, metric.name, do.upper.co=FALSE, do.lower.co=FALSE,
                          lower.bound=10E10, upper.bound=-10E10) {
  passed.qc <- vector(mode="logical", length=length(colnames(data)))
  names(passed.qc) <- colnames(data)

  df.qc[[metric.name]] <- data[[metric.name]][[metric.name]]
  df.qc[[paste0(metric.name, ".upper.co")]] <- NaN
  df.qc[[paste0(metric.name, ".lower.co")]] <- NaN
  df.qc[[paste0(metric.name, ".passed.qc")]] <- FALSE

  for (cl in levels(data$seurat_clusters)) {
    idx <- data$seurat_clusters == cl
    values = data[ ,idx][[metric.name]][[metric.name]]

    median.v <- median(values)
    mad.v <- mad(values)
    lower.co <- min(median.v - param * mad.v, lower.bound)
    upper.co <- max(median.v + param * mad.v, upper.bound)

    qc.pass.cl <- vector(mode="logical", length=length(values))
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


#' Calculate which cells are passing ddqc
#'
#' This function takes a Seurat object after InitialQC, and then performs ddqc on it
#' Returns a data.frame that tells which cells have passed ddqc and additional information
#'
#' @param data Seurat object
#' @param n.pcs number of principal componets for clustering. 50 by default
#' @param k.param k for FindNeighbors. 20 by default
#' @param res clustering resolution. 1 by default
#' @param threshold MAD multiplier for ddqc. 2 by default
#' @param do.counts whether to consider nCount_RNA for ddqc. TRUE by default
#' @param do.genes whether to consider nFeature_RNA for ddqc. TRUE by default
#' @param do.mito whether to consider percent_mt for ddqc. TRUE by default
#' @param do.ribo whether to consider percent_rb for ddqc. TRUE by default
#' @param n.genes.lower.bound bound for lower nFeature_RNA cluster-level threshold. 200 by default
#' @param percent.mito.upper.bound bound for upper percent_mt cluster-level threshold. 30 by default
#' @param random.state random seed for clustering results reproducibility. 42 by default
#'
#' @return data.frame with ddqc statistics
#' @export
# ddqc.metrics <- function(data, n.pcs=50, k.param=20, res=1, threshold=2, do.plots = TRUE, do.counts=TRUE, do.genes=TRUE, do.mito=TRUE, do.ribo=TRUE,
#                          n.genes.lower.bound=200, percent.mito.upper.bound=30, random.state=42) {
  
#   data <- .clusterData(data, res=res, n.pcs=n.pcs, k.param=k.param, random.seed = random.state)
#   sample_id <- data[[]]$orig.ident[1]  
#   df.qc <- data.frame("cluster_labels"=data$seurat_clusters, row.names=colnames(data))
#   passed.qc <- vector(mode="logical", length=length(data$seurat_clusters))
#   passed.qc <- TRUE

#   if (do.counts) {
#     df.qc <- .metricFilter(data, df.qc, threshold, "nCount_RNA", do.lower.co=TRUE)
#     if (do.plots) {.ddqcBoxplot(df.qc, sample_id, "nCount_RNA")}
#     passed.qc <- passed.qc & df.qc$nCount_RNA.passed.qc
#   }
#   if (do.genes) {
#     df.qc <- .metricFilter(data, df.qc, threshold, "nFeature_RNA", do.lower.co=TRUE,
#                            lower.bound=n.genes.lower.bound)
#     if (do.plots) {.ddqcBoxplot(df.qc, sample_id, "nFeature_RNA", log2(n.genes.lower.bound), TRUE)}
#     passed.qc <- passed.qc & df.qc$nFeature_RNA.passed.qc
#   }
#   if (do.mito) {
#     df.qc <- .metricFilter(data, df.qc, threshold, "percent_mt", do.upper.co=TRUE,
#                            upper.bound=percent.mito.upper.bound)
#     if (do.plots) {.ddqcBoxplot(df.qc, sample_id, "percent_mt", percent.mito.upper.bound, FALSE)}

#     passed.qc <- passed.qc & df.qc$percent_mt.passed.qc
#   }
#   if (do.ribo) {
#     df.qc <- .metricFilter(data, df.qc, threshold, "percent_rb", do.upper.co=TRUE)
#     if (do.plots) {.ddqcBoxplot(df.qc, sample_id, "percent_rb")}
#     passed.qc <- passed.qc & df.qc$percent_rb.passed.qc
#   }
#   df.qc[["passed.qc"]] <- passed.qc

#   return(df.qc)
# }


#' Filter the Seurat object
#'
#' This function filters Seurat object based on df.qc
#'
#' @param data Seurat object
#' @param df.qc result of ddqc.metrics

#' @return Filtered Seurat object
#' @export
filterData <- function(data, df.qc) {
  data[["passed.qc"]] <- df.qc$passed.qc
  data <- subset(data, subset = passed.qc)
  data[["passed.qc"]] <- NULL
  return(data)
}


.ddqcBoxplot <- function(df.qc, metric.name, h.line.x=0, do.log=FALSE) {
  plt.data <- data.frame(metric=df.qc[[metric.name]], clusters=df.qc$cluster_labels)
  colnames(plt.data) <- c("metric", "clusters")
  
  plt.data$failed_qc <- !df.qc[[paste0(metric.name, ".passed.qc")]]
  
  if (do.log) {
    plt.data$metric <- log2(plt.data$metric)
    axis.labels <- labs(y=paste0("log2(", metric.name, ")"))
  } else {
    axis.labels <- labs(y=metric.name)
  }
  plt.data$clusters = with(plt.data, reorder(clusters, -metric, mean))
  
  horizontal_line <- NULL
  if (h.line.x > 0){
    horizontal_line <- geom_hline(yintercept=h.line.x, color="red", size=0.5)
  }
  
  boxplot <- ggplot(plt.data, aes(x=clusters, y=metric)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color=failed_qc), alpha=0.5, size=0.5, width=0.2) +
    axis.labels +
    horizontal_line +
    ggtitle(metric.name) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values=c("black", "red"))
}



ddqc.metrics <- function(data, n.pcs=50, k.param=20, res=1, threshold=2, do.plots = TRUE, do.counts=TRUE, do.genes=TRUE, do.mito=TRUE, do.ribo=TRUE,
                         n.genes.lower.bound=200, percent.mito.upper.bound=15, random.state=42) {
  
  data <- .clusterData(data, res=res, n.pcs=n.pcs, k.param=k.param, random.seed = random.state)
  sample_id <- data[[]]$orig.ident[1]  
  df.qc <- data.frame("cluster_labels"=data$seurat_clusters, row.names=colnames(data))
  passed.qc <- vector(mode="logical", length=length(data$seurat_clusters))
  passed.qc <- TRUE
  
  plots <- list() # Initialize an empty list to store plots
  
  if (do.counts) {
    df.qc <- .metricFilter(data, df.qc, threshold, "nCount_RNA", do.lower.co=TRUE)
    if (do.plots) {
      plots[["nCount_RNA"]] <- .ddqcBoxplot(df.qc, "nCount_RNA")
    }
    passed.qc <- passed.qc & df.qc$nCount_RNA.passed.qc
  }
  if (do.genes) {
    df.qc <- .metricFilter(data, df.qc, threshold, "nFeature_RNA", do.lower.co=TRUE,
                           lower.bound=n.genes.lower.bound)
    if (do.plots) {
      plots[["nFeature_RNA"]] <- .ddqcBoxplot(df.qc, "nFeature_RNA", log2(n.genes.lower.bound), TRUE)
    }
    passed.qc <- passed.qc & df.qc$nFeature_RNA.passed.qc
  }
  if (do.mito) {
    df.qc <- .metricFilter(data, df.qc, threshold, "percent_mt", do.upper.co=TRUE,
                           upper.bound=percent.mito.upper.bound)
    if (do.plots) {
      plots[["percent_mt"]] <- .ddqcBoxplot(df.qc, "percent_mt", percent.mito.upper.bound, FALSE)
    }
    passed.qc <- passed.qc & df.qc$percent_mt.passed.qc
  }
  if (do.ribo) {
    df.qc <- .metricFilter(data, df.qc, threshold, "percent_rb", do.upper.co=TRUE)
    if (do.plots) {
      plots[["percent_rb"]] <- .ddqcBoxplot(df.qc, "percent_rb")
    }
    passed.qc <- passed.qc & df.qc$percent_rb.passed.qc
  }
  df.qc[["passed.qc"]] <- passed.qc
  
if (do.plots) {
  # Remove legends from individual plots
  plots <- lapply(plots, function(plot) plot + theme(legend.position = "none"))
    
  # Combine plots using cowplot
  combined_plot <- cowplot::plot_grid(plotlist = plots)
  
  combined_plot_path <- glue::glue("{analysis_cache}/ddqc_out/{sample_id}_combined_boxplots.png")
  dir.create(dirname(combined_plot_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(combined_plot_path, combined_plot, width = 10, height = 10)
}
  return(df.qc)
}
