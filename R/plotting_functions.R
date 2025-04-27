#' Create a Volcano Plot for Differential Expression Analysis
#'
#' Creates a volcano plot showing log2 fold changes versus significance levels,
#' highlighting genes that pass both fold change and significance thresholds.
#' Significant genes are labeled on the plot.
#'
#' @param .data A data frame containing differential expression results
#' @param .transcripts Column name for gene/transcript identifiers (default: ".feature")
#' @param pvalue Column name for p-values (default: "pvalue")
#' @param FDR Column name for adjusted p-values (default: "padj")
#' @param FDR_cutoff Significance threshold for FDR (default: 0.05)
#' @param log2FoldChange Column name for log2 fold changes (default: "log2FoldChange")
#' @param log2FoldChange_cutoff Minimum absolute log2 fold change threshold (default: 1)
#' @param title Plot title (default: NULL)
#' @param subtitle Plot subtitle (default: "Upregulated in <- control \\ case ->")
#' @param xlim Vector of length 2 for x-axis limits (default: NULL)
#' @param ylim Vector of length 2 for y-axis limits (default: NULL)
#'
#' @return A ggplot object representing the volcano plot
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggtext element_textbox_simple
#' @importFrom dplyr arrange filter mutate
#' @export 
volcano_plot <- function(.data,
                         .transcripts = ".feature",
                         pvalue = "pvalue",
                         FDR = "padj",
                         FDR_cutoff = 0.05,
                         log2FoldChange = "log2FoldChange",
                         log2FoldChange_cutoff = 1,
                         title = NULL,
                         subtitle = "Upregulated in <- control \\ case ->",
                         xlim = NULL,
                         ylim = NULL) {
    require(tidyverse)
    require(ggrepel)
    require(ggtext)

    custom_theme <- theme_bw() +
        theme(
            panel.spacing.x = unit(0.5, "lines"),
            strip.text = element_text(size = rel(1.2)),
            strip.text.x = element_text(size = rel(1.2)),
            strip.text.y = element_text(size = rel(1.2)),
            legend.title = element_blank(),
            plot.title = element_text(size = rel(1.2), face = "bold"),
            axis.title = element_text(
                margin = margin(t = 0L, r = 0L, b = 0L, l = 0L)
            ),
            text = element_text(family = "Roboto")
        )
    # Get significant genes
    sig_genes_symbols <- .data %>%
        arrange(!!sym(FDR)) %>%
        dplyr::filter(!!sym(FDR) < FDR_cutoff & abs(!!sym(log2FoldChange)) >= log2FoldChange_cutoff) %>%
        pull(!!sym(.transcripts))

    # Prepare the data
    plot_data <- .data %>%
        mutate(significant = !!sym(FDR) < FDR_cutoff & abs(!!sym(log2FoldChange)) >= log2FoldChange_cutoff) %>%
        mutate(symbol = ifelse(!!sym(.transcripts) %in% sig_genes_symbols, as.character(!!sym(.transcripts)), ""))

    # Create the plot
    p <- plot_data %>%
        ggplot(aes(
            x = !!sym(log2FoldChange),
            y = -log10(!!sym(pvalue)),
            label = symbol
        )) +
        geom_point(aes(
            color = significant,
            size = significant,
            alpha = significant
        )) +
        geom_text_repel(min.segment.length = 1L, size = rel(2L)) +
        custom_theme +
        scale_color_manual(values = c("black", "#e11f28")) +
        scale_size_discrete(range = c(0.5, 2)) +
        labs(
            title = title,
            subtitle = subtitle,
            x = "Log2 Fold Change",
            y = "-Log10 p-value"
        ) +
        theme(legend.position = "none", plot.title = element_textbox_simple())

    # Add x and y limits if specified
    if (!is.null(xlim)) {
        p <- p + xlim(xlim)
    }
    if (!is.null(ylim)) {
        p <- p + ylim(ylim)
    }

    return(p)
}

# ---------
#' Plot Significant Genes from RNA-seq Data
#'
#' This function creates a boxplot of the most significant genes from RNA-seq data,
#' allowing for optional normalization and log transformation of the abundance values.
#'
#' @param .data A data frame or tibble containing the RNA-seq data.
#' @param .abundance Character string specifying the column name for abundance values. Default is "RNA_counts".
#' @param .transcripts Character string specifying the column name for transcript identifiers. Default is ".feature".
#' @param .sample Character string specifying the column name for sample identifiers. Default is "orig.ident".
#' @param FDR Character string specifying the column name for adjusted p-values. Default is "padj".
#' @param logFC Character string specifying the column name for log fold change values. Default is "logFC".
#' @param group_by Character string specifying the column name for the grouping variable. Default is "diabetes_status".
#' @param fill_palette Character string specifying the name of a color palette in the global environment. Default is "diabetes_palette".
#' @param title Character string for the plot title. Default is NULL.
#' @param n_genes Integer specifying the number of top genes to plot. Default is 12.
#' @param normalization Character string specifying the normalization method. Options are "none", "CPM", or "library_size". Default is "none".
#' @param log_transform Logical indicating whether to apply log2 transformation. Default is TRUE.
#'
#' @details
#' The function first selects the top genes based on significance and fold change.
#' It then applies normalization (if specified) and log transformation (if log_transform = TRUE) to the abundance values.
#'
#' Normalization options:
#' - "none": Uses raw counts
#' - "CPM": Counts Per Million
#' - "library_size": Divides counts by library size
#'
#' Log transformation: log2(abundance + 1)
#'
#' The function also extends the y-axis slightly to accommodate significance annotations without overlapping with facet labels.
#'
#' @return A ggplot object representing the boxplot of significant genes.
#'
#' @import dplyr
#' @import tidybulk
#' @import ggtext
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # Assuming 'rna_seq_data' is your dataset
#' sig_genes_plot(rna_seq_data,
#'     group_by = "treatment_group",
#'     fill_palette = "treatment_colors"
#' )
#'
#' # With TPM normalization
#' sig_genes_plot(rna_seq_data, normalization = "TPM")
#'
#' # Without log transformation
#' sig_genes_plot(rna_seq_data, log_transform = FALSE)
#'
#' # Plotting top 20 genes with CPM normalization
#' sig_genes_plot(rna_seq_data, n_genes = 20, normalization = "CPM")
#' }
#'
#' @export
sig_genes_plot <- function(.data,
                           .abundance = "RNA_counts",
                           .transcripts = ".feature",
                           .sample = "orig.ident",
                           FDR = "padj",
                           logFC = "logFC",
                           group_by = "diabetes_status",
                           fill_palette = "diabetes_palette",
                           title = NULL,
                           n_genes = 12,
                           normalization = c("none", "CPM", "library_size"),
                           log_transform = TRUE) {
    require(dplyr)
    require(tidybulk)
    require(ggtext)
    require(ggplot2)

    normalization <- match.arg(normalization)

    # Get top genes
    top_genes_symbols <- .data %>%
        tidybulk::tidybulk(.transcript = !!sym(.transcripts), .sample = !!sym(.sample), .abundance = !!sym(.abundance)) %>%
        tidybulk::pivot_transcript() %>%
        dplyr::mutate(significant = p_val_adj < 0.05 & abs(avg_log2FC) >= 1) %>%
        dplyr::arrange(desc(significant), desc(abs(avg_log2FC))) %>%
        dplyr::pull(!!sym(.transcripts)) %>%
        head(n_genes)

    # Prepare the data
    plot_data <- .data %>%
        dplyr::filter(!!sym(.transcripts) %in% top_genes_symbols) %>%
        dplyr::mutate(!!sym(.transcripts) := factor(!!sym(.transcripts),
            levels = top_genes_symbols,
            ordered = TRUE
        ))

    # Apply normalization if specified
    if (normalization != "none") {
        plot_data <- plot_data %>%
            group_by(!!sym(.sample)) %>%
            mutate(
                library_size = sum(!!sym(.abundance)),
                normalized_counts = case_when(
                    normalization == "CPM" ~ !!sym(.abundance) / library_size * 1e6,
                    normalization == "library_size" ~ !!sym(.abundance) / library_size
                )
            ) %>%
            ungroup()
    } else {
        plot_data <- plot_data %>%
            mutate(normalized_counts = !!sym(.abundance))
    }

    # Apply log transformation if specified
    if (log_transform) {
        plot_data <- plot_data %>%
            mutate(transformed_counts = log2(normalized_counts + 1))
    } else {
        plot_data <- plot_data %>%
            mutate(transformed_counts = normalized_counts)
    }

    # Calculate significance and prepare annotation data
    sig_data <- plot_data %>%
        group_by(!!sym(.transcripts)) %>%
        summarise(
            p_value = min(!!sym(FDR)),
            y_min = min(transformed_counts),
            y_max = max(transformed_counts),
            y_position = y_max + 0.05 * (y_max - y_min), # 5% above max value
            significance = case_when(
                p_value < 0.001 ~ "***",
                p_value < 0.01 ~ "**",
                p_value < 0.05 ~ "*",
                TRUE ~ "ns"
            )
        ) %>%
        mutate(
            xmin = 1.33,
            xmax = 1.67,
            ymin = y_position,
            ymax = y_position
        )

    # Extend y-axis to ensure annotation doesn't overlap with facet label
    y_range <- max(sig_data$y_position) - min(plot_data$transformed_counts)
    y_extension <- y_range * 0.0333 # Add about 3.33% to the top of the plot (1/3 of previous 10%)

    # Create the plot
    p <- plot_data %>%
        ggplot(aes(
            !!sym(group_by),
            transformed_counts,
            fill = !!sym(group_by),
            label = !!sym(.sample)
        )) +
        geom_boxplot(aes(alpha = 0.5), outlier.shape = NA) +
        geom_jitter(width = 0.2, shape = 21, size = 3, color = "ivory") +
        geom_text(
            data = sig_data,
            inherit.aes = FALSE,
            aes(x = 1.5, y = y_position, label = significance),
            vjust = -0.5
        ) +
        geom_segment(
            data = sig_data,
            inherit.aes = FALSE,
            aes(x = xmin, xend = xmax, y = ymin, yend = ymax),
            size = 0.3
        ) +
        facet_wrap(as.formula(paste("~", .transcripts)), scales = "free_y") +
        ggtitle(title) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_textbox_simple(
                margin = margin(b = 10),
                fill = NA,
                box.color = NA
            ),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(1, "lines"),
            axis.line = element_line(size = 0.3),
            panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            panel.background = element_blank()
        ) +
        labs(
            y = case_when(
                log_transform ~
                    paste0(
                        "log2(",
                        case_when(
                            normalization == "none" ~ "counts",
                            normalization == "TPM" ~ "TPM",
                            normalization == "CPM" ~ "CPM",
                            normalization == "library_size" ~
                                "library size normalized counts"
                        ),
                        " + 1)"
                    ),
                !log_transform ~
                    case_when(
                        normalization == "none" ~ "Raw counts",
                        normalization == "TPM" ~ "TPM",
                        normalization == "CPM" ~ "CPM",
                        normalization == "library_size" ~
                            "Library size normalized counts"
                    )
            ),
            x = NULL
        ) +
        coord_cartesian(
            clip = "off", # Allows drawing outside the plot area
            ylim = c(NA, max(sig_data$y_position) + y_extension)
        ) # Extend y-axis

    # Apply custom fill palette if provided
    if (!is.null(fill_palette)) {
        if (exists(fill_palette, envir = .GlobalEnv)) {
            palette_values <- get(fill_palette, envir = .GlobalEnv)
            p <- p + scale_fill_manual(values = palette_values)
        } else {
            warning(paste("The specified fill_palette '", fill_palette, "' does not exist in the global environment. Using default palette."))
        }
    }

    return(p)
}
