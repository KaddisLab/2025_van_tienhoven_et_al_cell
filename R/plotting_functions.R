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
            strip.text = element_text(size = 16),
            strip.text.x = element_text(size = 20),
            strip.text.y = element_text(size = 20),
            legend.key.size = unit(12, "mm"),
            legend.text = element_text(size = 20),
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(size = 26, face = "bold"),
            plot.subtitle = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
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
        ggplot(aes(x = !!sym(log2FoldChange), y = -log10(!!sym(pvalue)), label = symbol)) +
        geom_point(aes(color = significant, size = significant, alpha = significant)) +
        geom_text_repel(min.segment.length = 1, size = 4) +
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
sig_genes_plot <- function(.data,
                           .abundance = "RNA_counts_scaled",
                           .transcripts = ".feature",
                           .sample = "orig.ident",
                           FDR = "padj",
                           logFC = "logFC",
                           group_by = "diabetes_status",
                           fill_palette = "diabetes_palette",
                           title = NULL,
                           log = TRUE,
                           n_genes = 12) {
    require(dplyr)
    require(tidybulk)
    require(ggtext)


    # Get top genes
    top_genes_symbols <- .data %>%
        tidybulk::pivot_transcript() %>%
        dplyr::arrange(desc(abs(!!sym(logFC))), FDR) %>%
        dplyr::pull(!!sym(.transcripts)) %>%
        head(n_genes)

    # Prepare the data
    plot_data <- .data %>%
        dplyr::filter(!!sym(.transcripts) %in% top_genes_symbols) |>
        dplyr::mutate(!!sym(.transcripts) := factor(!!sym(.transcripts),
            levels = top_genes_symbols,
            ordered = TRUE
        ))


    # Apply log transformation if specified
    if (log) {
        plot_data <- plot_data %>%
            mutate(!!sym(.abundance) := log(!!sym(.abundance) + 1))
    }

    # Create the plot
    p <- plot_data %>%
        ggplot(aes(!!sym(group_by), !!sym(.abundance), fill = !!sym(group_by), label = !!sym(.sample))) +
        geom_boxplot(aes(alpha = 0.5), outlier.shape = NA) +
        geom_jitter(width = 0.2, shape = 21, size = 3, color = "ivory") +
        facet_wrap(as.formula(paste("~", .transcripts))) +
        ggtitle(title) +
        custom_theme +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_textbox_simple()
        ) +
        labs(y = ifelse(log, "log(Scaled counts + 1)", "Scaled counts"), x = NULL)

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
#---------------------------------------------------------
perform_gsea <- function(.data,
                         title = NULL,
                         .transcripts = ".feature",
                         .sample = "orig.ident",
                         .log2fc = "log2FoldChange",
                         log2fc_threshold = 0,
                         analysis_type = "GO",
                         ont = "BP") {
    require(dplyr)
    require(clusterProfiler)
    require(org.Hs.eg.db)
    require(ggplot2)
    require(patchwork)
    require(tidybulk)
    require(showtext)
    font_add_google(name = "Roboto", family = "Roboto")
    font_add_google(name = "Roboto Condensed", family = "Roboto Condensed")
    showtext_auto()

    # Validate ont parameter
    if (analysis_type == "GO" && !(ont %in% c("BP", "CC", "MF"))) {
        stop("Invalid ont parameter. Use 'BP' for biological processes, 'CC' for cellular components, or 'MF' for molecular functions.")
    }

    # Extract gene list for GSEA
    gsea_gene_list <- .data %>%
        dplyr::filter(!is.na(!!sym(.log2fc)), abs(!!sym(.log2fc)) >= log2fc_threshold) %>%
        dplyr::arrange(desc(!!sym(.log2fc))) %>%
        dplyr::select(!!sym(.transcripts), !!sym(.log2fc)) %>%
        dplyr::left_join(
            clusterProfiler::bitr(pull(.data, !!sym(.transcripts)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
            by = setNames("SYMBOL", .transcripts)
        ) %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(ENTREZID)) %>%
        {
            stats::setNames(pull(., !!sym(.log2fc)), .$ENTREZID)
        }

    # Perform GSEA with error handling
    gsea_result <- tryCatch(
        {
            if (analysis_type == "GO") {
                clusterProfiler::gseGO(
                    geneList = gsea_gene_list,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = ont,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    seed = 42, eps = 1e-10
                )
            } else if (analysis_type == "WP") {
                clusterProfiler::gseWP(
                    geneList = gsea_gene_list,
                    organism = "Homo sapiens",
                    pvalueCutoff = 0.05,
                    minGSSize = 10,
                    maxGSSize = 500,
                    verbose = TRUE,
                    seed = 42, eps = 1e-10
                )
            } else if (analysis_type == "KEGG") {
                clusterProfiler::gseKEGG(
                    geneList = gsea_gene_list,
                    organism = "hsa",
                    keyType = "kegg",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    seed = 42, eps = 1e-10
                )
            } else {
                stop("Invalid analysis_type. Use 'GO', 'WP', or 'KEGG'.")
            }
        },
        error = function(e) {
            warning("Error in GSEA analysis: ", e$message)
            return("error")
        }
    )

    # Define custom theme
    custom_theme <- theme_minimal() +
        theme(
            text = element_text(family = "Roboto"),
            plot.title = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.65, "lines"),
            legend.spacing.x = unit(0.3, "lines"),
            legend.spacing.y = unit(0.3, "lines"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
        )

    # Custom function to wrap text and capitalize only the first letter of the first word
    custom_label_wrap <- function(width) {
        function(labels) {
            wrapped_labels <- scales::label_wrap(width)(labels)
            sapply(wrapped_labels, function(label) {
                stringr::str_c(stringr::str_to_upper(stringr::str_sub(label, 1, 1)), stringr::str_sub(label, 2))
            })
        }
    }

    # Create separate plots for activated and suppressed gene sets
    create_plot <- function(data, direction) {
        if (is.null(data) || nrow(data) == 0) {
            ggplot() +
                geom_text(aes(x = 0.5, y = 0.5, label = paste0("No ", ifelse(direction == "Positive", "activated", "suppressed"), " gene sets")), size = 4) +
                theme_void() +
                labs(title = ifelse(direction == "Positive", "Activated gene sets", "Suppressed gene sets"))
        } else {
            clusterProfiler::dotplot(data, showCategory = 10) +
                custom_theme +
                theme(axis.text.y = element_text(size = 12, family = "Roboto Condensed", hjust = 1, lineheight = 0.5)) +
                scale_y_discrete(labels = custom_label_wrap(40)) +
                labs(title = ifelse(direction == "Positive", "Activated gene sets", "Suppressed gene sets")) +
                theme(plot.title = element_text(face = "plain", size = 14)) +
                guides(
                    fill = guide_colorbar(
                        title.position = "top",
                        title.hjust = 0.5,
                        label.theme = element_text(angle = 45, hjust = 1)
                    ),
                    size = guide_legend(
                        title.position = "top", nrow = 1,
                        title.hjust = 0.5,
                        label.position = "bottom",
                        label.spacing = unit(0.1, "lines"),
                        override.aes = list(label.size = 3)
                    )
                )
        }
    }

    # Handle GSEA results
    if (identical(gsea_result, "error")) {
        error_plot <- ggplot() +
            geom_text(aes(x = 0.5, y = 0.5, label = "Error in GSEA analysis. Check logs."), size = 4) +
            theme_void() +
            labs(title = "GSEA Error")

        final_plot <- error_plot +
            plot_annotation(
                title = title,
                subtitle = switch(analysis_type,
                    "GO" = paste0("GO: ", switch(ont,
                        "BP" = "Biological Processes",
                        "CC" = "Cellular Components",
                        "MF" = "Molecular Functions"
                    )),
                    "WP" = "WikiPathways",
                    "KEGG" = "KEGG Pathways"
                ),
                theme = theme(
                    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Roboto"),
                    plot.subtitle = element_text(size = 16, hjust = 0.5, family = "Roboto")
                )
            )

        return(list(gsea_results = NULL, plot = final_plot))
    } else {
        # Split results into activated and suppressed gene sets
        activated_sets <- gsea_result %>%
            dplyr::filter(NES > 0) %>%
            dplyr::arrange(desc(abs(NES)))
        suppressed_sets <- gsea_result %>%
            dplyr::filter(NES < 0) %>%
            dplyr::arrange(desc(abs(NES)))

        # Create plots
        activated_plot <- create_plot(activated_sets, "Positive")
        suppressed_plot <- create_plot(suppressed_sets, "Negative")

        # Combine plots using patchwork with explicit dimensions
        combined_plot <- (activated_plot | suppressed_plot) +
            plot_layout(widths = c(1, 1)) +
            plot_annotation(
                title = title,
                subtitle = switch(analysis_type,
                    "GO" = paste0("GO: ", switch(ont,
                        "BP" = "Biological Processes",
                        "CC" = "Cellular Components",
                        "MF" = "Molecular Functions"
                    )),
                    "WP" = "WikiPathways",
                    "KEGG" = "KEGG Pathways"
                ),
                theme = theme(
                    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Roboto"),
                    plot.subtitle = element_text(size = 16, hjust = 0.5, family = "Roboto")
                )
            )

        # Set explicit dimensions for the plot
        final_plot <- combined_plot & theme(plot.background = element_rect(fill = "white", color = NA))

        return(list(gsea_results = gsea_result, plot = final_plot))
    }
}
