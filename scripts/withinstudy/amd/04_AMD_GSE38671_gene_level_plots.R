# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# AMD gene-level variability, MD, heatmap and individual-array DEG figures.
# Use the fitted results from script 03 without repeating modelling or testing.
# Expression displays reuse saved log2 microarray measurements, not CPM.
stopifnot(exists("amd_gene_contrast_fit"), exists("amd_gene_results"),
          requireNamespace("ggplot2", quietly = TRUE), requireNamespace("ggrepel", quietly = TRUE),
          requireNamespace("pheatmap", quietly = TRUE))

# Use consistent study colours and readable sample/group labels in every expression display.
amd_group_colours <- c(Young_WT = "#0072B2", Young_CfhKO = "#D55E00", Aged_WT = "#009E73", Aged_CfhKO = "#CC79A7")
amd_group_labels <- c(Young_WT = "Young WT", Young_CfhKO = "Young CfhKO", Aged_WT = "Aged WT", Aged_CfhKO = "Aged CfhKO")
amd_comparison_labels <- c(Young_CfhKO_vs_WT = "Young: CfhKO vs WT", Aged_CfhKO_vs_WT = "Aged: CfhKO vs WT")
amd_figure_files <- character()

# Save the correct microarray counterpart of the RNA-seq QL diagnostic.
# Blue highlights concern unusual gene variances, not DEGs; the dark-blue curve is the prior trend.
amd_sa_plot_file <- file.path(amd_within_dir, "AMD_GSE38671_gene_level_mean_variance_diagnostic.png")
grDevices::png(amd_sa_plot_file, width = 2400, height = 1800, res = 300)
limma::plotSA(amd_gene_contrast_fit, col = c("black", "#5E97F3"),
  main = "CFH-deficiency model: gene-level mean-variance diagnostic",
  xlab = "Average log2 microarray expression", ylab = "Quarter-root residual variance")
grDevices::dev.off()
amd_figure_files <- c(amd_figure_files, amd_sa_plot_file)

# Match installed limma 3.66.0 plotSA's display rule; the DEG threshold remains 0.05.
amd_sa_variance_outliers <- rep(FALSE, nrow(amd_gene_expr))
if (length(amd_gene_contrast_fit$df.prior) > 1L) {
  amd_sa_ratio <- amd_gene_contrast_fit$sigma^2 / amd_gene_contrast_fit$s2.prior
  amd_sa_tail_p <- 2 * pmin(
    stats::pf(amd_sa_ratio, df1 = amd_gene_contrast_fit$df.residual, df2 = max(amd_gene_contrast_fit$df.prior)),
    stats::pf(amd_sa_ratio, df1 = amd_gene_contrast_fit$df.residual, df2 = max(amd_gene_contrast_fit$df.prior), lower.tail = FALSE))
  amd_sa_variance_outliers <- stats::p.adjust(amd_sa_tail_p, method = "BH") <= 0.5
}

# Build both MD panels from all tested genes, selecting labels by P value within DEG direction.
amd_md_data <- do.call(rbind, lapply(names(amd_gene_results), function(comparison) {
  x <- amd_gene_results[[comparison]]
  x$Comparison <- comparison
  x$Status <- ifelse(x$adj.P.Val < 0.05, ifelse(x$logFC > 0, "Higher", "Lower"), "Not significant")
  x
}))
amd_md_data$Comparison <- factor(amd_md_data$Comparison, levels = names(amd_gene_results), labels = amd_comparison_labels)
amd_md_labels <- do.call(rbind, lapply(levels(amd_md_data$Comparison), function(comparison) {
  x <- amd_md_data[amd_md_data$Comparison == comparison & amd_md_data$Status != "Not significant", , drop = FALSE]
  rbind(head(x[x$Status == "Higher", , drop = FALSE], 5), head(x[x$Status == "Lower", , drop = FALSE], 5))
}))
amd_md_plot <- ggplot2::ggplot(amd_md_data, ggplot2::aes(AveExpr, logFC)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(data = subset(amd_md_data, Status == "Not significant"), colour = "grey65", alpha = 0.3, size = 0.55) +
  ggplot2::geom_point(data = subset(amd_md_data, Status != "Not significant"), ggplot2::aes(colour = Status), size = 2) +
  ggrepel::geom_label_repel(data = amd_md_labels, ggplot2::aes(label = gene_label, colour = Status),
    fill = "white", size = 3.3, seed = 1, box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, max.overlaps = Inf, max.time = 3, show.legend = FALSE) +
  ggplot2::facet_wrap(~ Comparison, nrow = 1) +
  ggplot2::scale_colour_manual(values = c(Higher = "#D55E00", Lower = "#0072B2"), breaks = c("Higher", "Lower")) +
  ggplot2::labs(title = "CFH-deficiency model: gene-level differential expression",
    subtitle = "Age-matched CfhKO versus WT | 4 arrays per group | 19,825 genes per comparison",
    x = "Average log2 expression (all 16 arrays)", y = "Log2 fold change (CfhKO / age-matched WT)", colour = "BH FDR < 0.05",
    caption = "Labels: up to five lowest-P DEGs per direction. Grey points do not meet FDR < 0.05.\nNo separate fold-change cutoff. More DEGs at one age is not itself an age-by-genotype interaction test.") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5), plot.caption = ggplot2::element_text(hjust = 0),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"))
amd_md_file <- file.path(amd_within_dir, "AMD_GSE38671_gene_level_MD_plots_by_age.png")
ggplot2::ggsave(amd_md_file, amd_md_plot, width = 13, height = 7, dpi = 300, bg = "white")
write.csv(amd_md_data, file.path(amd_within_dir, "AMD_GSE38671_gene_level_MD_plot_data.csv"), row.names = FALSE)
write.csv(amd_md_labels, file.path(amd_within_dir, "AMD_GSE38671_gene_level_MD_label_key.csv"), row.names = FALSE)
amd_figure_files <- c(amd_figure_files, amd_md_file)

# Show all DEGs when there are at most 30; otherwise select up to 15 per direction by P value.
# Each comparison gets its own heatmap, displaying the same 16 arrays in fixed group order.
amd_heat_samples <- amd_within_metadata[order(amd_within_metadata$group, amd_within_metadata$geo_accession), , drop = FALSE]
amd_heat_annotations <- data.frame(Group = amd_heat_samples$group, row.names = amd_heat_samples$geo_accession)
amd_heatmaps <- list()
for (comparison in names(amd_gene_degs)) {
  hits <- amd_gene_degs[[comparison]]
  higher <- hits[hits$logFC > 0, , drop = FALSE]
  lower <- hits[hits$logFC < 0, , drop = FALSE]
  if (nrow(hits) > 30L) { higher <- head(higher, 15); lower <- head(lower, 15) }
  higher$Direction <- rep("Higher", nrow(higher))
  lower$Direction <- rep("Lower", nrow(lower))
  genes <- rbind(higher, lower)
  stopifnot(nrow(genes) > 0L)
  expression <- amd_gene_expr[genes$entrez_id, amd_heat_samples$geo_accession, drop = FALSE]
  zscores <- t(scale(t(expression)))
  stopifnot(all(is.finite(zscores)))
  row_annotation <- data.frame(Direction = factor(genes$Direction, levels = c("Higher", "Lower")), row.names = genes$entrez_id)
  row_gap <- if (nrow(higher) > 0L && nrow(lower) > 0L) nrow(higher) else NULL
  heatmap <- pheatmap::pheatmap(zscores,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    breaks = seq(-max(abs(zscores)), max(abs(zscores)), length.out = 102),
    cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
    annotation_col = amd_heat_annotations, annotation_row = row_annotation,
    annotation_colors = list(Group = amd_group_colours,
      Direction = c(Higher = "#D55E00", Lower = "#0072B2")[unique(genes$Direction)]),
    annotation_names_row = FALSE, labels_row = genes$gene_label, labels_col = amd_heat_samples$geo_accession,
    gaps_row = row_gap, gaps_col = c(4, 8, 12), border_color = "white",
    fontsize = 10, fontsize_row = if (nrow(genes) > 10) 10 else 12, fontsize_col = 8, angle_col = "45",
    main = paste("CFH-deficiency model: relative DEG expression", amd_comparison_labels[[comparison]],
      paste0("All ", nrow(genes), " DEGs | Row z-scores of log2 expression | All 16 arrays"), sep = "\n"), silent = TRUE)
  prefix <- file.path(amd_within_dir, paste0("AMD_GSE38671_gene_level_", comparison, "_DEG_heatmap"))
  figure <- paste0(prefix, ".png")
  ggplot2::ggsave(figure, heatmap$gtable, width = 12, height = max(5, 3.1 + 0.23 * nrow(genes)), dpi = 300, bg = "white")
  write.csv(genes, paste0(prefix, "_genes.csv"), row.names = FALSE)
  write.csv(expression, paste0(prefix, "_log2_expression.csv"))
  write.csv(zscores, paste0(prefix, "_row_zscores.csv"))
  write.csv(amd_heat_annotations, paste0(prefix, "_sample_key.csv"))
  amd_heatmaps[[comparison]] <- list(genes = genes, expression = expression, zscores = zscores, plot = heatmap, file = figure)
  amd_figure_files <- c(amd_figure_files, figure)
}

# Select all young DEGs plus up to five higher/five lower aged DEGs by P value, showing shared genes once.
amd_dot_key <- do.call(rbind, lapply(names(amd_gene_degs), function(comparison) {
  x <- amd_gene_degs[[comparison]]
  if (nrow(x) > 10L) x <- rbind(head(x[x$logFC > 0, , drop = FALSE], 5), head(x[x$logFC < 0, , drop = FALSE], 5))
  x$Comparison <- comparison
  x
}))
amd_dot_genes <- amd_dot_key[!duplicated(amd_dot_key$entrez_id), , drop = FALSE]
amd_dot_expr <- amd_gene_expr[amd_dot_genes$entrez_id, , drop = FALSE]
amd_dot_data <- as.data.frame(as.table(amd_dot_expr), responseName = "Expression", stringsAsFactors = FALSE)
names(amd_dot_data)[1:2] <- c("entrez_id", "sample_id")
amd_dot_data$group <- amd_within_metadata$group[match(amd_dot_data$sample_id, amd_within_metadata$geo_accession)]
amd_dot_data$gene <- factor(amd_dot_data$entrez_id, levels = amd_dot_genes$entrez_id)
amd_dot_labels <- stats::setNames(amd_dot_genes$gene_label, amd_dot_genes$entrez_id)
amd_dotplot <- ggplot2::ggplot(amd_dot_data, ggplot2::aes(group, Expression, colour = group)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1), size = 2.1, alpha = 0.9) +
  ggplot2::facet_wrap(~ gene, ncol = 3, scales = "free_y", labeller = ggplot2::as_labeller(amd_dot_labels)) +
  ggplot2::scale_colour_manual(values = amd_group_colours, guide = "none") +
  ggplot2::scale_x_discrete(labels = amd_group_labels) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.16))) +
  ggplot2::labs(title = "CFH-deficiency model: individual-array expression of selected DEGs",
    subtitle = "All young DEGs + five higher/five lower aged DEGs by lowest P value\nShared genes shown once; significance comes from the corresponding age-matched test",
    x = NULL, y = "Log2 microarray expression",
    caption = "One point per array (four per group). Horizontal offsets keep points visible; expression values are unchanged.\nEach gene has its own y-axis range. Colours identify age/genotype groups; no sex shapes are inferred.") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(face = "italic", size = 11),
    plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0))
amd_dot_file <- file.path(amd_within_dir, "AMD_GSE38671_gene_level_DEG_expression_dotplots_selected.png")
ggplot2::ggsave(amd_dot_file, amd_dotplot, width = 12, height = 2.7 + 2.5 * ceiling(nrow(amd_dot_genes) / 3), dpi = 300, bg = "white")
write.csv(amd_dot_key, file.path(amd_within_dir, "AMD_GSE38671_gene_level_DEG_dotplot_selection_by_age.csv"), row.names = FALSE)
write.csv(amd_dot_genes, file.path(amd_within_dir, "AMD_GSE38671_gene_level_DEG_dotplot_gene_key.csv"), row.names = FALSE)
write.csv(amd_dot_data, file.path(amd_within_dir, "AMD_GSE38671_gene_level_DEG_dotplot_data.csv"), row.names = FALSE)
amd_figure_files <- c(amd_figure_files, amd_dot_file)
stopifnot(all(file.exists(amd_figure_files)))
print(c(Variance_display_outliers = sum(amd_sa_variance_outliers), Selected_dotplot_genes = nrow(amd_dot_genes)))
print(amd_figure_files)
