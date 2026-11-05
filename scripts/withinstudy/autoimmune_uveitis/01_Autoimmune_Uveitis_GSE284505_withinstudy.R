# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "autoimmune_uveitis"), recursive = TRUE, showWarnings = FALSE)

# Autoimmune uveitis GSE284505: within-study gene-level analysis, 5 November 2026.
# Continue from the reviewed filtered counts and saved TMM factors.
# Method: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Locate the output folder and completed precheck without using the live R session.
uveitis_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "autoimmune_uveitis")
dir.create(uveitis_within_dir, recursive = TRUE, showWarnings = FALSE)
uveitis_precheck_file <- Sys.getenv("UVEITIS_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "autoimmune_uveitis", "Autoimmune_Uveitis_GSE284505_precheck_checkpoint_2026-11-02.RData"))
stopifnot(file.exists(uveitis_precheck_file))
uveitis_precheck <- new.env()
uveitis_loaded <- load(uveitis_precheck_file, envir = uveitis_precheck)
uveitis_within_dge <- uveitis_precheck$uveitis_dge_qc
uveitis_within_metadata <- uveitis_precheck$uveitis_metadata_selected

# Check the reviewed dimensions, sample order, gene order and saved scaling factors.
stopifnot(
  identical(dim(uveitis_within_dge), c(18187L, 6L)),
  identical(colnames(uveitis_within_dge), uveitis_within_metadata$sample_id),
  identical(rownames(uveitis_within_dge), uveitis_within_dge$genes$AccID),
  anyDuplicated(rownames(uveitis_within_dge)) == 0L,
  all(is.finite(uveitis_within_dge$samples$norm.factors)),
  all(uveitis_within_dge$samples$norm.factors > 0)
)
print(dim(uveitis_within_dge))

# Set control as reference; matching replicate numbers do not establish pairing.
# Sex is unreported in the reviewed metadata, so no sex category is inferred.
uveitis_within_metadata$group <- factor(
  uveitis_within_metadata$genotype, levels = c("CON", "EAU")
)
stopifnot(!anyNA(uveitis_within_metadata$group))
print(table(uveitis_within_metadata$group))
uveitis_design <- model.matrix(~ group, data = uveitis_within_metadata)
rownames(uveitis_design) <- uveitis_within_metadata$sample_id
stopifnot(qr(uveitis_design)$rank == ncol(uveitis_design))
print(uveitis_design)
print(nrow(uveitis_design) - qr(uveitis_design)$rank)
write.csv(uveitis_within_metadata, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_sample_key.csv"), row.names = FALSE)
write.csv(uveitis_design, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_design.csv"))

# Fit the same robust current edgeR QL pipeline used for the prior count datasets.
# Raw filtered counts and their saved TMM offsets enter the model, not log-CPM.
stopifnot(as.character(utils::packageVersion("edgeR")) == "4.8.2")
uveitis_fit <- edgeR::glmQLFit(uveitis_within_dge, design = uveitis_design, robust = TRUE)
stopifnot(identical(dim(uveitis_fit$coefficients), c(18187L, 2L)))

# Save the QL dispersion diagnostic: raw estimates, moderated estimates and trend.
uveitis_ql_plot_file <- file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_QL_dispersion.png")
png(uveitis_ql_plot_file, width = 2400, height = 1800, res = 300)
edgeR::plotQLDisp(uveitis_fit, main = "Autoimmune uveitis: QL dispersion")
dev.off()

# Test disease/injury versus control and retain all tested genes with BH FDR.
uveitis_test <- edgeR::glmQLFTest(uveitis_fit, coef = "groupEAU")
uveitis_results <- edgeR::topTags(uveitis_test, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
uveitis_results$gene_id <- rownames(uveitis_results)
stopifnot(
  identical(uveitis_results$gene_id, as.character(uveitis_results$AccID)),
  all(is.finite(uveitis_results$logFC)), all(is.finite(uveitis_results$F)),
  all(is.finite(uveitis_results$PValue)), all(is.finite(uveitis_results$FDR)),
  all(uveitis_results$FDR >= 0 & uveitis_results$FDR <= 1)
)

# Use an available name for display; keep original source IDs as unique row keys.
# Fallback labels are source identifiers, not newly verified gene-symbol mappings.
uveitis_results$display_name <- as.character(uveitis_results$gene_name)
uveitis_missing_name <- is.na(uveitis_results$display_name) | !nzchar(uveitis_results$display_name)
uveitis_results$display_name[uveitis_missing_name] <- uveitis_results$gene_id[uveitis_missing_name]
uveitis_duplicate_name <- duplicated(uveitis_results$display_name) | duplicated(uveitis_results$display_name, fromLast = TRUE)
uveitis_results$display_name[uveitis_duplicate_name] <- paste0(
  uveitis_results$display_name[uveitis_duplicate_name], " [", uveitis_results$gene_id[uveitis_duplicate_name], "]"
)

# Apply the established FDR < 0.05 criterion without an extra fold-change cutoff.
uveitis_degs <- subset(uveitis_results, FDR < 0.05)
uveitis_summary <- data.frame(
  Comparison = "EAU_vs_CON", Tested = nrow(uveitis_results),
  DEGs = nrow(uveitis_degs), Higher = sum(uveitis_degs$logFC > 0),
  Lower = sum(uveitis_degs$logFC < 0)
)
stopifnot(uveitis_summary$Higher + uveitis_summary$Lower == uveitis_summary$DEGs)
print(uveitis_summary)
print(head(uveitis_results[, c("gene_id", "display_name", "logFC", "F", "PValue", "FDR")], 10))
write.csv(uveitis_results, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_EAU_vs_CON_all_genes.csv"), row.names = FALSE)
write.csv(uveitis_degs, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_EAU_vs_CON_DEGs_FDR005.csv"), row.names = FALSE)
write.csv(uveitis_summary, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_DEG_summary.csv"), row.names = FALSE)

# Rank every tested gene using the established signed square-root F score.
# The tolerance handles only negligible numerical negatives and leaves original F intact.
uveitis_rank_F_check <- c(Minimum_F = min(uveitis_results$F), Negative_F = sum(uveitis_results$F < 0), Nonfinite_F = sum(!is.finite(uveitis_results$F)))
uveitis_rank_tolerance <- sqrt(.Machine$double.eps)
stopifnot(all(uveitis_results$F >= -uveitis_rank_tolerance))
uveitis_ranked <- transform(uveitis_results, rank_score = sign(logFC) * sqrt(pmax(F, 0)))
stopifnot(all(is.finite(uveitis_ranked$rank_score)))
uveitis_ranked <- uveitis_ranked[order(uveitis_ranked$rank_score, decreasing = TRUE), , drop = FALSE]
print(uveitis_rank_F_check)
write.csv(uveitis_ranked, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_EAU_vs_CON_ranked_all_genes.csv"), row.names = FALSE)

# Plot all tested genes; colour significant increases/decreases and label up to five each.
uveitis_md_data <- uveitis_results
uveitis_md_data$Status <- "Not significant"
uveitis_md_data$Status[uveitis_md_data$FDR < 0.05 & uveitis_md_data$logFC > 0] <- "Higher"
uveitis_md_data$Status[uveitis_md_data$FDR < 0.05 & uveitis_md_data$logFC < 0] <- "Lower"
uveitis_md_hits <- subset(uveitis_md_data, Status != "Not significant")
uveitis_md_labels <- rbind(head(subset(uveitis_md_hits, Status == "Higher"), 5), head(subset(uveitis_md_hits, Status == "Lower"), 5))
uveitis_md_plot <- ggplot2::ggplot(uveitis_md_data, ggplot2::aes(logCPM, logFC)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(data = subset(uveitis_md_data, Status == "Not significant"), colour = "grey65", alpha = 0.35, size = 0.6) +
  ggplot2::geom_point(data = uveitis_md_hits, ggplot2::aes(colour = Status), size = 1, alpha = 0.7) +
  ggrepel::geom_label_repel(data = uveitis_md_labels, ggplot2::aes(label = display_name, colour = Status),
    size = 3, seed = 1, box.padding = 0.45, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 5, max.iter = 20000,
    nudge_y = ifelse(uveitis_md_labels$logFC > 0, 1.2, -0.4), force_pull = 0.05, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(Higher = "#D55E00", Lower = "#0072B2"), breaks = c("Higher", "Lower")) +
  ggplot2::labs(title = "Autoimmune uveitis: genome-wide differential expression", subtitle = "EAU versus CON | all retained libraries",
    x = "Average log2 CPM (all 6 libraries)", y = "Log2 fold change (EAU versus CON)", colour = "FDR < 0.05",
    caption = "All tested genes shown. Labels: up to five lowest-P-value DEGs per direction.\nGrey: FDR criterion not met. No separate fold-change cutoff.") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5), plot.caption = ggplot2::element_text(hjust = 0))
uveitis_md_file <- file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_MD_plot.png")
ggplot2::ggsave(uveitis_md_file, uveitis_md_plot, width = 9, height = 6.5, units = "in", dpi = 300, bg = "white")

# Select up to 15 significant genes per direction by P value for the descriptive heatmap.
# If a direction has fewer significant genes, use all available; do not relabel non-DEGs.
uveitis_heat_genes <- rbind(
  transform(head(subset(uveitis_degs, logFC > 0), 15), Direction = rep("Higher", min(15L, sum(uveitis_degs$logFC > 0)))),
  transform(head(subset(uveitis_degs, logFC < 0), 15), Direction = rep("Lower", min(15L, sum(uveitis_degs$logFC < 0))))
)
uveitis_heat_file <- file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_DEG_expression_heatmap_top30.png")
uveitis_heat_samples <- uveitis_within_metadata[order(uveitis_within_metadata$group, uveitis_within_metadata$sample_id), , drop = FALSE]
if (nrow(uveitis_heat_genes) > 0L) {
  # Calculate normalized log-CPM from the full filtered object before selecting rows.
  uveitis_heat_log_cpm <- edgeR::cpm(uveitis_within_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)[uveitis_heat_genes$gene_id, uveitis_heat_samples$sample_id, drop = FALSE]
  uveitis_heat_z <- t(scale(t(uveitis_heat_log_cpm)))
  stopifnot(all(is.finite(uveitis_heat_z)))
  uveitis_heat_limit <- max(abs(uveitis_heat_z))
  uveitis_heat_annotations <- data.frame(Group = uveitis_heat_samples$group, row.names = uveitis_heat_samples$sample_id)
  uveitis_heat_row_annotations <- data.frame(Direction = factor(uveitis_heat_genes$Direction, levels = c("Higher", "Lower")), row.names = uveitis_heat_genes$gene_id)
  uveitis_heat_annotation_colours <- list(Group = c(CON = "#0072B2", EAU = "#D55E00"), Direction = c(Higher = "#D55E00", Lower = "#0072B2"))
  uveitis_heat_nhigher <- sum(uveitis_heat_genes$Direction == "Higher")
  uveitis_heatmap <- pheatmap::pheatmap(uveitis_heat_z,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    breaks = seq(-uveitis_heat_limit, uveitis_heat_limit, length.out = 102), scale = "none",
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = uveitis_heat_annotations, annotation_row = uveitis_heat_row_annotations,
    annotation_colors = uveitis_heat_annotation_colours, annotation_names_row = FALSE,
    labels_row = uveitis_heat_genes$display_name, labels_col = uveitis_heat_samples$sample_id,
    gaps_row = if (uveitis_heat_nhigher > 0L && uveitis_heat_nhigher < nrow(uveitis_heat_genes)) uveitis_heat_nhigher else NULL,
    gaps_col = sum(uveitis_heat_samples$group == "CON"),
    border_color = "white", fontsize = 10, fontsize_row = 9, fontsize_col = 9, angle_col = "45",
    main = paste("Autoimmune uveitis DEGs: relative expression", paste0(uveitis_heat_nhigher, " higher + ", sum(uveitis_heat_genes$Direction == "Lower"), " lower by lowest P value"), "Row z-scores of TMM log2 CPM", sep = "\n"), silent = TRUE)
  ggplot2::ggsave(uveitis_heat_file, uveitis_heatmap$gtable, width = 11, height = 10, units = "in", dpi = 300, bg = "white")
  write.csv(uveitis_heat_log_cpm, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_top30_DEG_heatmap_logCPM.csv"))
  write.csv(uveitis_heat_z, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_top30_DEG_heatmap_row_zscores.csv"))
  write.csv(uveitis_heat_annotations, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_top30_DEG_heatmap_sample_key.csv"))
} else {
  # Preserve a truthful figure if no gene meets the established DEG criterion.
  png(uveitis_heat_file, width = 2400, height = 1800, res = 300)
  plot.new(); text(0.5, 0.5, "No DEGs at FDR < 0.05; no DEG heatmap can be drawn.")
  dev.off()
}
write.csv(uveitis_heat_genes, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_top30_DEG_heatmap_selected_genes.csv"), row.names = FALSE)

# Select up to five DEGs per direction, using the same P-value ordering as MD labels.
uveitis_dot_genes <- rbind(
  transform(head(subset(uveitis_degs, logFC > 0), 5), Direction = rep("Higher", min(5L, sum(uveitis_degs$logFC > 0)))),
  transform(head(subset(uveitis_degs, logFC < 0), 5), Direction = rep("Lower", min(5L, sum(uveitis_degs$logFC < 0))))
)
uveitis_dot_file <- file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_DEG_expression_dotplots_top10.png")
if (nrow(uveitis_dot_genes) > 0L) {
  # Display observed non-log CPM, retaining one point per measured library.
  uveitis_dot_cpm <- edgeR::cpm(uveitis_within_dge, normalized.lib.sizes = TRUE, log = FALSE)[uveitis_dot_genes$gene_id, , drop = FALSE]
  uveitis_dot_data <- as.data.frame(as.table(uveitis_dot_cpm), responseName = "CPM", stringsAsFactors = FALSE)
  names(uveitis_dot_data)[1:2] <- c("gene_id", "sample_id")
  uveitis_dot_sample_index <- match(uveitis_dot_data$sample_id, uveitis_within_metadata$sample_id)
  uveitis_dot_gene_index <- match(uveitis_dot_data$gene_id, uveitis_dot_genes$gene_id)
  stopifnot(!anyNA(uveitis_dot_sample_index), !anyNA(uveitis_dot_gene_index))
  uveitis_dot_data$group <- uveitis_within_metadata$group[uveitis_dot_sample_index]
  uveitis_dot_data$display_name <- uveitis_dot_genes$display_name[uveitis_dot_gene_index]
  uveitis_dot_data$Direction <- uveitis_dot_genes$Direction[uveitis_dot_gene_index]
  uveitis_dot_data$gene <- factor(uveitis_dot_data$gene_id, levels = uveitis_dot_genes$gene_id)
  uveitis_dot_labels <- stats::setNames(uveitis_dot_genes$display_name, uveitis_dot_genes$gene_id)
  uveitis_dotplot <- ggplot2::ggplot(uveitis_dot_data, ggplot2::aes(group, CPM, colour = group)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1), size = 2.8) +
    ggplot2::facet_wrap(~ gene, ncol = 2, dir = "v", scales = "free_y", labeller = ggplot2::as_labeller(uveitis_dot_labels)) +
    ggplot2::scale_colour_manual(values = c(CON = "#0072B2", EAU = "#D55E00"), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
    ggplot2::labs(title = "Autoimmune uveitis: individual-library expression of selected DEGs",
      subtitle = "Up to five higher DEGs on left; up to five lower on right\nLowest P values within direction among genes with FDR < 0.05",
      x = NULL, y = "TMM-normalized CPM",
      caption = "One point per library; horizontal offsets keep points visible. Each gene has its own y-axis range.\nPoints show observed expression. Sex was not recorded in the reviewed metadata.") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
      strip.text = ggplot2::element_text(face = "italic", size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5), plot.caption = ggplot2::element_text(hjust = 0))
  ggplot2::ggsave(uveitis_dot_file, uveitis_dotplot, width = 11, height = 12, units = "in", dpi = 300, bg = "white")
  write.csv(uveitis_dot_data, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_DEG_dotplot_top10_data.csv"), row.names = FALSE)
} else {
  # Do not present nonsignificant genes as significant in an empty-result analysis.
  png(uveitis_dot_file, width = 2400, height = 1800, res = 300)
  plot.new(); text(0.5, 0.5, "No DEGs at FDR < 0.05; no DEG dot plots can be drawn.")
  dev.off()
}
write.csv(uveitis_dot_genes, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_top10_DEG_dotplot_selected_genes.csv"), row.names = FALSE)

# Save only this study's working objects, excluding the precheck environment.
uveitis_session <- sessionInfo()
uveitis_final_checkpoint <- file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_withinstudy_final_checkpoint_2026-11-05.RData")
uveitis_save_objects <- setdiff(ls(pattern = "^uveitis_"), c("uveitis_precheck", "uveitis_loaded", "uveitis_save_objects"))
save(list = uveitis_save_objects, file = uveitis_final_checkpoint)
writeLines(capture.output(uveitis_session), file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_sessionInfo_2026-11-05.txt"))

# Record the input identity and check that all four requested images were created.
uveitis_input_manifest <- data.frame(path = uveitis_precheck_file, md5 = unname(tools::md5sum(uveitis_precheck_file)))
write.csv(uveitis_input_manifest, file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_input_manifest.csv"), row.names = FALSE)
uveitis_plot_files <- c(uveitis_ql_plot_file, uveitis_md_file, uveitis_heat_file, uveitis_dot_file)
stopifnot(all(file.exists(uveitis_plot_files)), all(file.info(uveitis_plot_files)$size > 0), file.exists(uveitis_final_checkpoint))
print(uveitis_plot_files)
print(file.info(uveitis_plot_files)$size)
