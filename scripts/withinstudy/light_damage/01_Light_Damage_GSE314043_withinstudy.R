# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "light_damage"), recursive = TRUE, showWarnings = FALSE)

# Light damage GSE314043: within-study gene-level analysis, 5 November 2026.
# Continue from the reviewed filtered counts and saved TMM factors.
# Method: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Locate the output folder and completed precheck without using the live R session.
ld_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "light_damage")
dir.create(ld_within_dir, recursive = TRUE, showWarnings = FALSE)
ld_precheck_file <- Sys.getenv("LIGHT_DAMAGE_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "light_damage", "Light_Damage_GSE314043_precheck_checkpoint_2026-11-03.RData"))
stopifnot(file.exists(ld_precheck_file))
ld_precheck <- new.env()
ld_loaded <- load(ld_precheck_file, envir = ld_precheck)
ld_within_dge <- ld_precheck$ld_dge_qc
ld_within_metadata <- ld_precheck$ld_metadata_selected

# Check the reviewed dimensions, sample order, gene order and saved scaling factors.
stopifnot(
  identical(dim(ld_within_dge), c(16852L, 10L)),
  identical(colnames(ld_within_dge), ld_within_metadata$sample_id),
  identical(rownames(ld_within_dge), ld_within_dge$genes$Geneid),
  anyDuplicated(rownames(ld_within_dge)) == 0L,
  all(is.finite(ld_within_dge$samples$norm.factors)),
  all(ld_within_dge$samples$norm.factors > 0)
)
print(dim(ld_within_dge))

# Set control as reference; matching replicate numbers do not establish pairing.
# Sex is unreported in the reviewed metadata, so no sex category is inferred.
ld_within_metadata$group <- factor(
  ld_within_metadata$condition, levels = c("CTRL", "LD")
)
stopifnot(!anyNA(ld_within_metadata$group))
print(table(ld_within_metadata$group))
ld_design <- model.matrix(~ group, data = ld_within_metadata)
rownames(ld_design) <- ld_within_metadata$sample_id
stopifnot(qr(ld_design)$rank == ncol(ld_design))
print(ld_design)
print(nrow(ld_design) - qr(ld_design)$rank)
write.csv(ld_within_metadata, file.path(ld_within_dir, "Light_Damage_GSE314043_sample_key.csv"), row.names = FALSE)
write.csv(ld_design, file.path(ld_within_dir, "Light_Damage_GSE314043_design.csv"))

# Fit the same robust current edgeR QL pipeline used for the prior count datasets.
# Raw filtered counts and their saved TMM offsets enter the model, not log-CPM.
stopifnot(as.character(utils::packageVersion("edgeR")) == "4.8.2")
ld_fit <- edgeR::glmQLFit(ld_within_dge, design = ld_design, robust = TRUE)
stopifnot(identical(dim(ld_fit$coefficients), c(16852L, 2L)))

# Save the QL dispersion diagnostic: raw estimates, moderated estimates and trend.
ld_ql_plot_file <- file.path(ld_within_dir, "Light_Damage_GSE314043_QL_dispersion.png")
png(ld_ql_plot_file, width = 2400, height = 1800, res = 300)
edgeR::plotQLDisp(ld_fit, main = "Light damage: QL dispersion")
dev.off()

# Test disease/injury versus control and retain all tested genes with BH FDR.
ld_test <- edgeR::glmQLFTest(ld_fit, coef = "groupLD")
ld_results <- edgeR::topTags(ld_test, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
ld_results$gene_id <- rownames(ld_results)
stopifnot(
  identical(ld_results$gene_id, as.character(ld_results$Geneid)),
  all(is.finite(ld_results$logFC)), all(is.finite(ld_results$F)),
  all(is.finite(ld_results$PValue)), all(is.finite(ld_results$FDR)),
  all(ld_results$FDR >= 0 & ld_results$FDR <= 1)
)

# Use an available name for display; keep original source IDs as unique row keys.
# Fallback labels are source identifiers, not newly verified gene-symbol mappings.
ld_results$display_name <- as.character(ld_results$gene_name)
ld_missing_name <- is.na(ld_results$display_name) | !nzchar(ld_results$display_name)
ld_results$display_name[ld_missing_name] <- ld_results$gene_id[ld_missing_name]
ld_duplicate_name <- duplicated(ld_results$display_name) | duplicated(ld_results$display_name, fromLast = TRUE)
ld_results$display_name[ld_duplicate_name] <- paste0(
  ld_results$display_name[ld_duplicate_name], " [", ld_results$gene_id[ld_duplicate_name], "]"
)

# Apply the established FDR < 0.05 criterion without an extra fold-change cutoff.
ld_degs <- subset(ld_results, FDR < 0.05)
ld_summary <- data.frame(
  Comparison = "LD_vs_CTRL", Tested = nrow(ld_results),
  DEGs = nrow(ld_degs), Higher = sum(ld_degs$logFC > 0),
  Lower = sum(ld_degs$logFC < 0)
)
stopifnot(ld_summary$Higher + ld_summary$Lower == ld_summary$DEGs)
print(ld_summary)
print(head(ld_results[, c("gene_id", "display_name", "logFC", "F", "PValue", "FDR")], 10))
write.csv(ld_results, file.path(ld_within_dir, "Light_Damage_GSE314043_LD_vs_CTRL_all_genes.csv"), row.names = FALSE)
write.csv(ld_degs, file.path(ld_within_dir, "Light_Damage_GSE314043_LD_vs_CTRL_DEGs_FDR005.csv"), row.names = FALSE)
write.csv(ld_summary, file.path(ld_within_dir, "Light_Damage_GSE314043_DEG_summary.csv"), row.names = FALSE)

# Rank every tested gene using the established signed square-root F score.
# The tolerance handles only negligible numerical negatives and leaves original F intact.
ld_rank_F_check <- c(Minimum_F = min(ld_results$F), Negative_F = sum(ld_results$F < 0), Nonfinite_F = sum(!is.finite(ld_results$F)))
ld_rank_tolerance <- sqrt(.Machine$double.eps)
stopifnot(all(ld_results$F >= -ld_rank_tolerance))
ld_ranked <- transform(ld_results, rank_score = sign(logFC) * sqrt(pmax(F, 0)))
stopifnot(all(is.finite(ld_ranked$rank_score)))
ld_ranked <- ld_ranked[order(ld_ranked$rank_score, decreasing = TRUE), , drop = FALSE]
print(ld_rank_F_check)
write.csv(ld_ranked, file.path(ld_within_dir, "Light_Damage_GSE314043_LD_vs_CTRL_ranked_all_genes.csv"), row.names = FALSE)

# Plot all tested genes; colour significant increases/decreases and label up to five each.
ld_md_data <- ld_results
ld_md_data$Status <- "Not significant"
ld_md_data$Status[ld_md_data$FDR < 0.05 & ld_md_data$logFC > 0] <- "Higher"
ld_md_data$Status[ld_md_data$FDR < 0.05 & ld_md_data$logFC < 0] <- "Lower"
ld_md_hits <- subset(ld_md_data, Status != "Not significant")
ld_md_labels <- rbind(head(subset(ld_md_hits, Status == "Higher"), 5), head(subset(ld_md_hits, Status == "Lower"), 5))
ld_md_plot <- ggplot2::ggplot(ld_md_data, ggplot2::aes(logCPM, logFC)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(data = subset(ld_md_data, Status == "Not significant"), colour = "grey65", alpha = 0.35, size = 0.6) +
  ggplot2::geom_point(data = ld_md_hits, ggplot2::aes(colour = Status), size = 1, alpha = 0.7) +
  ggrepel::geom_text_repel(data = ld_md_labels, ggplot2::aes(label = display_name, colour = Status),
    size = 3, seed = 1, box.padding = 0.45, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 5, max.iter = 20000,
    nudge_y = ifelse(ld_md_labels$logFC > 0, 1.2, -0.4), force_pull = 0.05, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(Higher = "#D55E00", Lower = "#0072B2"), breaks = c("Higher", "Lower")) +
  ggplot2::labs(title = "Light damage: genome-wide differential expression", subtitle = "LD versus CTRL | all retained libraries",
    x = "Average log2 CPM (all 10 libraries)", y = "Log2 fold change (LD versus CTRL)", colour = "FDR < 0.05",
    caption = "All tested genes shown. Labels: up to five lowest-P-value DEGs per direction.\nGrey: FDR criterion not met. No separate fold-change cutoff.") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5), plot.caption = ggplot2::element_text(hjust = 0))
ld_md_file <- file.path(ld_within_dir, "Light_Damage_GSE314043_MD_plot.png")
ggplot2::ggsave(ld_md_file, ld_md_plot, width = 9, height = 6.5, units = "in", dpi = 300, bg = "white")

# Select up to 15 significant genes per direction by P value for the descriptive heatmap.
# If a direction has fewer significant genes, use all available; do not relabel non-DEGs.
ld_heat_genes <- rbind(
  transform(head(subset(ld_degs, logFC > 0), 15), Direction = rep("Higher", min(15L, sum(ld_degs$logFC > 0)))),
  transform(head(subset(ld_degs, logFC < 0), 15), Direction = rep("Lower", min(15L, sum(ld_degs$logFC < 0))))
)
ld_heat_file <- file.path(ld_within_dir, "Light_Damage_GSE314043_DEG_expression_heatmap_top30.png")
ld_heat_samples <- ld_within_metadata[order(ld_within_metadata$group, ld_within_metadata$sample_id), , drop = FALSE]
if (nrow(ld_heat_genes) > 0L) {
  # Calculate normalized log-CPM from the full filtered object before selecting rows.
  ld_heat_log_cpm <- edgeR::cpm(ld_within_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)[ld_heat_genes$gene_id, ld_heat_samples$sample_id, drop = FALSE]
  ld_heat_z <- t(scale(t(ld_heat_log_cpm)))
  stopifnot(all(is.finite(ld_heat_z)))
  ld_heat_limit <- max(abs(ld_heat_z))
  ld_heat_annotations <- data.frame(Group = ld_heat_samples$group, row.names = ld_heat_samples$sample_id)
  ld_heat_row_annotations <- data.frame(Direction = factor(ld_heat_genes$Direction, levels = c("Higher", "Lower")), row.names = ld_heat_genes$gene_id)
  ld_heat_annotation_colours <- list(Group = c(CTRL = "#0072B2", LD = "#D55E00"), Direction = c(Higher = "#D55E00", Lower = "#0072B2"))
  ld_heat_nhigher <- sum(ld_heat_genes$Direction == "Higher")
  ld_heatmap <- pheatmap::pheatmap(ld_heat_z,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    breaks = seq(-ld_heat_limit, ld_heat_limit, length.out = 102), scale = "none",
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = ld_heat_annotations, annotation_row = ld_heat_row_annotations,
    annotation_colors = ld_heat_annotation_colours, annotation_names_row = FALSE,
    labels_row = ld_heat_genes$display_name, labels_col = ld_heat_samples$sample_id,
    gaps_row = if (ld_heat_nhigher > 0L && ld_heat_nhigher < nrow(ld_heat_genes)) ld_heat_nhigher else NULL,
    gaps_col = sum(ld_heat_samples$group == "CTRL"),
    border_color = "white", fontsize = 10, fontsize_row = 9, fontsize_col = 9, angle_col = "45",
    main = paste("Light damage DEGs: relative expression", paste0(ld_heat_nhigher, " higher + ", sum(ld_heat_genes$Direction == "Lower"), " lower by lowest P value"), "Row z-scores of TMM log2 CPM", sep = "\n"), silent = TRUE)
  ggplot2::ggsave(ld_heat_file, ld_heatmap$gtable, width = 11, height = 10, units = "in", dpi = 300, bg = "white")
  write.csv(ld_heat_log_cpm, file.path(ld_within_dir, "Light_Damage_GSE314043_top30_DEG_heatmap_logCPM.csv"))
  write.csv(ld_heat_z, file.path(ld_within_dir, "Light_Damage_GSE314043_top30_DEG_heatmap_row_zscores.csv"))
  write.csv(ld_heat_annotations, file.path(ld_within_dir, "Light_Damage_GSE314043_top30_DEG_heatmap_sample_key.csv"))
} else {
  # Preserve a truthful figure if no gene meets the established DEG criterion.
  png(ld_heat_file, width = 2400, height = 1800, res = 300)
  plot.new(); text(0.5, 0.5, "No DEGs at FDR < 0.05; no DEG heatmap can be drawn.")
  dev.off()
}
write.csv(ld_heat_genes, file.path(ld_within_dir, "Light_Damage_GSE314043_top30_DEG_heatmap_selected_genes.csv"), row.names = FALSE)

# Select up to five DEGs per direction, using the same P-value ordering as MD labels.
ld_dot_genes <- rbind(
  transform(head(subset(ld_degs, logFC > 0), 5), Direction = rep("Higher", min(5L, sum(ld_degs$logFC > 0)))),
  transform(head(subset(ld_degs, logFC < 0), 5), Direction = rep("Lower", min(5L, sum(ld_degs$logFC < 0))))
)
ld_dot_file <- file.path(ld_within_dir, "Light_Damage_GSE314043_DEG_expression_dotplots_top10.png")
if (nrow(ld_dot_genes) > 0L) {
  # Display observed non-log CPM, retaining one point per measured library.
  ld_dot_cpm <- edgeR::cpm(ld_within_dge, normalized.lib.sizes = TRUE, log = FALSE)[ld_dot_genes$gene_id, , drop = FALSE]
  ld_dot_data <- as.data.frame(as.table(ld_dot_cpm), responseName = "CPM", stringsAsFactors = FALSE)
  names(ld_dot_data)[1:2] <- c("gene_id", "sample_id")
  ld_dot_sample_index <- match(ld_dot_data$sample_id, ld_within_metadata$sample_id)
  ld_dot_gene_index <- match(ld_dot_data$gene_id, ld_dot_genes$gene_id)
  stopifnot(!anyNA(ld_dot_sample_index), !anyNA(ld_dot_gene_index))
  ld_dot_data$group <- ld_within_metadata$group[ld_dot_sample_index]
  ld_dot_data$display_name <- ld_dot_genes$display_name[ld_dot_gene_index]
  ld_dot_data$Direction <- ld_dot_genes$Direction[ld_dot_gene_index]
  ld_dot_data$gene <- factor(ld_dot_data$gene_id, levels = ld_dot_genes$gene_id)
  ld_dot_labels <- stats::setNames(ld_dot_genes$display_name, ld_dot_genes$gene_id)
  ld_dotplot <- ggplot2::ggplot(ld_dot_data, ggplot2::aes(group, CPM, colour = group)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1), size = 2.8) +
    ggplot2::facet_wrap(~ gene, ncol = 2, dir = "v", scales = "free_y", labeller = ggplot2::as_labeller(ld_dot_labels)) +
    ggplot2::scale_colour_manual(values = c(CTRL = "#0072B2", LD = "#D55E00"), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
    ggplot2::labs(title = "Light damage: individual-library expression of selected DEGs",
      subtitle = "Up to five higher DEGs on left; up to five lower on right\nLowest P values within direction among genes with FDR < 0.05",
      x = NULL, y = "TMM-normalized CPM",
      caption = "One point per library; horizontal offsets keep points visible. Each gene has its own y-axis range.\nPoints show observed expression. Sex was not recorded in the reviewed metadata.") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
      strip.text = ggplot2::element_text(face = "italic", size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5), plot.caption = ggplot2::element_text(hjust = 0))
  ggplot2::ggsave(ld_dot_file, ld_dotplot, width = 11, height = 12, units = "in", dpi = 300, bg = "white")
  write.csv(ld_dot_data, file.path(ld_within_dir, "Light_Damage_GSE314043_DEG_dotplot_top10_data.csv"), row.names = FALSE)
} else {
  # Do not present nonsignificant genes as significant in an empty-result analysis.
  png(ld_dot_file, width = 2400, height = 1800, res = 300)
  plot.new(); text(0.5, 0.5, "No DEGs at FDR < 0.05; no DEG dot plots can be drawn.")
  dev.off()
}
write.csv(ld_dot_genes, file.path(ld_within_dir, "Light_Damage_GSE314043_top10_DEG_dotplot_selected_genes.csv"), row.names = FALSE)

# Save only this study's working objects, excluding the precheck environment.
ld_session <- sessionInfo()
ld_final_checkpoint <- file.path(ld_within_dir, "Light_Damage_GSE314043_withinstudy_final_checkpoint_2026-11-05.RData")
ld_save_objects <- setdiff(ls(pattern = "^ld_"), c("ld_precheck", "ld_loaded", "ld_save_objects"))
save(list = ld_save_objects, file = ld_final_checkpoint)
writeLines(capture.output(ld_session), file.path(ld_within_dir, "Light_Damage_GSE314043_sessionInfo_2026-11-05.txt"))

# Record the input identity and check that all four requested images were created.
ld_input_manifest <- data.frame(path = ld_precheck_file, md5 = unname(tools::md5sum(ld_precheck_file)))
write.csv(ld_input_manifest, file.path(ld_within_dir, "Light_Damage_GSE314043_input_manifest.csv"), row.names = FALSE)
ld_plot_files <- c(ld_ql_plot_file, ld_md_file, ld_heat_file, ld_dot_file)
stopifnot(all(file.exists(ld_plot_files)), all(file.info(ld_plot_files)$size > 0), file.exists(ld_final_checkpoint))
print(ld_plot_files)
print(file.info(ld_plot_files)$size)
