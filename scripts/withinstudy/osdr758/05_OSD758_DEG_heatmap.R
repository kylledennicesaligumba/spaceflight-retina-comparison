# OSD-758: relative expression of 15 higher and 15 lower DEGs across 17 libraries.
# Reuse the RP display selection; complete results and the DEG criterion stay unchanged.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the plotting packages and reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("pheatmap", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("osd758_within_dge"), exists("osd758_degs"),
  exists("osd758_within_metadata"), dir.exists(osd758_within_dir)
)

# Select the 15 lowest-P DEGs per direction, keeping higher genes first.
osd758_heat_candidates <- osd758_degs[order(osd758_degs$PValue), , drop = FALSE]
osd758_heat_genes <- rbind(
  transform(head(subset(osd758_heat_candidates, logFC > 0), n = 15), Direction = "Higher"),
  transform(head(subset(osd758_heat_candidates, logFC < 0), n = 15), Direction = "Lower")
)
stopifnot(
  nrow(osd758_heat_genes) == 30L,
  anyDuplicated(osd758_heat_genes$gene_id) == 0L,
  all(osd758_heat_genes$gene_id %in% rownames(osd758_within_dge$counts))
)

# Calculate log2 CPM from the full DGEList's saved TMM factors, then select rows.
# A positive prior count keeps the logarithmic display finite for zero counts.
osd758_heat_log_cpm <- edgeR::cpm(
  osd758_within_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2
)[osd758_heat_genes$gene_id, , drop = FALSE]

# Place HGC libraries first, then flight 0G, in sample-ID order within each group.
osd758_heat_samples <- osd758_within_metadata[
  match(colnames(osd758_heat_log_cpm), osd758_within_metadata$sample), , drop = FALSE
]
osd758_heat_samples <- osd758_heat_samples[
  order(osd758_heat_samples$group, osd758_heat_samples$sample), , drop = FALSE
]
osd758_heat_log_cpm <- osd758_heat_log_cpm[, osd758_heat_samples$sample, drop = FALSE]
stopifnot(identical(colnames(osd758_heat_log_cpm), osd758_heat_samples$sample))
osd758_heat_samples$short_label <- sub("^.*_", "", osd758_heat_samples$sample)

# Centre and scale each gene across all libraries to show relative expression.
osd758_heat_z <- t(scale(t(osd758_heat_log_cpm)))
stopifnot(all(is.finite(osd758_heat_z)))
osd758_heat_limit <- max(abs(osd758_heat_z))

# Annotate columns by group and rows by the tested change direction.
osd758_heat_annotations <- data.frame(
  Group = osd758_heat_samples$group,
  row.names = osd758_heat_samples$sample
)
osd758_heat_row_annotations <- data.frame(
  Direction = factor(osd758_heat_genes$Direction, levels = c("Higher", "Lower")),
  row.names = osd758_heat_genes$gene_id
)
osd758_heat_annotation_colours <- list(
  Group = c(HGC = "#0072B2", FLT_0G = "#D55E00"),
  Direction = c(Higher = "#D55E00", Lower = "#0072B2")
)

# Keep rows and columns in the declared order, using original Ensembl labels.
osd758_heatmap <- pheatmap::pheatmap(
  osd758_heat_z,
  color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-osd758_heat_limit, osd758_heat_limit, length.out = 102),
  scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = osd758_heat_annotations,
  annotation_row = osd758_heat_row_annotations,
  annotation_colors = osd758_heat_annotation_colours,
  annotation_names_row = FALSE,
  labels_row = osd758_heat_genes$gene_id,
  labels_col = osd758_heat_samples$short_label,
  gaps_row = 15, gaps_col = sum(osd758_heat_samples$group == "HGC"),
  border_color = "white", fontsize = 10, fontsize_row = 8,
  fontsize_col = 10, angle_col = "45",
  main = paste(
    "OSD-758: flight 0G versus HGC",
    "15 higher + 15 lower DEGs by lowest P value",
    "Row z-scores of TMM log2 CPM", sep = "\n"
  ),
  silent = TRUE
)
print(osd758_heatmap)

# Set paths for the figure, gene key, plotted matrices and full sample key.
osd758_heatmap_files <- file.path(
  osd758_within_dir,
  c(
    "OSD758_DEG_expression_heatmap_top30.png",
    "OSD758_top30_DEG_heatmap_selected_genes.csv",
    "OSD758_top30_DEG_heatmap_logCPM.csv",
    "OSD758_top30_DEG_heatmap_row_zscores.csv",
    "OSD758_top30_DEG_heatmap_sample_key.csv"
  )
)

# Save the heatmap and its exact plotted data in the within-study folder.
ggplot2::ggsave(
  osd758_heatmap_files[1], plot = osd758_heatmap$gtable,
  width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)
write.csv(osd758_heat_genes, osd758_heatmap_files[2], row.names = FALSE)
write.csv(osd758_heat_log_cpm, osd758_heatmap_files[3], row.names = TRUE)
write.csv(osd758_heat_z, osd758_heatmap_files[4], row.names = TRUE)
write.csv(osd758_heat_samples, osd758_heatmap_files[5], row.names = FALSE)
print(file.exists(osd758_heatmap_files))
