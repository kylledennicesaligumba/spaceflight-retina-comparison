# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# LHON-linked model: selected combined-genotype DEGs across all 12 libraries.
# Select 15 higher and 15 lower DEGs by lowest P value in YARS2/COI versus WT.
# This figure selection uses the same rule as RP and preserves all full result tables.
# Method: https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html#heat-map-clustering

# Confirm that plotting packages and the reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("pheatmap", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("lhon_within_dge"), exists("lhon_YARS2_COI_degs"),
  exists("lhon_within_metadata"), dir.exists(lhon_within_dir)
)

# Select the 15 lowest-P-value DEGs in each direction, keeping higher genes first.
lhon_heat_candidates <- lhon_YARS2_COI_degs[
  order(lhon_YARS2_COI_degs$PValue), , drop = FALSE
]
lhon_heat_genes <- rbind(
  transform(head(subset(lhon_heat_candidates, logFC > 0), n = 15), Direction = "Higher"),
  transform(head(subset(lhon_heat_candidates, logFC < 0), n = 15), Direction = "Lower")
)

# LHON's supplied unique gene_name labels are the count-row keys.
stopifnot(
  nrow(lhon_heat_genes) == 30L,
  anyDuplicated(lhon_heat_genes$gene_name) == 0L,
  all(lhon_heat_genes$gene_name %in% rownames(lhon_within_dge$counts))
)

# Calculate log2 CPM using the full DGEList's saved TMM factors, then select these genes.
# The positive prior count makes the logarithmic display finite for zero counts.
lhon_heat_log_cpm <- edgeR::cpm(
  lhon_within_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2
)[lhon_heat_genes$gene_name, , drop = FALSE]

# Order all libraries as WT, YARS2, COI and YARS2/COI, with sample IDs within group.
lhon_heat_samples <- lhon_within_metadata[
  match(colnames(lhon_heat_log_cpm), lhon_within_metadata$sample_id), , drop = FALSE
]
lhon_heat_samples <- lhon_heat_samples[
  order(lhon_heat_samples$group, lhon_heat_samples$sample_id), , drop = FALSE
]
lhon_heat_log_cpm <- lhon_heat_log_cpm[, lhon_heat_samples$sample_id, drop = FALSE]
stopifnot(identical(colnames(lhon_heat_log_cpm), lhon_heat_samples$sample_id))

# Centre and scale each gene across the 12 libraries to show relative expression.
lhon_heat_z <- t(scale(t(lhon_heat_log_cpm)))
stopifnot(all(is.finite(lhon_heat_z)))
lhon_heat_limit <- max(abs(lhon_heat_z))

# Annotate columns by genotype and rows by direction in the combined-versus-WT test.
# All libraries are male, so genotype provides the varying sample annotation.
lhon_heat_annotations <- data.frame(
  Group = lhon_heat_samples$group, row.names = lhon_heat_samples$sample_id
)
lhon_heat_row_annotations <- data.frame(
  Direction = factor(lhon_heat_genes$Direction, levels = c("Higher", "Lower")),
  row.names = lhon_heat_genes$gene_name
)
lhon_heat_annotation_colours <- list(
  Group = c(
    WT = "#0072B2", YARS2 = "#D55E00", COI = "#009E73", "YARS2/COI" = "#CC79A7"
  ),
  Direction = c(Higher = "#D55E00", Lower = "#0072B2")
)

# Draw the familiar blue-white-red row-z-score display, with fixed gene/sample order.
# Gaps separate higher/lower genes and the four groups of three libraries.
lhon_heatmap <- pheatmap::pheatmap(
  lhon_heat_z,
  color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-lhon_heat_limit, lhon_heat_limit, length.out = 102),
  scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = lhon_heat_annotations, annotation_row = lhon_heat_row_annotations,
  annotation_colors = lhon_heat_annotation_colours,
  annotation_names_row = FALSE,
  labels_row = lhon_heat_genes$gene_name,
  labels_col = lhon_heat_samples$sample_id,
  gaps_row = 15, gaps_col = c(3, 6, 9),
  border_color = "white", fontsize = 10, fontsize_row = 9,
  fontsize_col = 10, angle_col = "45",
  main = paste(
    "LHON-linked model: selected DEG expression",
    "15 higher + 15 lower DEGs in YARS2/COI vs WT, by lowest P value",
    "Row z-scores of TMM log2 CPM | all 12 male libraries", sep = "\n"
  ),
  silent = TRUE
)
print(lhon_heatmap)

# Set paths for the figure, selected-gene results and plotted expression/sample tables.
lhon_heatmap_files <- file.path(
  lhon_within_dir,
  c(
    "LHON_GSE304707_DEG_expression_heatmap_top30.png",
    "LHON_GSE304707_top30_DEG_heatmap_selected_genes.csv",
    "LHON_GSE304707_top30_DEG_heatmap_logCPM.csv",
    "LHON_GSE304707_top30_DEG_heatmap_row_zscores.csv",
    "LHON_GSE304707_top30_DEG_heatmap_sample_key.csv"
  )
)

# Save the heatmap and its underlying data in the LHON output folder.
ggplot2::ggsave(
  lhon_heatmap_files[1], plot = lhon_heatmap$gtable,
  width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)
write.csv(lhon_heat_genes, file = lhon_heatmap_files[2], row.names = FALSE)
write.csv(lhon_heat_log_cpm, file = lhon_heatmap_files[3], row.names = TRUE)
write.csv(lhon_heat_z, file = lhon_heatmap_files[4], row.names = TRUE)
write.csv(lhon_heat_samples, file = lhon_heatmap_files[5], row.names = FALSE)
print(file.exists(lhon_heatmap_files))
