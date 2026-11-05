# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# DR: relative expression of all three DEGs across the 12 selected arrays.
# Run in the current within-study session after reviewing the MD plot.
# This figure uses existing DEGs and processed microarray expression values.
# Row scaling: https://stat.ethz.ch/R-manual/R-devel/library/base/html/scale.html
# Plot options: https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf

# Confirm that the plotting packages and reviewed input objects are available.
stopifnot(
  requireNamespace("pheatmap", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("dr_degs"), exists("dr_primary_expr"),
  exists("dr_primary_metadata"), dir.exists(dr_within_dir)
)

# Include every DEG, placing higher genes first and ordering by P value within direction.
dr_heat_candidates <- dr_degs[order(dr_degs$P.Value), , drop = FALSE]
dr_heat_genes <- rbind(
  transform(subset(dr_heat_candidates, logFC > 0), Direction = "Higher"),
  transform(subset(dr_heat_candidates, logFC < 0), Direction = "Lower")
)
stopifnot(
  nrow(dr_heat_genes) == 3L,
  anyDuplicated(dr_heat_genes$ID) == 0L,
  all(dr_heat_genes$ID %in% rownames(dr_primary_expr))
)

# Place Normal arrays first, followed by Diabetic arrays, using GSM IDs within groups.
dr_heat_samples <- dr_primary_metadata[
  order(dr_primary_metadata$group, dr_primary_metadata$geo_accession),
  , drop = FALSE
]

# Select the three gene rows and 12 ordered arrays from the existing processed matrix.
# No extra logarithm, CPM calculation or expression normalization is introduced.
dr_heat_expr <- dr_primary_expr[
  as.character(dr_heat_genes$ID), dr_heat_samples$geo_accession, drop = FALSE
]
stopifnot(
  identical(colnames(dr_heat_expr), dr_heat_samples$geo_accession),
  all(is.finite(dr_heat_expr))
)

# Centre and scale each gene across arrays to show relative expression within that gene.
dr_heat_z <- t(scale(t(dr_heat_expr)))
stopifnot(all(is.finite(dr_heat_z)))
dr_heat_limit <- max(abs(dr_heat_z))

# Annotate columns by sample group and rows by the tested expression-change direction.
dr_heat_annotations <- data.frame(
  Group = dr_heat_samples$group,
  row.names = dr_heat_samples$geo_accession
)
dr_heat_row_annotations <- data.frame(
  Direction = factor(dr_heat_genes$Direction, levels = c("Higher", "Lower")),
  row.names = as.character(dr_heat_genes$ID)
)
dr_heat_annotation_colours <- list(
  Group = c(Normal = "#0072B2", Diabetic = "#D55E00"),
  Direction = c(Higher = "#D55E00", Lower = "#0072B2")
)

# Draw the familiar blue-white-red row-z-score heatmap, with fixed gene and sample order.
# Red is above a gene's own average; blue is below its own average across these arrays.
dr_heatmap <- pheatmap::pheatmap(
  dr_heat_z,
  color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-dr_heat_limit, dr_heat_limit, length.out = 102),
  scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = dr_heat_annotations, annotation_row = dr_heat_row_annotations,
  annotation_colors = dr_heat_annotation_colours,
  annotation_names_row = FALSE,
  labels_row = dr_heat_genes$Symbol,
  labels_col = dr_heat_samples$geo_accession,
  gaps_row = sum(dr_heat_genes$Direction == "Higher"),
  gaps_col = sum(dr_heat_samples$group == "Normal"),
  border_color = "white", fontsize = 10, fontsize_row = 12,
  fontsize_col = 9, angle_col = "45",
  main = paste(
    "DR DEGs: relative expression",
    "All 3 DEGs | Diabetic versus Normal",
    "Row z-scores of processed expression", sep = "\n"
  ),
  silent = TRUE
)
print(dr_heatmap)

# Set paths for the heatmap and the exact gene, expression, z-score and sample tables.
dr_heatmap_files <- file.path(dr_within_dir, c(
  "DR_GSE87433_DEG_expression_heatmap.png",
  "DR_GSE87433_DEG_heatmap_genes.csv",
  "DR_GSE87433_DEG_heatmap_processed_expression.csv",
  "DR_GSE87433_DEG_heatmap_row_zscores.csv",
  "DR_GSE87433_DEG_heatmap_sample_key.csv"
))

# Save the figure and supporting tables in the DR within-study folder.
ggplot2::ggsave(
  dr_heatmap_files[1], plot = dr_heatmap$gtable,
  width = 10, height = 5.5, units = "in", dpi = 300, bg = "white"
)
write.csv(dr_heat_genes, dr_heatmap_files[2], row.names = FALSE)
write.csv(dr_heat_expr, dr_heatmap_files[3], row.names = TRUE)
write.csv(dr_heat_z, dr_heatmap_files[4], row.names = TRUE)
write.csv(dr_heat_annotations, dr_heatmap_files[5], row.names = TRUE)
print(file.exists(dr_heatmap_files))
