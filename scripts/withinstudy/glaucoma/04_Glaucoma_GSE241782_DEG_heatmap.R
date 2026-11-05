# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Glaucoma: relative expression of the nine DEGs across the eight libraries.
# Run in the current within-study R session after reviewing the dot plot.
# Heatmap input follows the edgeR guide's moderated log-CPM recommendation.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the packages and reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("pheatmap", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("glaucoma_within_dge"), exists("glaucoma_deg_counts"),
  exists("glaucoma_within_metadata"), dir.exists(glaucoma_within_dir)
)

# Calculate observed log2 CPM with the saved TMM factors and a positive prior count.
# The prior count makes the logarithmic display finite when original counts are zero.
glaucoma_heat_log_cpm <- edgeR::cpm(
  glaucoma_within_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2
)[rownames(glaucoma_deg_counts), , drop = FALSE]

# Match the native gene names and arrange rows alphabetically for easy comparison.
glaucoma_heat_symbols <- glaucoma_within_dge$genes$gene_symbol[
  match(rownames(glaucoma_heat_log_cpm), rownames(glaucoma_within_dge$counts))
]
stopifnot(!anyNA(glaucoma_heat_symbols))
glaucoma_heat_row_order <- order(glaucoma_heat_symbols)
glaucoma_heat_log_cpm <- glaucoma_heat_log_cpm[glaucoma_heat_row_order, , drop = FALSE]
glaucoma_heat_symbols <- glaucoma_heat_symbols[glaucoma_heat_row_order]

# Keep controls first, followed by 3D/2W/6W, with Female then Male within each group.
glaucoma_heat_samples <- glaucoma_within_metadata[
  match(colnames(glaucoma_heat_log_cpm), glaucoma_within_metadata$sample_id),
  , drop = FALSE
]
glaucoma_heat_samples <- glaucoma_heat_samples[
  order(glaucoma_heat_samples$group, glaucoma_heat_samples$sex), , drop = FALSE
]
glaucoma_heat_log_cpm <- glaucoma_heat_log_cpm[
  , glaucoma_heat_samples$sample_id, drop = FALSE
]
stopifnot(identical(colnames(glaucoma_heat_log_cpm), glaucoma_heat_samples$sample_id))

# Centre and scale each gene across libraries so its own relative pattern sets the colours.
# This display transformation is separate from the already completed gene tests.
glaucoma_heat_z <- t(scale(t(glaucoma_heat_log_cpm)))
stopifnot(all(is.finite(glaucoma_heat_z)))
glaucoma_heat_limit <- max(abs(glaucoma_heat_z))

# Annotate each library by its stage and sex, retaining the chronological order.
glaucoma_heat_annotations <- data.frame(
  Stage = factor(
    glaucoma_heat_samples$group,
    levels = c("Naïve_0D", "Glaucoma_3D", "Glaucoma_2W", "Glaucoma_6W"),
    labels = c("Naive", "3 days", "2 weeks", "6 weeks")
  ),
  Sex = glaucoma_heat_samples$sex,
  row.names = glaucoma_heat_samples$sample_id
)
glaucoma_heat_annotation_colours <- list(
  Stage = c(
    "Naive" = "#0072B2", "3 days" = "#D55E00",
    "2 weeks" = "#009E73", "6 weeks" = "#CC79A7"
  ),
  Sex = c(Female = "#A6A6A6", Male = "#404040")
)

# Draw blue-to-red relative expression with a symmetric colour scale around zero.
# Fixed rows/columns preserve the same gene order and chronological sample order.
glaucoma_heatmap <- pheatmap::pheatmap(
  glaucoma_heat_z,
  color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-glaucoma_heat_limit, glaucoma_heat_limit, length.out = 102),
  scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = glaucoma_heat_annotations,
  annotation_colors = glaucoma_heat_annotation_colours,
  labels_row = glaucoma_heat_symbols,
  labels_col = paste(glaucoma_heat_samples$sample_id, substr(glaucoma_heat_samples$sex, 1, 1)),
  gaps_col = c(2, 4, 6), border_color = "white",
  fontsize = 10, fontsize_row = 11, fontsize_col = 10, angle_col = "45",
  main = "Glaucoma DEGs: relative expression\nRow z-scores of TMM log2 CPM",
  silent = TRUE
)
print(glaucoma_heatmap)

# Save the heatmap and its expression, colour-scale and annotation tables.
glaucoma_heatmap_file <- file.path(
  glaucoma_within_dir, "Glaucoma_GSE241782_DEG_expression_heatmap.png"
)
ggplot2::ggsave(
  glaucoma_heatmap_file, plot = glaucoma_heatmap$gtable,
  width = 10, height = 7, units = "in", dpi = 300, bg = "white"
)
write.csv(
  glaucoma_heat_log_cpm,
  file.path(glaucoma_within_dir, "Glaucoma_GSE241782_DEG_heatmap_logCPM.csv"),
  row.names = TRUE
)
write.csv(
  glaucoma_heat_z,
  file.path(glaucoma_within_dir, "Glaucoma_GSE241782_DEG_heatmap_row_zscores.csv"),
  row.names = TRUE
)
write.csv(
  glaucoma_heat_annotations,
  file.path(glaucoma_within_dir, "Glaucoma_GSE241782_DEG_heatmap_sample_key.csv"),
  row.names = TRUE
)
print(file.exists(glaucoma_heatmap_file))
