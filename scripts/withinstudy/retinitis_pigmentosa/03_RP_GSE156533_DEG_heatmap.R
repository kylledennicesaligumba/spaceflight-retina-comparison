# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# RP: relative expression of 15 higher and 15 lower DEGs across seven libraries.
# Gene selection is for this descriptive figure; full result tables remain available.
# Moderated log-CPM display follows the edgeR guide's heatmap recommendation.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the plotting packages and reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("pheatmap", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("rp_within_dge"), exists("rp_degs"),
  exists("rp_within_metadata"), dir.exists(rp_within_dir)
)

# Select the 15 DEGs with lowest P values in each direction, keeping higher genes first.
rp_heat_candidates <- rp_degs[order(rp_degs$PValue), , drop = FALSE]
rp_heat_genes <- rbind(
  transform(head(subset(rp_heat_candidates, logFC > 0), n = 15), Direction = "Higher"),
  transform(head(subset(rp_heat_candidates, logFC < 0), n = 15), Direction = "Lower")
)
stopifnot(
  nrow(rp_heat_genes) == 30L,
  anyDuplicated(rp_heat_genes$gene_id) == 0L,
  all(rp_heat_genes$gene_id %in% rownames(rp_within_dge$counts))
)

# Calculate log2 CPM using the full DGEList's saved TMM factors, then select these rows.
# The positive prior count makes the logarithmic display finite for zero counts.
rp_heat_log_cpm <- edgeR::cpm(
  rp_within_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2
)[rp_heat_genes$gene_id, , drop = FALSE]

# Place WT libraries first, then P23H, with Female before Male within each group.
rp_heat_samples <- rp_within_metadata[
  match(colnames(rp_heat_log_cpm), rp_within_metadata$sample_id), , drop = FALSE
]
rp_heat_samples <- rp_heat_samples[
  order(rp_heat_samples$group, rp_heat_samples$sex, rp_heat_samples$sample_id),
  , drop = FALSE
]
rp_heat_log_cpm <- rp_heat_log_cpm[, rp_heat_samples$sample_id, drop = FALSE]
stopifnot(identical(colnames(rp_heat_log_cpm), rp_heat_samples$sample_id))

# Centre and scale each gene across libraries to show its relative expression pattern.
rp_heat_z <- t(scale(t(rp_heat_log_cpm)))
stopifnot(all(is.finite(rp_heat_z)))
rp_heat_limit <- max(abs(rp_heat_z))

# Annotate columns by genotype and sex, and rows by the tested change direction.
rp_heat_annotations <- data.frame(
  Group = rp_heat_samples$group, Sex = rp_heat_samples$sex,
  row.names = rp_heat_samples$sample_id
)
rp_heat_row_annotations <- data.frame(
  Direction = factor(rp_heat_genes$Direction, levels = c("Higher", "Lower")),
  row.names = rp_heat_genes$gene_id
)
rp_heat_annotation_colours <- list(
  Group = c(WT = "#0072B2", P23H = "#D55E00"),
  Sex = c(Female = "#A6A6A6", Male = "#404040"),
  Direction = c(Higher = "#D55E00", Lower = "#0072B2")
)

# Draw the same blue-white-red row-z-score display used for glaucoma.
# Fixed rows and columns preserve direction groups and the WT/P23H sample order.
rp_heatmap <- pheatmap::pheatmap(
  rp_heat_z,
  color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-rp_heat_limit, rp_heat_limit, length.out = 102),
  scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = rp_heat_annotations, annotation_row = rp_heat_row_annotations,
  annotation_colors = rp_heat_annotation_colours,
  annotation_names_row = FALSE,
  labels_row = rp_heat_genes$gene_name,
  labels_col = paste(rp_heat_samples$sample_id, substr(rp_heat_samples$sex, 1, 1)),
  gaps_row = 15, gaps_col = sum(rp_heat_samples$group == "WT"),
  border_color = "white", fontsize = 10, fontsize_row = 9,
  fontsize_col = 10, angle_col = "45",
  main = paste(
    "RP DEGs: relative expression",
    "15 higher + 15 lower DEGs by lowest P value",
    "Row z-scores of TMM log2 CPM", sep = "\n"
  ),
  silent = TRUE
)
print(rp_heatmap)

# Set paths for the figure, selected-gene key and the plotted expression/annotation tables.
rp_heatmap_files <- file.path(
  rp_within_dir,
  c(
    "RP_GSE156533_DEG_expression_heatmap_top30.png",
    "RP_GSE156533_top30_DEG_heatmap_selected_genes.csv",
    "RP_GSE156533_top30_DEG_heatmap_logCPM.csv",
    "RP_GSE156533_top30_DEG_heatmap_row_zscores.csv",
    "RP_GSE156533_top30_DEG_heatmap_sample_key.csv"
  )
)

# Save the heatmap and its data in the RP within-study folder.
ggplot2::ggsave(
  rp_heatmap_files[1], plot = rp_heatmap$gtable,
  width = 10, height = 10, units = "in", dpi = 300, bg = "white"
)
write.csv(rp_heat_genes, file = rp_heatmap_files[2], row.names = FALSE)
write.csv(rp_heat_log_cpm, file = rp_heatmap_files[3], row.names = TRUE)
write.csv(rp_heat_z, file = rp_heatmap_files[4], row.names = TRUE)
write.csv(rp_heat_annotations, file = rp_heatmap_files[5], row.names = TRUE)
print(file.exists(rp_heatmap_files))
