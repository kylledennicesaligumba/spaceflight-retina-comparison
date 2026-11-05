# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# RP: observed expression of five higher and five lower DEGs in all seven libraries.
# Select genes by lowest P value within direction among the existing FDR-qualified DEGs.
# This display selection does not change the complete gene results or DEG criterion.

# Confirm that the plotting package and reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("rp_degs"), exists("rp_within_dge"),
  exists("rp_within_metadata"), dir.exists(rp_within_dir)
)

# Select five DEGs per direction in the same order used for the MD labels.
rp_dot_candidates <- rp_degs[order(rp_degs$PValue), , drop = FALSE]
rp_dot_genes <- rbind(
  transform(head(subset(rp_dot_candidates, logFC > 0), n = 5), Direction = "Higher"),
  transform(head(subset(rp_dot_candidates, logFC < 0), n = 5), Direction = "Lower")
)
stopifnot(
  nrow(rp_dot_genes) == 10L,
  anyDuplicated(rp_dot_genes$gene_id) == 0L,
  all(rp_dot_genes$gene_id %in% rownames(rp_within_dge$counts))
)

# Calculate CPM from the full DGEList with its saved TMM factors, then select these genes.
rp_dot_cpm <- edgeR::cpm(
  rp_within_dge, normalized.lib.sizes = TRUE, log = FALSE
)[rp_dot_genes$gene_id, , drop = FALSE]

# Reshape to one row per gene-library measurement, retaining gene IDs as keys.
rp_dot_data <- as.data.frame(
  as.table(rp_dot_cpm), responseName = "CPM", stringsAsFactors = FALSE
)
names(rp_dot_data)[1:2] <- c("gene_id", "sample_id")
rp_dot_sample_index <- match(rp_dot_data$sample_id, rp_within_metadata$sample_id)
rp_dot_gene_index <- match(rp_dot_data$gene_id, rp_dot_genes$gene_id)
stopifnot(!anyNA(rp_dot_sample_index), !anyNA(rp_dot_gene_index))
rp_dot_data$group <- rp_within_metadata$group[rp_dot_sample_index]
rp_dot_data$sex <- rp_within_metadata$sex[rp_dot_sample_index]
rp_dot_data$gene_name <- rp_dot_genes$gene_name[rp_dot_gene_index]
rp_dot_data$Direction <- rp_dot_genes$Direction[rp_dot_gene_index]
rp_dot_data$gene <- factor(rp_dot_data$gene_id, levels = rp_dot_genes$gene_id)
rp_dot_gene_labels <- stats::setNames(rp_dot_genes$gene_name, rp_dot_genes$gene_id)

# Show each library, separating sexes and slightly spreading same-sex replicates horizontally.
# Vertical CPM values are unchanged; each gene has its own y-axis range.
rp_dotplot <- ggplot2::ggplot(
  rp_dot_data,
  ggplot2::aes(x = group, y = CPM, colour = group, shape = sex, group = sex)
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.12, jitter.height = 0, dodge.width = 0.45, seed = 1
    ),
    size = 2.8
  ) +
  ggplot2::facet_wrap(
    ~ gene, ncol = 2, dir = "v", scales = "free_y",
    labeller = ggplot2::as_labeller(rp_dot_gene_labels)
  ) +
  ggplot2::scale_colour_manual(
    values = c(WT = "#0072B2", P23H = "#D55E00"), guide = "none"
  ) +
  ggplot2::scale_shape_manual(values = c(Female = 16, Male = 17)) +
  ggplot2::scale_y_continuous(
    limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.03, 0.12))
  ) +
  ggplot2::labs(
    title = "RP: individual-library expression of selected DEGs",
    subtitle = paste(
      "Left: five higher DEGs | Right: five lower DEGs",
      "Lowest P values within direction among genes with FDR < 0.05", sep = "\n"
    ),
    x = NULL, y = "TMM-normalized CPM", shape = "Sex",
    caption = paste(
      "Each panel has its own y-axis range. One point per library; horizontal offsets keep points visible.",
      "Points show observed expression; the statistical test accounts for sex.", sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(face = "italic", size = 11),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )
print(rp_dotplot)

# Save the figure and the exact measurements and labels used for its points.
rp_dotplot_file <- file.path(
  rp_within_dir, "RP_GSE156533_DEG_expression_dotplots_top10.png"
)
rp_dotplot_data_file <- file.path(
  rp_within_dir, "RP_GSE156533_DEG_dotplot_top10_data.csv"
)
ggplot2::ggsave(
  rp_dotplot_file, plot = rp_dotplot,
  width = 11, height = 12, units = "in", dpi = 300, bg = "white"
)
write.csv(rp_dot_data, file = rp_dotplot_data_file, row.names = FALSE)
print(file.exists(c(rp_dotplot_file, rp_dotplot_data_file)))
