# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Glaucoma: observed expression of the nine DEGs across all eight libraries.
# Run in the current within-study R session, after the reviewed CPM block.
# The colour mapping follows the current group-coloured precheck boxplot.
# https://ggplot2.tidyverse.org/reference/facet_wrap.html

# Confirm that the plotting package and reviewed input objects are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  exists("glaucoma_deg_cpm"), exists("glaucoma_within_dge"),
  exists("glaucoma_within_metadata"), dir.exists(glaucoma_within_dir)
)

# Reshape the matrix into one row per gene-library measurement.
glaucoma_dot_data <- as.data.frame(
  as.table(glaucoma_deg_cpm), responseName = "CPM", stringsAsFactors = FALSE
)
names(glaucoma_dot_data)[1:2] <- c("gene_label", "sample_id")

# Match each measurement to the correct sample information and native gene symbol.
glaucoma_dot_sample_index <- match(
  glaucoma_dot_data$sample_id, glaucoma_within_metadata$sample_id
)
glaucoma_dot_gene_index <- match(
  glaucoma_dot_data$gene_label, rownames(glaucoma_within_dge$counts)
)
stopifnot(!anyNA(glaucoma_dot_sample_index), !anyNA(glaucoma_dot_gene_index))
glaucoma_dot_data$group <- glaucoma_within_metadata$group[glaucoma_dot_sample_index]
glaucoma_dot_data$sex <- glaucoma_within_metadata$sex[glaucoma_dot_sample_index]
glaucoma_dot_data$gene <- glaucoma_within_dge$genes$gene_symbol[glaucoma_dot_gene_index]
glaucoma_dot_data$gene <- factor(
  glaucoma_dot_data$gene, levels = sort(unique(glaucoma_dot_data$gene))
)

# Use chronological group order and the established boxplot colours.
glaucoma_plot_colours <- c(
  "Naïve_0D" = "#0072B2", "Glaucoma_3D" = "#D55E00",
  "Glaucoma_2W" = "#009E73", "Glaucoma_6W" = "#CC79A7"
)
glaucoma_plot_labels <- c(
  "Naïve_0D" = "Naive", "Glaucoma_3D" = "3 days",
  "Glaucoma_2W" = "2 weeks", "Glaucoma_6W" = "6 weeks"
)

# Show every library, separating the two sexes slightly along the horizontal axis.
# Each gene has its own vertical range, starting at zero.
glaucoma_dotplot <- ggplot2::ggplot(
  glaucoma_dot_data,
  ggplot2::aes(x = group, y = CPM, colour = group, shape = sex, group = sex)
) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.45), size = 2.8
  ) +
  ggplot2::facet_wrap(~ gene, ncol = 3, scales = "free_y") +
  ggplot2::scale_colour_manual(values = glaucoma_plot_colours, guide = "none") +
  ggplot2::scale_shape_manual(values = c(Female = 16, Male = 17)) +
  ggplot2::scale_x_discrete(labels = glaucoma_plot_labels) +
  ggplot2::scale_y_continuous(
    limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.03, 0.12))
  ) +
  ggplot2::labs(
    title = "Glaucoma: individual-library expression of DEGs",
    subtitle = "FDR < 0.05 at one or more stages | one point per pooled library",
    x = NULL, y = "TMM-normalized CPM", shape = "Sex",
    caption = paste(
      "Each panel has its own y-axis range. Points show observed expression;",
      "the statistical tests account for sex."
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(face = "italic", size = 11),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )
print(glaucoma_dotplot)

# Save the figure and the exact measurements and labels used to draw its points.
glaucoma_dotplot_file <- file.path(
  glaucoma_within_dir, "Glaucoma_GSE241782_DEG_expression_dotplots.png"
)
glaucoma_dotplot_data_file <- file.path(
  glaucoma_within_dir, "Glaucoma_GSE241782_DEG_dotplot_data.csv"
)
ggplot2::ggsave(
  glaucoma_dotplot_file, plot = glaucoma_dotplot,
  width = 11, height = 9, units = "in", dpi = 300, bg = "white"
)
write.csv(glaucoma_dot_data, glaucoma_dotplot_data_file, row.names = FALSE)
print(file.exists(glaucoma_dotplot_file))
