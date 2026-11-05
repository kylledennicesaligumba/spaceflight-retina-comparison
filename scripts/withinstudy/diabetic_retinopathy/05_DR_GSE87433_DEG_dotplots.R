# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# DR: observed expression of all three DEGs in individual selected arrays.
# Use the existing processed expression values and FDR-qualified DEG list.
# One point per array; no model refitting, retesting or additional normalization.

# Confirm that the plotting package and reviewed input objects are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  exists("dr_degs"), exists("dr_primary_expr"),
  exists("dr_primary_metadata"), dir.exists(dr_within_dir)
)

# Include all three DEGs, in ascending P-value order, with their supplied IDs and symbols.
dr_dot_genes <- dr_degs[order(dr_degs$P.Value), , drop = FALSE]
stopifnot(
  nrow(dr_dot_genes) == 3L,
  anyDuplicated(dr_dot_genes$ID) == 0L,
  all(dr_dot_genes$ID %in% rownames(dr_primary_expr))
)

# Select the original processed expression measurements for these genes and arrays.
dr_dot_expr <- dr_primary_expr[
  as.character(dr_dot_genes$ID), dr_primary_metadata$geo_accession, drop = FALSE
]
stopifnot(all(is.finite(dr_dot_expr)))

# Reshape to one row per gene-array measurement, retaining gene IDs as keys.
dr_dot_data <- as.data.frame(
  as.table(dr_dot_expr), responseName = "Expression", stringsAsFactors = FALSE
)
names(dr_dot_data)[1:2] <- c("gene_id", "sample_id")
dr_dot_sample_index <- match(dr_dot_data$sample_id, dr_primary_metadata$geo_accession)
dr_dot_gene_index <- match(dr_dot_data$gene_id, as.character(dr_dot_genes$ID))
stopifnot(!anyNA(dr_dot_sample_index), !anyNA(dr_dot_gene_index))
dr_dot_data$group <- dr_primary_metadata$group[dr_dot_sample_index]
dr_dot_data$Symbol <- dr_dot_genes$Symbol[dr_dot_gene_index]
dr_dot_data$gene <- factor(dr_dot_data$gene_id, levels = as.character(dr_dot_genes$ID))
dr_dot_gene_labels <- stats::setNames(dr_dot_genes$Symbol, dr_dot_genes$ID)

# Show all six arrays per group, spreading points horizontally without changing expression.
# Each gene has its own y-axis range; no zero baseline is forced on processed expression.
dr_dotplot <- ggplot2::ggplot(
  dr_dot_data, ggplot2::aes(x = group, y = Expression, colour = group)
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1),
    size = 2.8, alpha = 0.9
  ) +
  ggplot2::facet_wrap(
    ~ gene, ncol = 3, scales = "free_y",
    labeller = ggplot2::as_labeller(dr_dot_gene_labels)
  ) +
  ggplot2::scale_colour_manual(
    values = c(Normal = "#0072B2", Diabetic = "#D55E00"), guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.12, 0.16))
  ) +
  ggplot2::labs(
    title = "DR: expression of DEGs in individual samples",
    subtitle = "All 3 DEGs with FDR < 0.05 | Diabetic versus Normal",
    x = NULL, y = "Processed expression",
    caption = paste(
      "One point per array (six per group); horizontal offsets keep points visible.",
      "Each gene has its own y-axis range. Values are the supplied processed expression.",
      sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(face = "italic", size = 12),
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )
print(dr_dotplot)

# Save the dot plots and the exact measurements and labels used for their points.
dr_dotplot_file <- file.path(dr_within_dir, "DR_GSE87433_DEG_expression_dotplots.png")
dr_dotplot_data_file <- file.path(dr_within_dir, "DR_GSE87433_DEG_dotplot_data.csv")
ggplot2::ggsave(
  dr_dotplot_file, plot = dr_dotplot,
  width = 11, height = 5, units = "in", dpi = 300, bg = "white"
)
write.csv(dr_dot_data, dr_dotplot_data_file, row.names = FALSE)
print(file.exists(c(dr_dotplot_file, dr_dotplot_data_file)))
