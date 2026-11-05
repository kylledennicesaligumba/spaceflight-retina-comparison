# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Glaucoma: three genome-wide mean-difference (MD/MA) plots.
# Run in the current within-study R session after reviewing the heatmap.
# All tested genes are included; significance uses the existing BH FDR < 0.05 rule.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the plotting packages and complete result tables are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE),
  exists("glaucoma_3D_results"), exists("glaucoma_2W_results"),
  exists("glaucoma_6W_results"), dir.exists(glaucoma_within_dir)
)

# Stack the three full tables for plotting, retaining a separate panel for each contrast.
glaucoma_md_data <- rbind(
  transform(glaucoma_3D_results, Stage = "3D"),
  transform(glaucoma_2W_results, Stage = "2W"),
  transform(glaucoma_6W_results, Stage = "6W")
)
glaucoma_md_data$Stage <- factor(
  glaucoma_md_data$Stage,
  levels = c("3D", "2W", "6W"),
  labels = c("3 days vs naive", "2 weeks vs naive", "6 weeks vs naive")
)

# Stop for review if any required plotted statistic is missing or nonfinite.
stopifnot(
  all(is.finite(glaucoma_md_data$logFC)),
  all(is.finite(glaucoma_md_data$logCPM)),
  all(is.finite(glaucoma_md_data$FDR))
)

# Use the same FDR threshold and change directions already reviewed in the result tables.
glaucoma_md_data$Status <- "Not significant"
glaucoma_md_data$Status[
  glaucoma_md_data$FDR < 0.05 & glaucoma_md_data$logFC > 0
] <- "Higher"
glaucoma_md_data$Status[
  glaucoma_md_data$FDR < 0.05 & glaucoma_md_data$logFC < 0
] <- "Lower"
glaucoma_md_hits <- subset(glaucoma_md_data, Status != "Not significant")

# Plot the full background first, then highlight and label significant genes.
# Fixed axis ranges make the three contrasts directly readable on the same scales.
glaucoma_md_plot <- ggplot2::ggplot(
  glaucoma_md_data, ggplot2::aes(x = logCPM, y = logFC)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(
    data = subset(glaucoma_md_data, Status == "Not significant"),
    colour = "grey65", alpha = 0.35, size = 0.6
  ) +
  ggplot2::geom_point(
    data = glaucoma_md_hits, ggplot2::aes(colour = Status), size = 2
  ) +
  ggrepel::geom_text_repel(
    data = glaucoma_md_hits,
    ggplot2::aes(label = gene_symbol, colour = Status),
    size = 3, seed = 1, box.padding = 0.45, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 2,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Stage, nrow = 1) +
  ggplot2::scale_colour_manual(values = c(Higher = "#D55E00", Lower = "#0072B2")) +
  ggplot2::labs(
    title = "Glaucoma: genome-wide differential expression",
    subtitle = "Sex-adjusted stage-versus-naive comparisons | all tested genes shown",
    x = "Average log2 CPM (all eight libraries)",
    y = "Log2 fold change versus naive", colour = "FDR < 0.05",
    caption = "Grey points do not meet FDR < 0.05. No separate fold-change cutoff is applied."
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )
print(glaucoma_md_plot)

# Save the three-panel overview in the within-study folder.
glaucoma_md_file <- file.path(
  glaucoma_within_dir, "Glaucoma_GSE241782_MD_plots_by_stage.png"
)
ggplot2::ggsave(
  glaucoma_md_file, plot = glaucoma_md_plot,
  width = 13, height = 5.5, units = "in", dpi = 300, bg = "white"
)
print(file.exists(glaucoma_md_file))
