# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# DR: mean-difference (MD) plot for Diabetic versus Normal retina.
# Run in the current within-study session after reviewing the DEG summary.
# Use the existing limma results; this script does not refit or retest the model.
# MD axes: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

# Confirm that the plotting packages, complete results and output folder are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE),
  exists("dr_results"), dir.exists(dr_within_dir)
)

# Make a plotting copy, keeping every tested gene and its supplied annotation.
dr_md_data <- dr_results
stopifnot(
  all(is.finite(dr_md_data$AveExpr)),
  all(is.finite(dr_md_data$logFC)),
  all(is.finite(dr_md_data$adj.P.Val))
)

# Apply the same FDR cutoff and directions used in the reviewed DEG summary.
dr_md_data$Status <- "Not significant"
dr_md_data$Status[
  dr_md_data$adj.P.Val < 0.05 & dr_md_data$logFC > 0
] <- "Higher"
dr_md_data$Status[
  dr_md_data$adj.P.Val < 0.05 & dr_md_data$logFC < 0
] <- "Lower"
dr_md_hits <- subset(dr_md_data, Status != "Not significant")

# Plot average processed expression against the signed Diabetic-minus-Normal coefficient.
# Both axes use the supplied microarray expression scale; no CPM or extra log is computed.
dr_md_plot <- ggplot2::ggplot(
  dr_md_data, ggplot2::aes(x = AveExpr, y = logFC)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(
    data = subset(dr_md_data, Status == "Not significant"),
    colour = "grey65", alpha = 0.35, size = 0.6
  ) +
  ggplot2::geom_point(
    data = dr_md_hits, ggplot2::aes(colour = Status), size = 2.2
  ) +
  ggrepel::geom_text_repel(
    data = dr_md_hits,
    ggplot2::aes(label = Symbol, colour = Status),
    size = 4, seed = 1, box.padding = 0.6, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 2,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(Higher = "#D55E00", Lower = "#0072B2"),
    breaks = c("Higher", "Lower")
  ) +
  ggplot2::labs(
    title = "DR: genome-wide differential expression",
    subtitle = "Diabetic versus Normal | 6 arrays per group",
    x = "Average processed expression (12 selected arrays)",
    y = "Expression difference (Diabetic - Normal)",
    colour = "FDR < 0.05",
    caption = paste(
      "All tested genes shown. All DEGs labelled.",
      "Grey points do not meet FDR < 0.05. No separate fold-change cutoff is applied.",
      sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )

# Display the plot, highlighting and labelling all three reviewed DEGs.
print(dr_md_plot)

# Save the MD plot in the DR within-study output folder.
dr_md_file <- file.path(dr_within_dir, "DR_GSE87433_MD_plot.png")
ggplot2::ggsave(
  dr_md_file, plot = dr_md_plot,
  width = 9, height = 6.5, units = "in", dpi = 300, bg = "white"
)

# Confirm that R can find the saved figure.
print(file.exists(dr_md_file))
