# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# LHON-linked model: three genome-wide mean-difference (MD) plots.
# Run in the current within-study R session after reviewing the genotype summary.
# All tested genes are included; significance uses the existing BH FDR < 0.05 rule.
# Method: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the plotting packages and complete result tables are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE),
  exists("lhon_YARS2_results"), exists("lhon_COI_results"),
  exists("lhon_YARS2_COI_results"), dir.exists(lhon_within_dir)
)

# Stack the full result tables for plotting, keeping each comparison identified.
lhon_md_data <- rbind(
  transform(lhon_YARS2_results, Comparison = "YARS2"),
  transform(lhon_COI_results, Comparison = "COI"),
  transform(lhon_YARS2_COI_results, Comparison = "YARS2/COI")
)
lhon_md_data$Comparison <- factor(
  lhon_md_data$Comparison,
  levels = c("YARS2", "COI", "YARS2/COI"),
  labels = c("YARS2 vs WT", "COI vs WT", "YARS2/COI vs WT")
)

# Stop for review if any statistic needed for plotting or label selection is nonfinite.
stopifnot(
  all(is.finite(lhon_md_data$logFC)),
  all(is.finite(lhon_md_data$logCPM)),
  all(is.finite(lhon_md_data$FDR)),
  all(is.finite(lhon_md_data$PValue))
)

# Apply the same FDR threshold and directions used in the reviewed summaries.
lhon_md_data$Status <- "Not significant"
lhon_md_data$Status[
  lhon_md_data$FDR < 0.05 & lhon_md_data$logFC > 0
] <- "Higher"
lhon_md_data$Status[
  lhon_md_data$FDR < 0.05 & lhon_md_data$logFC < 0
] <- "Lower"
lhon_md_hits <- subset(lhon_md_data, Status != "Not significant")

# Choose up to five gene labels per direction within each comparison, by lowest P value.
# All DEGs remain coloured; only the printed names are limited for readability.
lhon_md_labels <- do.call(
  rbind,
  lapply(
    split(
      lhon_md_hits,
      list(lhon_md_hits$Comparison, lhon_md_hits$Status),
      drop = TRUE
    ),
    function(x) head(x[order(x$PValue), , drop = FALSE], n = 5)
  )
)

# Draw all tested genes, then highlight DEGs and add the selected gene names.
# Keep fixed axis ranges so all three panels use the same expression/change scales.
lhon_md_plot <- ggplot2::ggplot(
  lhon_md_data, ggplot2::aes(x = logCPM, y = logFC)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(
    data = subset(lhon_md_data, Status == "Not significant"),
    colour = "grey65", alpha = 0.35, size = 0.6
  ) +
  ggplot2::geom_point(
    data = lhon_md_hits, ggplot2::aes(colour = Status), size = 1, alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = lhon_md_labels,
    ggplot2::aes(label = gene_name, colour = Status),
    size = 3, seed = 1, box.padding = 0.45, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 2,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Comparison, nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c(Higher = "#D55E00", Lower = "#0072B2"),
    breaks = c("Higher", "Lower")
  ) +
  ggplot2::labs(
    title = "LHON-linked model: genome-wide differential expression",
    subtitle = "Genotype-versus-WT comparisons | all tested genes shown",
    x = "Average log2 CPM (all 12 libraries)",
    y = "Log2 fold change versus WT", colour = "FDR < 0.05",
    caption = paste(
      "Labels: up to five DEGs with lowest P values per direction in each comparison.",
      "Grey points do not meet FDR < 0.05. No separate fold-change cutoff is applied.",
      sep = "\n"
    )
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
print(lhon_md_plot)

# Save the three-panel MD overview in the LHON within-study output folder.
lhon_md_file <- file.path(
  lhon_within_dir, "LHON_GSE304707_MD_plots_by_genotype.png"
)
ggplot2::ggsave(
  lhon_md_file, plot = lhon_md_plot,
  width = 13, height = 5.5, units = "in", dpi = 300, bg = "white"
)
print(file.exists(lhon_md_file))
