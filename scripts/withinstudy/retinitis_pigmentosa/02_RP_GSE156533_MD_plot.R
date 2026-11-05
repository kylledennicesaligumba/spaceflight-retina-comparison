# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# RP: genome-wide mean-difference (MD) plot for P23H mutant versus WT retina.
# Run in the current within-study session after reviewing the DEG count summary.
# All tested genes are included; significance uses the existing BH FDR < 0.05 rule.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the plotting packages and complete results are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE),
  exists("rp_results"), dir.exists(rp_within_dir)
)

# Make a plotting copy while retaining all tested genes and their annotations.
rp_md_data <- rp_results
stopifnot(
  all(is.finite(rp_md_data$logFC)),
  all(is.finite(rp_md_data$logCPM)),
  all(is.finite(rp_md_data$FDR)),
  all(is.finite(rp_md_data$PValue))
)

# Apply the same FDR threshold and change directions used in the reviewed summary.
rp_md_data$Status <- "Not significant"
rp_md_data$Status[rp_md_data$FDR < 0.05 & rp_md_data$logFC > 0] <- "Higher"
rp_md_data$Status[rp_md_data$FDR < 0.05 & rp_md_data$logFC < 0] <- "Lower"
rp_md_hits <- subset(rp_md_data, Status != "Not significant")

# Choose five labels per direction by lowest P value to keep the figure readable.
# All DEGs are still coloured; this selection only controls printed gene names.
rp_md_hits <- rp_md_hits[order(rp_md_hits$PValue), , drop = FALSE]
rp_md_labels <- rbind(
  head(subset(rp_md_hits, Status == "Higher"), n = 5),
  head(subset(rp_md_hits, Status == "Lower"), n = 5)
)

# Draw the full expression-change overview using the glaucoma MD colour scheme.
rp_md_plot <- ggplot2::ggplot(
  rp_md_data, ggplot2::aes(x = logCPM, y = logFC)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(
    data = subset(rp_md_data, Status == "Not significant"),
    colour = "grey65", alpha = 0.35, size = 0.6
  ) +
  ggplot2::geom_point(
    data = rp_md_hits, ggplot2::aes(colour = Status), size = 1, alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = rp_md_labels,
    ggplot2::aes(label = gene_name, colour = Status),
    size = 3, seed = 1, box.padding = 0.45, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 2,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(Higher = "#D55E00", Lower = "#0072B2"),
    breaks = c("Higher", "Lower")
  ) +
  ggplot2::labs(
    title = "RP: genome-wide differential expression",
    subtitle = "P23H mutant versus WT | sex-adjusted comparison",
    x = "Average log2 CPM (all seven libraries)",
    y = "Log2 fold change (P23H versus WT)", colour = "FDR < 0.05",
    caption = paste(
      "All tested genes shown. Labels: five DEGs with lowest P values per direction.",
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
print(rp_md_plot)

# Save the MD overview in the RP within-study folder.
rp_md_file <- file.path(rp_within_dir, "RP_GSE156533_MD_plot.png")
ggplot2::ggsave(
  rp_md_file, plot = rp_md_plot,
  width = 9, height = 6.5, units = "in", dpi = 300, bg = "white"
)
print(file.exists(rp_md_file))
