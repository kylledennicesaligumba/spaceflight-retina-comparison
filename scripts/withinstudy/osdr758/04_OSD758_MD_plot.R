# OSD-758: genome-wide mean-difference plot for flight 0G versus HGC retina.
# Run after reviewing the DEG summary; use the familiar RP/glaucoma plot style.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Confirm that the plotting packages and complete results are available.
stopifnot(
  requireNamespace("ggplot2", quietly = TRUE),
  requireNamespace("ggrepel", quietly = TRUE),
  exists("osd758_results"), dir.exists(osd758_within_dir)
)

# Make a plotting copy containing every tested gene and its original identifier.
osd758_md_data <- osd758_results
stopifnot(
  all(is.finite(osd758_md_data$logFC)),
  all(is.finite(osd758_md_data$logCPM)),
  all(is.finite(osd758_md_data$FDR)),
  all(is.finite(osd758_md_data$PValue))
)

# Apply the reviewed FDR threshold and flight-versus-HGC change directions.
osd758_md_data$Status <- "Not significant"
osd758_md_data$Status[
  osd758_md_data$FDR < 0.05 & osd758_md_data$logFC > 0
] <- "Higher"
osd758_md_data$Status[
  osd758_md_data$FDR < 0.05 & osd758_md_data$logFC < 0
] <- "Lower"
osd758_md_hits <- subset(osd758_md_data, Status != "Not significant")

# Label five lowest-P DEGs per direction while colouring every significant gene.
# Use the original Ensembl IDs because gene symbols have not been mapped yet.
osd758_md_hits <- osd758_md_hits[order(osd758_md_hits$PValue), , drop = FALSE]
osd758_md_labels <- rbind(
  head(subset(osd758_md_hits, Status == "Higher"), n = 5),
  head(subset(osd758_md_hits, Status == "Lower"), n = 5)
)

# Plot average expression against the estimated flight-versus-HGC log2 change.
osd758_md_plot <- ggplot2::ggplot(
  osd758_md_data, ggplot2::aes(x = logCPM, y = logFC)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_point(
    data = subset(osd758_md_data, Status == "Not significant"),
    colour = "grey65", alpha = 0.35, size = 0.6
  ) +
  ggplot2::geom_point(
    data = osd758_md_hits, ggplot2::aes(colour = Status), size = 1, alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = osd758_md_labels,
    ggplot2::aes(label = gene_id, colour = Status),
    size = 2.6, seed = 1, box.padding = 0.45, point.padding = 0.3,
    min.segment.length = 0, max.overlaps = Inf, max.time = 2,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(Higher = "#D55E00", Lower = "#0072B2"),
    breaks = c("Higher", "Lower")
  ) +
  ggplot2::labs(
    title = "OSD-758: genome-wide differential expression",
    subtitle = "Flight 0G versus habitat ground control | mouse retina",
    x = "Average log2 CPM (all 17 libraries)",
    y = "Log2 fold change (flight 0G versus HGC)", colour = "FDR < 0.05",
    caption = paste(
      "All 21,335 tested genes shown. Labels: five lowest-P DEGs per direction (Ensembl IDs).",
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
print(osd758_md_plot)

# Save the MD figure in the OSD-758 within-study folder and check its presence.
osd758_md_file <- file.path(osd758_within_dir, "OSD758_MD_plot.png")
ggplot2::ggsave(
  osd758_md_file, plot = osd758_md_plot,
  width = 10, height = 7, units = "in", dpi = 300, bg = "white"
)
print(file.exists(osd758_md_file))
