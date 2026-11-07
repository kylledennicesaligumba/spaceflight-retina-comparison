# OSD-758: individual-library expression for five higher and five lower DEGs.
# Use the familiar RP selection and TMM-normalized CPM display.

# Confirm that the plotting package and reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("osd758_degs"), exists("osd758_within_dge"),
  exists("osd758_within_metadata"), dir.exists(osd758_within_dir)
)

# Select five lowest-P DEGs per direction, matching the MD label selection.
osd758_dot_candidates <- osd758_degs[order(osd758_degs$PValue), , drop = FALSE]
osd758_dot_genes <- rbind(
  transform(head(subset(osd758_dot_candidates, logFC > 0), n = 5), Direction = "Higher"),
  transform(head(subset(osd758_dot_candidates, logFC < 0), n = 5), Direction = "Lower")
)
stopifnot(
  nrow(osd758_dot_genes) == 10L,
  anyDuplicated(osd758_dot_genes$gene_id) == 0L,
  all(osd758_dot_genes$gene_id %in% rownames(osd758_within_dge$counts))
)

# Calculate CPM from the full DGEList with its saved TMM factors, then select rows.
osd758_dot_cpm <- edgeR::cpm(
  osd758_within_dge, normalized.lib.sizes = TRUE, log = FALSE
)[osd758_dot_genes$gene_id, , drop = FALSE]

# Keep one row per gene-library measurement and attach its group and direction.
osd758_dot_data <- as.data.frame(
  as.table(osd758_dot_cpm), responseName = "CPM", stringsAsFactors = FALSE
)
names(osd758_dot_data)[1:2] <- c("gene_id", "sample")
osd758_dot_sample_index <- match(osd758_dot_data$sample, osd758_within_metadata$sample)
osd758_dot_gene_index <- match(osd758_dot_data$gene_id, osd758_dot_genes$gene_id)
stopifnot(!anyNA(osd758_dot_sample_index), !anyNA(osd758_dot_gene_index))
osd758_dot_data$group <- osd758_within_metadata$group[osd758_dot_sample_index]
osd758_dot_data$Direction <- osd758_dot_genes$Direction[osd758_dot_gene_index]
osd758_dot_data$gene <- factor(osd758_dot_data$gene_id, levels = osd758_dot_genes$gene_id)

# Show each library with horizontal offsets only; each panel has its own CPM range.
osd758_dotplot <- ggplot2::ggplot(
  osd758_dot_data, ggplot2::aes(x = group, y = CPM, colour = group)
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
    size = 2.5, alpha = 0.85
  ) +
  ggplot2::facet_wrap(~ gene, ncol = 2, dir = "v", scales = "free_y") +
  ggplot2::scale_colour_manual(
    values = c(HGC = "#0072B2", FLT_0G = "#D55E00"), guide = "none"
  ) +
  ggplot2::scale_x_discrete(labels = c(HGC = "HGC", FLT_0G = "Flight 0G")) +
  ggplot2::scale_y_continuous(
    limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.03, 0.12))
  ) +
  ggplot2::labs(
    title = "OSD-758: individual-library expression of selected DEGs",
    subtitle = paste(
      "Left: five higher DEGs | Right: five lower DEGs",
      "Lowest P values within direction among genes with FDR < 0.05", sep = "\n"
    ),
    x = NULL, y = "TMM-normalized CPM",
    caption = paste(
      "One point per library: 12 HGC and 5 flight 0G. Horizontal offsets keep points visible.",
      "Each panel has its own y-axis range. Labels retain original Ensembl gene IDs.", sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )
print(osd758_dotplot)

# Save the figure, exact point measurements and selected-gene key.
osd758_dotplot_files <- file.path(
  osd758_within_dir,
  c(
    "OSD758_DEG_expression_dotplots_top10.png",
    "OSD758_DEG_dotplot_top10_data.csv",
    "OSD758_top10_DEG_dotplot_selected_genes.csv"
  )
)
ggplot2::ggsave(
  osd758_dotplot_files[1], plot = osd758_dotplot,
  width = 11, height = 12, units = "in", dpi = 300, bg = "white"
)
write.csv(osd758_dot_data, osd758_dotplot_files[2], row.names = FALSE)
write.csv(osd758_dot_genes, osd758_dotplot_files[3], row.names = FALSE)
print(file.exists(osd758_dotplot_files))
