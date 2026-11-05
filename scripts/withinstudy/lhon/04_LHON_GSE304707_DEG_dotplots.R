# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# LHON-linked model: individual-library expression of selected DEGs in all four groups.
# Select five higher/five lower combined-genotype DEGs by lowest P value, as for RP.
# Also include every YARS2-only/COI-only comparison DEG, showing shared genes once.

# Confirm that plotting packages and the reviewed input objects are available.
stopifnot(
  requireNamespace("edgeR", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("lhon_within_dge"), exists("lhon_within_metadata"),
  exists("lhon_YARS2_degs"), exists("lhon_COI_degs"),
  exists("lhon_YARS2_COI_degs"), dir.exists(lhon_within_dir)
)

# Select the five lowest-P-value combined-genotype DEGs in each direction.
lhon_dot_candidates <- lhon_YARS2_COI_degs[
  order(lhon_YARS2_COI_degs$PValue), , drop = FALSE
]
lhon_dot_higher <- head(subset(lhon_dot_candidates, logFC > 0), n = 5)
lhon_dot_lower <- head(subset(lhon_dot_candidates, logFC < 0), n = 5)

# Add the single-genotype comparison DEGs and retain each supplied gene name once.
lhon_dot_names <- unique(c(
  lhon_dot_higher$gene_name, lhon_dot_lower$gene_name,
  lhon_YARS2_degs$gene_name, lhon_COI_degs$gene_name
))
lhon_dot_genes <- data.frame(
  gene_name = lhon_dot_names,
  Top5_combined_higher = lhon_dot_names %in% lhon_dot_higher$gene_name,
  Top5_combined_lower = lhon_dot_names %in% lhon_dot_lower$gene_name,
  DEG_YARS2_vs_WT = lhon_dot_names %in% lhon_YARS2_degs$gene_name,
  DEG_COI_vs_WT = lhon_dot_names %in% lhon_COI_degs$gene_name,
  DEG_YARS2_COI_vs_WT = lhon_dot_names %in% lhon_YARS2_COI_degs$gene_name
)
stopifnot(
  nrow(lhon_dot_higher) == 5L, nrow(lhon_dot_lower) == 5L,
  anyDuplicated(lhon_dot_genes$gene_name) == 0L,
  all(lhon_dot_genes$gene_name %in% rownames(lhon_within_dge$counts))
)
print(lhon_dot_genes)

# Calculate observed CPM from the full DGEList using its saved TMM factors.
lhon_dot_cpm <- edgeR::cpm(
  lhon_within_dge, normalized.lib.sizes = TRUE, log = FALSE
)[lhon_dot_genes$gene_name, , drop = FALSE]

# Reshape to one row per gene-library measurement, retaining the supplied gene names.
lhon_dot_data <- as.data.frame(
  as.table(lhon_dot_cpm), responseName = "CPM", stringsAsFactors = FALSE
)
names(lhon_dot_data)[1:2] <- c("gene_name", "sample_id")
lhon_dot_sample_index <- match(lhon_dot_data$sample_id, lhon_within_metadata$sample_id)
stopifnot(!anyNA(lhon_dot_sample_index), all(is.finite(lhon_dot_data$CPM)))
lhon_dot_data$group <- lhon_within_metadata$group[lhon_dot_sample_index]
lhon_dot_data$sex <- lhon_within_metadata$sex[lhon_dot_sample_index]
lhon_dot_data$gene <- factor(lhon_dot_data$gene_name, levels = lhon_dot_genes$gene_name)

# Show one triangle per male library, with small horizontal offsets for visibility.
# Each gene has its own y-axis range; the vertical CPM measurements are unchanged.
lhon_dotplot <- ggplot2::ggplot(
  lhon_dot_data, ggplot2::aes(x = group, y = CPM, colour = group)
) +
  ggplot2::geom_point(
    shape = 17, size = 2.6,
    position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1)
  ) +
  ggplot2::facet_wrap(~ gene, ncol = 3, scales = "free_y") +
  ggplot2::scale_colour_manual(
    values = c(
      WT = "#0072B2", YARS2 = "#D55E00", COI = "#009E73", "YARS2/COI" = "#CC79A7"
    ),
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.03, 0.12))
  ) +
  ggplot2::labs(
    title = "LHON-linked model: individual-library expression of selected DEGs",
    subtitle = paste(
      "Five higher + five lower DEGs in YARS2/COI vs WT, by lowest P value",
      "Plus all DEGs from YARS2 vs WT and COI vs WT; shared genes shown once",
      sep = "\n"
    ),
    x = NULL, y = "TMM-normalized CPM",
    caption = paste(
      "One triangle per library (all male); horizontal offsets keep points visible.",
      "Each gene has its own y-axis range. Points show observed normalized expression.",
      "DEG membership comes from the corresponding genotype-versus-WT tests.",
      sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(face = "italic", size = 11),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0)
  )
print(lhon_dotplot)

# Set paths for the figure, selection key and exact plotted measurements.
lhon_dotplot_files <- file.path(
  lhon_within_dir,
  c(
    "LHON_GSE304707_DEG_expression_dotplots_selected.png",
    "LHON_GSE304707_selected_DEG_dotplot_gene_key.csv",
    "LHON_GSE304707_selected_DEG_dotplot_data.csv"
  )
)

# Save the figure and its data in the LHON output folder.
ggplot2::ggsave(
  lhon_dotplot_files[1], plot = lhon_dotplot,
  width = 13, height = 11, units = "in", dpi = 300, bg = "white"
)
write.csv(lhon_dot_genes, file = lhon_dotplot_files[2], row.names = FALSE)
write.csv(lhon_dot_data, file = lhon_dotplot_files[3], row.names = FALSE)
print(file.exists(lhon_dotplot_files))
