# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Plot the reviewed gene-level Spearman results without recalculating the scores.
# Show every contrast on the full correlation scale with its numerical estimate.

# Restore the saved score table and confirm that the plotting package is available.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
main_gene_plot_input_file <- file.path(
  main_analysis_dir, "Main_analysis_gene_Spearman_checkpoint_2026-11-06.rds"
)
stopifnot(file.exists(main_gene_plot_input_file))
main_gene_plot_context <- readRDS(main_gene_plot_input_file)
main_gene_similarity_plot_data <- main_gene_plot_context$results
stopifnot(
  nrow(main_gene_similarity_plot_data) == sum(main_gene_plot_context$contrast_inputs$Role == "Terrestrial"),
  all(is.finite(main_gene_similarity_plot_data$Spearman_rho)),
  all(abs(main_gene_similarity_plot_data$Spearman_rho) <= 1),
  all(main_gene_similarity_plot_data$Common_genes == length(main_gene_plot_context$common_ids))
)

# Give each model/stage/genotype/age a readable label while preserving its exact contrast key.
main_gene_plot_labels <- c(
  Glaucoma__3D_vs_Naive = "Glaucoma: 3 days vs naive",
  Glaucoma__2W_vs_Naive = "Glaucoma: 2 weeks vs naive",
  Glaucoma__6W_vs_Naive = "Glaucoma: 6 weeks vs naive",
  Autoimmune_uveitis__EAU_vs_CON = "Autoimmune uveitis: EAU vs control",
  DR__Diabetic_vs_Normal = "Diabetic retinal injury: diabetic vs normal",
  LHON__YARS2_vs_WT = "LHON-linked model: YARS2 vs WT",
  LHON__COI_vs_WT = "LHON-linked model: COI vs WT",
  LHON__YARS2_COI_vs_WT = "LHON-linked model: YARS2 + COI vs WT",
  AMD__Young_CfhKO_vs_WT = "AMD-related model: young CfhKO vs WT",
  AMD__Aged_CfhKO_vs_WT = "AMD-related model: aged CfhKO vs WT",
  RP__P23H_vs_WT = "RP model: P23H vs WT",
  Light_damage__LD_vs_CTRL = "Light damage: LD vs control"
)
main_gene_similarity_plot_data$Contrast_key <- paste(
  main_gene_similarity_plot_data$Study, main_gene_similarity_plot_data$Comparison, sep = "__"
)
stopifnot(
  all(main_gene_similarity_plot_data$Contrast_key %in% names(main_gene_plot_labels)),
  anyDuplicated(main_gene_similarity_plot_data$Contrast_key) == 0L
)

# Place the highest observed rho at the top; this is descriptive ordering.
main_gene_similarity_plot_data <- main_gene_similarity_plot_data[
  order(-main_gene_similarity_plot_data$Spearman_rho,
    main_gene_similarity_plot_data$Study, main_gene_similarity_plot_data$Comparison), , drop = FALSE
]
main_gene_similarity_plot_data$Plot_label <- unname(
  main_gene_plot_labels[main_gene_similarity_plot_data$Contrast_key]
)
main_gene_similarity_plot_data$Plot_label <- factor(
  main_gene_similarity_plot_data$Plot_label,
  levels = rev(main_gene_similarity_plot_data$Plot_label)
)
main_gene_similarity_plot_data$Rho_label <- sprintf(
  "%+.3f", main_gene_similarity_plot_data$Spearman_rho
)
main_gene_similarity_plot_data$Label_hjust <- ifelse(
  main_gene_similarity_plot_data$Spearman_rho >= 0, -0.35, 1.35
)
main_gene_plot_common_n <- unique(main_gene_similarity_plot_data$Common_genes)
stopifnot(length(main_gene_plot_common_n) == 1L)

# Draw point estimates from zero using a fixed -1 to +1 scale and a uniform colour.
main_gene_similarity_plot <- ggplot2::ggplot(
  main_gene_similarity_plot_data,
  ggplot2::aes(x = Spearman_rho, y = Plot_label)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Spearman_rho, yend = Plot_label),
    linewidth = 0.65, colour = "grey65"
  ) +
  ggplot2::geom_point(size = 2.8, colour = "#0072B2") +
  ggplot2::geom_text(
    ggplot2::aes(label = Rho_label, hjust = Label_hjust),
    size = 3.3, colour = "grey15"
  ) +
  ggplot2::scale_x_continuous(breaks = seq(-1, 1, by = 0.25)) +
  ggplot2::coord_cartesian(xlim = c(-1, 1), clip = "off") +
  ggplot2::labs(
    title = "Gene-score agreement with the OSD-758 flight response",
    subtitle = paste0(
      "Reference: flight 0G versus HGC | ",
      format(main_gene_plot_common_n, big.mark = ",", scientific = FALSE),
      " common mouse genes in every contrast"
    ),
    x = "Spearman rho (signed gene-score ordering)", y = NULL,
    caption = paste(
      "Full correlation scale (-1 to +1). All common genes included regardless of FDR.",
      "Descriptive point estimates; uncertainty intervals have not been estimated.",
      "Stages, genotypes and ages remain separate contrasts.", sep = "\n"
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    axis.ticks.y = ggplot2::element_blank(),
    plot.title.position = "plot", plot.caption.position = "plot",
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0, size = 9),
    plot.margin = ggplot2::margin(12, 24, 12, 12)
  )
print(main_gene_similarity_plot)

# Save the figure and its exact plotted data beside the original score table.
main_gene_similarity_plot_files <- file.path(
  main_analysis_dir,
  c("Main_analysis_OSD758_gene_Spearman_common.png", "Main_analysis_OSD758_gene_Spearman_plot_data.csv")
)
ggplot2::ggsave(
  main_gene_similarity_plot_files[1], plot = main_gene_similarity_plot,
  width = 11.5, height = 8, units = "in", dpi = 300, bg = "white"
)
utils::write.csv(main_gene_similarity_plot_data, main_gene_similarity_plot_files[2], row.names = FALSE)
print(file.exists(main_gene_similarity_plot_files))
