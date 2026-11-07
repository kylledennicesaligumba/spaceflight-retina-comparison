# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Locate the user's downloaded OSD-758 input files.
osd758_dir <- Sys.getenv("OSD758_DATA_DIR", unset = file.path(portfolio_repo_dir, "data", "raw", "osdr758"))

# Keep prequality-check scripts, plots, tables and checkpoints together.
osd758_output_dir <- file.path(portfolio_repo_dir, "results", "osdr758")
dir.create(osd758_output_dir, recursive = TRUE, showWarnings = FALSE)

# Identify the count table and the two downloaded sample-information files.
osd758_count_file <- file.path(osd758_dir, "GLDS-664_rna-seq_RSEM_Unnormalized_Counts.csv")
osd758_sample_file <- file.path(osd758_dir, "GLDS-664_rna-seq_SampleTable.csv")
osd758_metadata_file <- file.path(osd758_dir, "GLDS-664_rna-seq_MHU-8_RTN_metadata.csv")

# Confirm file locations before importing any counts or sample information.
c(
  input_folder = dir.exists(osd758_dir),
  output_folder = dir.exists(osd758_output_dir),
  counts = file.exists(osd758_count_file),
  sample_table = file.exists(osd758_sample_file),
  metadata = file.exists(osd758_metadata_file)
)

# User confirmed all five path checks TRUE on 6 November 2026.

# Preview five count rows while preserving the downloaded column names.
osd758_count_preview <- read.csv(
  osd758_count_file, nrows = 5, check.names = FALSE
)
print(osd758_count_preview)

# Read the small sample-information tables and preview their columns and labels.
osd758_sample_table <- read.csv(osd758_sample_file, check.names = FALSE)
osd758_metadata <- read.csv(osd758_metadata_file, check.names = FALSE)
head(osd758_sample_table)
head(osd758_metadata)

# User-reviewed previews on 6 November 2026 show gene IDs in the count table's
# first column, followed by named sample columns. Its first heading appears blank.
# The metadata preview has named sample and condition columns. FLT labels include
# several gravity conditions: FLT05 is FLT_0G in the displayed metadata.
# These previews do not establish full-table integrity or final sample selection.

# Import the complete downloaded count table without changing column names.
osd758_raw <- read.csv(osd758_count_file, check.names = FALSE)

# Show the total rows and columns, including the gene-identifier column.
dim(osd758_raw)

# Confirm the first column heading and the first five sample-column names.
names(osd758_raw)[1:6]

# Count the recorded conditions, including any missing condition labels.
table(osd758_metadata$condition, useNA = "ifany")

# User-confirmed import on 6 November 2026: 56,840 rows and 36 columns,
# comprising one unnamed identifier column and 35 sample columns.
# Recorded conditions: FLT_0.33G 6; FLT_0.67G 6; FLT_0G 5; FLT_1G 6; HGC 12.
# The intended primary comparison is FLT_0G versus HGC (17 samples before QC).
# Retain all source rows and samples here; no primary subset has been created.

# Preserve the original gene labels separately from the expression values.
osd758_gene_labels <- as.character(osd758_raw[[1]])

# Create the count matrix with original gene labels as row names.
osd758_counts_all <- as.matrix(osd758_raw[, -1, drop = FALSE])
rownames(osd758_counts_all) <- osd758_gene_labels

# Confirm matrix dimensions and check values, identifiers and sample order.
dim(osd758_counts_all)
c(
  numeric_counts = is.numeric(osd758_counts_all),
  finite_counts = all(is.finite(osd758_counts_all)),
  nonnegative_counts = all(osd758_counts_all >= 0),
  gene_ids_present = all(!is.na(osd758_gene_labels) & nzchar(osd758_gene_labels)),
  gene_ids_unique = anyDuplicated(osd758_gene_labels) == 0,
  sample_order = identical(colnames(osd758_counts_all), osd758_metadata$sample)
)

# User confirmed a 56,840 x 35 matrix and all six checks TRUE on 6 November 2026.
# These establish input integrity and ordering, not final sample quality or
# stable-ID annotation for every source row.

# Check that both sample-information files give the same assignments.
c(
  sample_ids_agree = identical(osd758_sample_table[[1]], osd758_metadata$sample),
  conditions_agree = identical(osd758_sample_table$condition, osd758_metadata$condition),
  sample_ids_unique = anyDuplicated(osd758_metadata$sample) == 0
)

# Make a separate metadata table for flight 0G and habitat ground controls.
osd758_metadata_selected <- osd758_metadata[
  osd758_metadata$condition %in% c("FLT_0G", "HGC"),
  , drop = FALSE
]

# Reset display row numbers and review all selected sample labels and groups.
rownames(osd758_metadata_selected) <- NULL
print(osd758_metadata_selected, row.names = FALSE)
table(osd758_metadata_selected$condition)

# User-reviewed selection on 6 November 2026: FLT05, FLT09, FLT13, FLT17,
# FLT21 are FLT_0G; HGC01-HGC12 are HGC, for 17 selected metadata rows.
# The three preceding consistency-check results were not shown in the screenshot;
# include them below with the selected-matrix alignment check for review.

# Select all expression rows for the 17 samples, in the metadata's sample order.
osd758_counts_selected <- osd758_counts_all[
  , osd758_metadata_selected$sample, drop = FALSE
]

# Confirm the selected matrix size, sample alignment and metadata consistency.
dim(osd758_counts_selected)
c(
  sample_order = identical(colnames(osd758_counts_selected), osd758_metadata_selected$sample),
  sample_ids_agree = identical(osd758_sample_table[[1]], osd758_metadata$sample),
  conditions_agree = identical(osd758_sample_table$condition, osd758_metadata$condition),
  sample_ids_unique = anyDuplicated(osd758_metadata$sample) == 0
)

# User confirmed selected dimensions 56,840 x 17 and all four checks TRUE
# on 6 November 2026: selected sample order, agreement of both sample files'
# IDs/conditions, and unique metadata sample IDs.

# Sum supplied counts across all current rows for each selected sample.
osd758_total_counts <- colSums(osd758_counts_selected)

# Display totals in metadata order, then from smallest to largest.
osd758_total_counts
sort(osd758_total_counts)

# User-reviewed totals on 6 November 2026 range from HGC09 = 40,215,185
# to HGC02 = 78,848,510 (about a twofold difference). All 17 samples remain
# for the remaining prequality checks; totals alone do not establish failure.

# Count source rows with positive counts in each selected sample.
osd758_detected_genes <- colSums(osd758_counts_selected > 0)

# Display the detection counts in metadata order and from smallest to largest.
osd758_detected_genes
sort(osd758_detected_genes)

# User-reviewed positive-row counts on 6 November 2026 range from
# HGC04 = 25,487 to HGC08 = 28,005. All 17 samples remain for the next checks.
# No expression filter has yet been applied; all 56,840 source rows remain.

# Calculate log2 CPM for the initial distribution plot, before filtering and TMM.
# The library-size-scaled prior count keeps zero counts finite on the log scale.
osd758_log_cpm <- edgeR::cpm(
  osd758_counts_selected, log = TRUE, prior.count = 2
)

# Define one plot function so the screen and saved figure use the same settings.
osd758_draw_boxplot <- function() {
  osd758_plot_order <- order(factor(
    osd758_metadata_selected$condition, levels = c("HGC", "FLT_0G")
  ))
  osd758_plot_labels <- sub(
    "^MHU-8_RTN_(FLT|HGC)_", "", colnames(osd758_log_cpm)
  )
  osd758_plot_colors <- c(HGC = "#0072B2", FLT_0G = "#D55E00")
  osd758_old_par <- par(mar = c(7, 4.5, 4.5, 1))
  on.exit(par(osd758_old_par))

  # Hide individual outlier dots for readability; all rows contribute to the boxes.
  boxplot(
    osd758_log_cpm[, osd758_plot_order, drop = FALSE],
    names = osd758_plot_labels[osd758_plot_order],
    col = osd758_plot_colors[osd758_metadata_selected$condition[osd758_plot_order]],
    las = 2, outline = FALSE, cex.axis = 0.8,
    ylab = "Log2 counts per million (CPM)"
  )
  title(main = "OSD-758: initial expression distributions", line = 2.5, cex.main = 1.1)
  mtext("Retina | HGC: habitat ground control; FLT: flight 0G",
        side = 3, line = 0.8, cex = 0.65)
}

# Display the boxplot in R.
osd758_draw_boxplot()

# Save the same plot in the OSD-758 prequality output folder.
osd758_boxplot_file <- file.path(osd758_output_dir, "OSD758_logCPM_boxplot.png")
png(osd758_boxplot_file, width = 2400, height = 1800, res = 300)
osd758_draw_boxplot()
dev.off()
file.exists(osd758_boxplot_file)

# User supplied the saved initial boxplot on 6 November 2026; the PNG is present.
# The 17 boxes have broadly similar distributions with modest upper-box and
# whisker differences. Medians overlap the lower box edges near log2 CPM -4.9.
# This agrees with more than half of source rows being zero within each sample,
# based on the reviewed positive-row counts. The prior count makes zeros finite.
# Similar marginal distributions do not establish absent gene-specific flight
# effects or identical sample profiles. Retain all 17 for the remaining checks.

# Estimate one density curve per sample from the same initial log2 CPM matrix.
osd758_density_curves <- lapply(
  seq_len(ncol(osd758_log_cpm)),
  function(i) density(osd758_log_cpm[, i])
)
names(osd758_density_curves) <- colnames(osd758_log_cpm)

# Give each individual sample its own colour, following the RP density plot.
osd758_density_colors <- hcl.colors(ncol(osd758_log_cpm), palette = "Dark 3")

# Define one drawing function for the screen and saved PNG.
osd758_draw_density <- function() {
  # Include the complete estimated range of every sample's curve.
  osd758_xlim <- range(sapply(osd758_density_curves, function(d) range(d$x)))
  osd758_ylim <- c(0, max(sapply(osd758_density_curves, function(d) max(d$y))))

  plot(
    osd758_density_curves[[1]],
    xlim = osd758_xlim, ylim = osd758_ylim,
    col = osd758_density_colors[1], lwd = 2,
    main = "OSD-758: initial expression density",
    xlab = "Log2 counts per million (CPM)", ylab = "Density"
  )
  for (i in 2:length(osd758_density_curves)) {
    lines(osd758_density_curves[[i]], col = osd758_density_colors[i], lwd = 2)
  }

  # Show shortened original sample IDs in the legend.
  legend(
    "topright",
    legend = sub("^MHU-8_RTN_(FLT|HGC)_", "", names(osd758_density_curves)),
    col = osd758_density_colors, lty = 1, lwd = 2,
    bty = "n", cex = 0.65
  )
}

# Display the density plot in R.
osd758_draw_density()

# Save the same plot in the OSD-758 prequality output folder.
osd758_density_file <- file.path(osd758_output_dir, "OSD758_logCPM_density.png")
png(osd758_density_file, width = 2400, height = 1800, res = 300)
osd758_draw_density()
dev.off()
file.exists(osd758_density_file)

# User supplied the saved density plot on 6 November 2026; the PNG is present.
# All 17 curves overlap closely, with a large peak near log2 CPM -5 and a
# smaller broad rise around 4-6. The low-end peak agrees with many zero-count
# source rows and their finite prior-adjusted log-CPM values. Density smoothing
# can extend beyond the observed values. No sample is excluded from this view.
# These distributions do not establish which specific genes agree or change.

# Count source labels that begin with the mouse Ensembl gene-ID prefix.
osd758_ensembl_prefix <- grepl("^ENSMUSG", osd758_gene_labels)
table(Mouse_Ensembl_prefix = osd758_ensembl_prefix)

# Preview other label formats before deciding how their rows should be handled.
head(osd758_gene_labels[!osd758_ensembl_prefix], n = 10)

# User reports 56,748 ENSMUSG-prefixed rows and 92 other rows; the first ten
# other labels are ERCC IDs. Verify all other labels before separating controls.

# Identify ERCC-labelled rows and any labels outside both known formats.
osd758_ercc_rows <- grepl("^ERCC-", osd758_gene_labels)
osd758_other_labels <- osd758_gene_labels[
  !osd758_ensembl_prefix & !osd758_ercc_rows
]

# Report row-type counts and preview any remaining unclassified labels.
c(
  Mouse_Ensembl = sum(osd758_ensembl_prefix),
  ERCC_controls = sum(osd758_ercc_rows),
  Other = length(osd758_other_labels)
)
head(osd758_other_labels, n = 10)

# User confirmed Mouse_Ensembl 56,748, ERCC_controls 92, Other 0 and
# character(0) for other labels on 6 November 2026.

# Make separate mouse-gene and ERCC matrices while preserving the full input.
osd758_counts_mouse <- osd758_counts_selected[
  osd758_ensembl_prefix, , drop = FALSE
]
osd758_counts_ercc <- osd758_counts_selected[
  osd758_ercc_rows, , drop = FALSE
]

# Confirm the size of each matrix after separating the technical controls.
dim(osd758_counts_mouse)
dim(osd758_counts_ercc)

# Check that original mouse-row labels and selected sample order are preserved.
c(
  mouse_row_order = identical(
    rownames(osd758_counts_mouse), osd758_gene_labels[osd758_ensembl_prefix]
  ),
  sample_order = identical(
    colnames(osd758_counts_mouse), osd758_metadata_selected$sample
  )
)

# User confirmed mouse dimensions 56,748 x 17, ERCC dimensions 92 x 17,
# and both row/sample-order checks TRUE on 6 November 2026.

# Define the sample groups with habitat ground controls as the reference level.
osd758_qc_group <- factor(
  osd758_metadata_selected$condition, levels = c("HGC", "FLT_0G")
)
table(osd758_qc_group)

# Store mouse-gene counts and groups in edgeR's input object.
osd758_dge <- edgeR::DGEList(
  counts = osd758_counts_mouse,
  group = osd758_qc_group
)

# Confirm the object dimensions and preservation of sample and gene-row order.
dim(osd758_dge)
c(
  sample_order = identical(colnames(osd758_dge), osd758_metadata_selected$sample),
  gene_order = identical(rownames(osd758_dge), rownames(osd758_counts_mouse))
)

# User confirmed HGC 12 / FLT_0G 5, dimensions 56,748 x 17, and both
# alignment checks TRUE on 6 November 2026.

# Apply the same expression-filter settings used for the other count datasets.
osd758_keep_genes <- edgeR::filterByExpr(
  osd758_dge,
  group = osd758_qc_group,
  min.count = 10,
  min.total.count = 15
)

# Report retained (TRUE) and removed (FALSE) mouse-gene rows.
table(Retain_for_QC = osd758_keep_genes)

# Create a filtered copy and recalculate library totals over its retained rows.
osd758_dge_qc <- osd758_dge[
  osd758_keep_genes, , keep.lib.sizes = FALSE
]

# Confirm that all 17 samples and the correct retained gene-row order remain.
dim(osd758_dge_qc)
c(
  sample_order = identical(colnames(osd758_dge_qc), osd758_metadata_selected$sample),
  gene_order = identical(
    rownames(osd758_dge_qc), rownames(osd758_dge)[osd758_keep_genes]
  )
)

# User confirmed 35,413 removed / 21,335 retained mouse-gene rows,
# filtered dimensions 21,335 x 17 and both alignment checks TRUE.
# The 92 ERCC controls were separated before this expression-filter input.

# Estimate TMM library-composition scaling factors for the filtered mouse counts.
osd758_dge_qc <- edgeR::normLibSizes(osd758_dge_qc, method = "TMM")

# Display each sample's group, retained-row library total and TMM factor.
osd758_dge_qc$samples[, c("group", "lib.size", "norm.factors")]

# Confirm that every estimated normalization factor is finite and positive.
c(
  finite_factors = all(is.finite(osd758_dge_qc$samples$norm.factors)),
  positive_factors = all(osd758_dge_qc$samples$norm.factors > 0)
)

# User confirmed all 17 TMM factors and both finite/positive checks TRUE.
# Displayed factors range from FLT13 = 0.8837160 to FLT17 = 1.1168971.
# These are modest library-composition adjustments, not biological effect sizes.

# Calculate filtered, TMM-adjusted log2 CPM for sample-relationship plots.
osd758_log_cpm_qc <- edgeR::cpm(
  osd758_dge_qc,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

# Confirm the new matrix size, finite values and selected sample order.
dim(osd758_log_cpm_qc)
c(
  finite_logCPM = all(is.finite(osd758_log_cpm_qc)),
  sample_order = identical(colnames(osd758_log_cpm_qc), osd758_metadata_selected$sample)
)

# User confirmed the finite-logCPM and sample-order checks TRUE on 6 November 2026.
# The log-CPM transformation preserves the filtered 21,335 x 17 matrix shape.

# Load the familiar point/label style with an outside condition legend.
source(file.path(portfolio_repo_dir, "scripts", "prechecks", "osdr758", "02_OSD758_projection_style.R"))

# Calculate pairwise top-500 MDS and common top-500 PCA from the same matrix.
osd758_mds <- limma::plotMDS(
  osd758_log_cpm_qc, top = 500, gene.selection = "pairwise", plot = FALSE
)
osd758_pca <- limma::plotMDS(
  osd758_log_cpm_qc, top = 500, gene.selection = "common", plot = FALSE
)

# Format the calculated coordinates, sample labels and axis percentages.
osd758_mds_plot <- osd758_projection_plot(
  osd758_mds, osd758_metadata_selected,
  title = "OSD-758: MDS after TMM",
  axis_names = c("Leading logFC dim 1", "Leading logFC dim 2")
)
osd758_pca_plot <- osd758_projection_plot(
  osd758_pca, osd758_metadata_selected,
  title = "OSD-758: PCA after TMM",
  axis_names = c("PC1", "PC2")
)

# Display both plots; each point represents one selected sample.
print(osd758_mds_plot)
print(osd758_pca_plot)

# Set separate filenames in the OSD-758 prequality output folder.
osd758_mds_file <- file.path(osd758_output_dir, "OSD758_TMM_MDS.png")
osd758_pca_file <- file.path(osd758_output_dir, "OSD758_TMM_PCA_top500.png")

# Save each plot and close its PNG device.
png(osd758_mds_file, width = 2400, height = 1800, res = 300)
print(osd758_mds_plot)
dev.off()
png(osd758_pca_file, width = 2400, height = 1800, res = 300)
print(osd758_pca_plot)
dev.off()

# Confirm both files exist before reviewing the figures together.
c(MDS = file.exists(osd758_mds_file), PCA = file.exists(osd758_pca_file))

# User supplied both saved plots on 6 November 2026; both PNGs are present.
# MDS axes show 20%/15%, and common-top-500 PCA axes show 25%/19%.
# FLT09/FLT17 lie furthest right in both views; the other flight points overlap
# the ground-sample range. Flight samples do not form one compact cluster.
# HGC05 is distinctly high on the second axis in both views. Retain it with a
# recorded review flag; these figures alone do not establish a failed sample.
# All 17 samples and 21,335 filtered mouse rows remain at this stage.
# Source study-design/provenance checks remain before within-study modelling.

# Save the reviewed stage, its sample flag and the current software versions.
source(file.path(portfolio_repo_dir, "scripts", "prechecks", "osdr758", "03_OSD758_save_current_stage.R"))
