# Portable repository copy; run from the repository root.
# Original working files are preserved. GitHub sync did not execute this code.
# LHON-linked mitochondrial optic neuropathy model: GSE304707
# Manual precheck; organized on 30 October 2026 to follow RP/glaucoma.
# Run one explained block at a time and review its output before continuing.
# The current session already has metadata and matching objects. Do not
# source this whole script to resume; use 06_LHON_GSE304707_reload.R.
# See docs/lhon-precheck.md for the completed-state record.
# Earlier chronological teaching code remains on original local storage.

# STEP 0: Set the input and output locations.
# Store paths so input files and output destinations stay explicit.
lhon_dir <- Sys.getenv("LHON_DATA_DIR", unset = file.path("data", "raw", "LHON_GSE304707"))
lhon_output_dir <- file.path("results", "lhon")
dir.create(lhon_output_dir, recursive = TRUE, showWarnings = FALSE)
lhon_count_file <- file.path(lhon_dir, "GSE304707_gene_count_matrix.csv")
lhon_metadata_archive <- file.path(lhon_dir, "GSE304707_family")

# Check that R can find the output folder and both downloaded inputs.
c(
  output_folder = dir.exists(lhon_output_dir),
  counts = file.exists(lhon_count_file),
  metadata = file.exists(lhon_metadata_archive)
)

# STEP 1: Import the table and inspect its structure.
# Keep original column labels and text fields as recorded in the CSV.
lhon_raw <- read.csv(
  lhon_count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the table size, column labels and first five rows.
dim(lhon_raw)
colnames(lhon_raw)
head(lhon_raw, 5)
# User-confirmed: 49,714 rows x 13 columns, comprising gene_name and 12 samples.

# STEP 2: Check missing values, gene labels and numeric count columns.
# Completed: user reports no missing or duplicated data and all sample
# columns numeric. No entries have been removed or replaced.
# Check for missing entries anywhere in the gene/count table.
anyNA(lhon_raw)

# Check for missing or blank gene labels.
sum(is.na(lhon_raw$gene_name) | trimws(lhon_raw$gene_name) == "")

# Return 0 if gene labels are distinct, or the first duplicate's position.
# If duplicates are present, review them before deciding how to handle labels.
anyDuplicated(lhon_raw$gene_name)

# Check that each sample-count column is stored as numeric data.
sapply(lhon_raw[, -1, drop = FALSE], is.numeric)

# STEP 3: Build the numeric count matrix and check the count values.
# Completed: user screenshot confirms 49,714 x 12, finite TRUE,
# range 0 to 1,058,428, whole-valued TRUE, and metadata/matrix order TRUE.
# Exclude gene_name from the numeric entries; retain it as row labels.
lhon_counts_all <- as.matrix(lhon_raw[, -1, drop = FALSE])
rownames(lhon_counts_all) <- lhon_raw$gene_name

# Check the matrix size. Expected from the imported table: 49,714 x 12.
dim(lhon_counts_all)

# Check for finite values and inspect the minimum and maximum counts.
all(is.finite(lhon_counts_all))
range(lhon_counts_all)

# Check whole-valuedness without rounding or changing stored counts.
all(lhon_counts_all == round(lhon_counts_all))
# Review actual outputs before total counts or further expression checks.

# STEP 4A: List and extract the metadata archive.
# These archive/extraction steps were completed together on 29 October.
# List the archive contents without extracting them.
untar(lhon_metadata_archive, list = TRUE)
# User-confirmed member: GSE304707_family.xml.

# Keep extracted metadata in its own output subfolder.
lhon_metadata_dir <- file.path(lhon_output_dir, "metadata")
dir.create(lhon_metadata_dir, recursive = TRUE, showWarnings = FALSE)
untar(
  lhon_metadata_archive,
  files = "GSE304707_family.xml",
  exdir = lhon_metadata_dir
)
lhon_metadata_file <- file.path(lhon_metadata_dir, "GSE304707_family.xml")
file.exists(lhon_metadata_file)

# STEP 4B: Load the XML and build the sample-identity data frame.
# Read and preview the metadata, using the same pattern as RP.
lhon_metadata_xml <- xml2::read_xml(lhon_metadata_file)
lhon_metadata_xml

# Use the XML namespace shown in the user's overview for searches.
lhon_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")
lhon_sample_nodes <- xml2::xml_find_all(
  lhon_metadata_xml,
  "/geo:MINiML/geo:Sample",
  ns = lhon_ns
)
length(lhon_sample_nodes)

# Create one row per sample, preserving the XML order.
# iid holds its GSM accession; Title supplies its descriptive name.
lhon_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(lhon_sample_nodes, "iid"),
  title = xml2::xml_text(
    xml2::xml_find_first(lhon_sample_nodes, "geo:Title", ns = lhon_ns)
  ),
  stringsAsFactors = FALSE
)
lhon_metadata_all

# STEP 4C: Inspect and add the recorded biological characteristics.
# Inventory the available field names before selecting their values.
lhon_characteristic_nodes <- xml2::xml_find_all(
  lhon_sample_nodes,
  ".//geo:Characteristics",
  ns = lhon_ns
)
lhon_characteristic_tags <- xml2::xml_attr(lhon_characteristic_nodes, "tag")
table(lhon_characteristic_tags, useNA = "ifany")
# User-confirmed labels: age, genotype, Sex, strain and tissue; 12 occurrences each.

# Define the same kind of reusable characteristic reader used for RP.
# Missing fields return NA; trim = TRUE removes surrounding whitespace.
lhon_get_characteristic <- function(tag) {
  lhon_path <- paste0(".//geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(lhon_sample_nodes, lhon_path, ns = lhon_ns),
    trim = TRUE
  )
}

# Add the five characteristic columns, preserving the metadata row order.
lhon_metadata_all$genotype <- lhon_get_characteristic("genotype")
lhon_metadata_all$sex <- lhon_get_characteristic("Sex")
lhon_metadata_all$age <- lhon_get_characteristic("age")
lhon_metadata_all$tissue <- lhon_get_characteristic("tissue")
lhon_metadata_all$strain <- lhon_get_characteristic("strain")
lhon_metadata_all

# STEP 4D: Review the extraction protocol for additional timing information.
# Preserve one description per sample and display distinct descriptions once.
lhon_extraction_protocol <- xml2::xml_text(
  xml2::xml_find_first(
    lhon_sample_nodes, ".//geo:Extract-Protocol", ns = lhon_ns
  ),
  trim = TRUE
)
names(lhon_extraction_protocol) <- xml2::xml_attr(lhon_sample_nodes, "iid")
cat(unique(lhon_extraction_protocol), sep = "\n\n")
cat("\n")
# User-confirmed limitation: collection/RNA-extraction timing is not reported
# in the inspected Characteristics and Extract-Protocol fields. Age is 12 months.
# See docs/lhon-metadata-notes.md for the scope of this declaration.

# STEP 4E: Obtain the count-column identifier recorded for each sample.
# Read the Description field and inspect it beside the sample title.
lhon_metadata_all$description <- xml2::xml_text(
  xml2::xml_find_first(lhon_sample_nodes, ".//geo:Description", ns = lhon_ns),
  trim = TRUE
)
lhon_metadata_all[, c("title", "description")]

# The user's output records the count-column name after the final colon-space.
# Extract that label as sample_id; retain the description and GSM accession.
lhon_metadata_all$sample_id <- trimws(
  sub(".*: ", "", lhon_metadata_all$description)
)

# STEP 5: Match metadata rows to count-matrix columns, as in RP.
# The metadata match was first built using CSV headers. The user has now
# confirmed its exact ID/order agreement against the numeric matrix (TRUE).
lhon_sample_ids <- colnames(lhon_counts_all)

# Require unique IDs and the same sample set in both sources.
stopifnot(
  !anyNA(lhon_metadata_all$sample_id),
  anyDuplicated(lhon_metadata_all$sample_id) == 0L,
  anyDuplicated(lhon_sample_ids) == 0L,
  setequal(lhon_metadata_all$sample_id, lhon_sample_ids)
)

# Arrange a metadata copy in the count-column order, retaining all fields.
lhon_metadata_matched <- lhon_metadata_all[
  match(lhon_sample_ids, lhon_metadata_all$sample_id),
  , drop = FALSE
]
rownames(lhon_metadata_matched) <- NULL

# Verify the number of samples and exact ID/order agreement.
dim(lhon_metadata_matched)
identical(lhon_metadata_matched$sample_id, colnames(lhon_counts_all))

# STEP 6: Review the matched biological characteristics and sample groups.
# Completed in the current session using the existing matched metadata object.
# Display sample identity and all five characteristics; all other fields
# remain stored in lhon_metadata_matched.
lhon_metadata_matched[, c(
  "geo_accession", "sample_id", "genotype", "sex", "age", "tissue", "strain"
)]

# Count the libraries in each genotype/sex combination.
# These counts describe measured libraries, not verified independent animals.
table(
  lhon_metadata_matched$genotype,
  lhon_metadata_matched$sex,
  useNA = "ifany"
)
# User-confirmed: all 12 libraries are male; WT, YARS2, COI and YARS2/COI
# each have three libraries. No female or missing-sex category is displayed.

# STEP 7A: Calculate total gene-matrix counts per sample.
# Sum every gene row within each sample column.
lhon_total_counts <- colSums(lhon_counts_all)

# Display the named totals in the original sample order.
lhon_total_counts

# Display them from smallest to largest for comparison.
# This sorts the displayed vector, not the count matrix or metadata.
sort(lhon_total_counts)

# User-confirmed: totals range from WT-2 = 55,648,273 to Y-3 = 65,678,112.
# Highest/lowest is approximately 1.18 (18% higher); totals alone do not
# establish sample quality or justify excluding a library.
# These totals summarize the deposited gene-count matrix. The earlier
# maximum of 1,058,428 was one gene/sample entry, not a sample's total.

# STEP 7B: Count detected gene rows per sample, as in RP/glaucoma.
# For this initial summary, a gene row is detected when its raw count is > 0.
# The comparison creates TRUE/FALSE entries; colSums() counts TRUEs per sample.
lhon_detected_genes <- colSums(lhon_counts_all > 0)

# Display detected-gene counts in the original sample order.
lhon_detected_genes

# Display them from smallest to largest without rearranging the data.
sort(lhon_detected_genes)

# User-confirmed: detected rows range from YC-2 = 15,851 to WT-1 = 18,414.
# All samples have fewer than half of the 49,714 rows detected, so most
# rows are zero in each sample. Interpret the initial all-gene plots accordingly.
# This counts positive gene rows; it does not filter genes or change counts.

# STEP 7C: Initial all-gene log-CPM boxplot, following the RP/glaucoma pattern.
# CPM scales by each sample's total counts; log2 compresses the wide range.
# prior.count = 2 is library-size scaled so zero counts have finite log values.
# This is the initial view, before expression filtering and TMM scaling.
lhon_log_cpm <- edgeR::cpm(
  lhon_counts_all,
  log = TRUE,
  prior.count = 2
)

# Define the plot once for consistent on-screen display and PNG export.
lhon_draw_boxplot <- function() {
  boxplot(
    lhon_log_cpm,
    las = 2,
    col = "lightblue",
    outline = FALSE,
    main = "LHON retinal libraries: expression distributions",
    xlab = "Sample",
    ylab = "Log2 counts per million (CPM)"
  )
}
# outline = FALSE hides individual outlier dots; it does not remove gene rows.

# Display the boxplot in R.
lhon_draw_boxplot()

# Set the output filename in the LHON folder.
lhon_boxplot_file <- file.path(
  lhon_output_dir,
  "LHON_GSE304707_logCPM_boxplot.png"
)

# Draw the same plot into a high-resolution PNG, close it and confirm saving.
png(lhon_boxplot_file, width = 2400, height = 1800, res = 300)
lhon_draw_boxplot()
dev.off()
file.exists(lhon_boxplot_file)

# Full target: results/lhon/LHON_GSE304707_logCPM_boxplot.png
# User supplied the saved PNG: the all-gene boxplot was generated.
# Median lines overlap the lower box edges at the zero-count log-CPM floor;
# upper box edges are modestly lower for Y-3 and the YC libraries, especially YC-2.
# These are descriptive observations; no library is excluded on this figure.
# CPM options were checked against installed edgeR 4.8.2 help and the official
# RNAseq123 transformation workflow; no user counts were loaded by the assistant.
# This preparation retains all 49,714 rows and all 12 samples. The saved image
# is not an R workspace checkpoint.

# STEP 7D: Density plot of the same all-gene log-CPM values, as in RP.
# Estimate one smooth expression distribution per sample.
# Use the existing lhon_log_cpm; no new filtering or TMM is performed here.
lhon_density_curves <- lapply(
  seq_len(ncol(lhon_log_cpm)),
  function(i) density(lhon_log_cpm[, i])
)
names(lhon_density_curves) <- colnames(lhon_log_cpm)

# Assign a distinct colour to each sample curve.
lhon_density_colors <- hcl.colors(ncol(lhon_log_cpm), palette = "Dark 3")

# Define the figure once for consistent screen display and saving.
lhon_draw_density <- function() {
  # Set axis limits using all curves so none is clipped.
  lhon_xlim <- range(sapply(lhon_density_curves, function(d) range(d$x)))
  lhon_ylim <- c(0, max(sapply(lhon_density_curves, function(d) max(d$y))))

  # Draw the first sample's curve and label the axes.
  plot(
    lhon_density_curves[[1]],
    xlim = lhon_xlim, ylim = lhon_ylim,
    col = lhon_density_colors[1], lwd = 2,
    main = "LHON retinal libraries: expression density",
    xlab = "Log2 counts per million (CPM)",
    ylab = "Density"
  )

  # Add the remaining samples to the same plot.
  for (i in 2:length(lhon_density_curves)) {
    lines(lhon_density_curves[[i]], col = lhon_density_colors[i], lwd = 2)
  }

  # Identify each sample's curve in the legend.
  legend(
    "topright",
    legend = names(lhon_density_curves),
    col = lhon_density_colors, lty = 1, lwd = 2,
    bty = "n", cex = 0.75
  )
}

# Display the density plot in R.
lhon_draw_density()

# Save the same figure in the LHON output folder.
lhon_density_file <- file.path(
  lhon_output_dir,
  "LHON_GSE304707_logCPM_density.png"
)
png(lhon_density_file, width = 2400, height = 1800, res = 300)
lhon_draw_density()
dev.off()
file.exists(lhon_density_file)

# Full target: results/lhon/LHON_GSE304707_logCPM_density.png
# User supplied the saved density PNG; figure creation is confirmed.
# Curves broadly overlap, with a dominant peak near the zero-count floor,
# a broad positive-expression hump near log-CPM 4-5, and a thin upper tail.
# YC-2 has a taller low-value peak, consistent with its fewer detected rows;
# smoothing also influences peak heights. No library is excluded on this plot.
# Each curve has approximately unit area; its height shows concentration
# along the log-CPM axis, not a raw gene count or total expression level.
# Expect a prominent low-value peak because many rows have count zero.
# Default per-sample smoothing matches RP; smoothed curves can extend below
# the observed minimum. Such extensions are not additional measured values.
# Reference: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# Interpretation uses the user's figure and prior summaries; no additional
# expression calculations or plots were performed by the assistant.
# Save a new checkpoint only when a separate save block is provided/requested.

# STEP 8A: Inspect gene labels before preparing gene annotations.
# The observed CSV has gene_name plus 12 sample columns. Unlike the RP
# workbook, it has no additional supplied gene-annotation columns.
# Follow glaucoma's label inventory before choosing how to annotate/map rows.

# Retrieve the labels already attached to the count matrix, in gene-row order.
lhon_gene_labels <- rownames(lhon_counts_all)

# Inspect ten label examples.
head(lhon_gene_labels, n = 10)

# Check which labels begin with the mouse Ensembl gene-ID prefix.
# ^ means "starts with"; this inventories a prefix, not database validity.
lhon_ensembl_prefix <- grepl("^ENSMUSG", lhon_gene_labels)
table(Ensembl_prefix = lhon_ensembl_prefix)

# Inspect ten examples that use another label format.
# A FALSE prefix result is a format observation, not a failed count check.
head(lhon_gene_labels[!lhon_ensembl_prefix], n = 10)

# Confirm matrix gene labels match the original gene_name column in order.
identical(lhon_gene_labels, lhon_raw$gene_name)

# User-confirmed: all 49,714 labels return FALSE for the ENSMUSG prefix,
# and exact label/order alignment to gene_name is TRUE. Examples include
# Gm10568, Gm7357, Gm6085, Sox17, Gm37587, Gm6101, Rp1, Gm37483,
# Mrpl15 and 4933401J01Rik. These are symbol-style examples; this prefix
# inventory does not validate every label against a reference database.
# Supplied labels are preserved in the gene-annotation table in Step 8B.
# Expression filtering uses numeric counts/library sizes/sample support, not
# identifier spelling. Do not exclude rows simply for lacking ENSMUSG.
# At this identifier-inspection stage no expression filter had been applied.
# Supplied labels stay unchanged throughout the subsequent QC steps.

# STEP 8B: Preserve the supplied gene names as the annotation table.
# Keep the deposited gene_name field exactly as supplied.
lhon_gene_annotation <- data.frame(
  gene_name = lhon_gene_labels,
  stringsAsFactors = FALSE
)

# Check the table size and inspect ten supplied names.
dim(lhon_gene_annotation)
head(lhon_gene_annotation, n = 10)

# Confirm the annotation follows the count matrix's gene-row order.
identical(lhon_gene_annotation$gene_name, rownames(lhon_counts_all))

# User-confirmed: 49,714 rows x 1 column and gene-row alignment TRUE.

# STEP 9: Assemble counts, sample metadata and gene annotations, as in RP.
# Check both existing alignments before combining the components.
stopifnot(
  identical(lhon_metadata_matched$sample_id, colnames(lhon_counts_all)),
  identical(lhon_gene_annotation$gene_name, rownames(lhon_counts_all))
)

# Record the four genotype categories in a consistent display order.
lhon_qc_group <- factor(
  lhon_metadata_matched$genotype,
  levels = c("WT", "YARS2", "COI", "YARS2/COI")
)
stopifnot(!anyNA(lhon_qc_group))
table(lhon_qc_group)

# Store raw counts and the two descriptive tables in one coordinated object.
lhon_dge <- edgeR::DGEList(
  counts = lhon_counts_all,
  samples = lhon_metadata_matched,
  group = lhon_qc_group,
  genes = lhon_gene_annotation
)

# Inspect the dimensions and sample information.
dim(lhon_dge)
lhon_dge$samples[, c("sample_id", "group", "sex", "lib.size", "norm.factors")]

# Check sample and gene order inside the assembled object.
c(
  sample_order = identical(lhon_dge$samples$sample_id, colnames(lhon_dge)),
  gene_order = identical(lhon_dge$genes$gene_name, rownames(lhon_dge))
)

# User-confirmed: 3 libraries per genotype, 49,714 x 12, all initial factors
# equal 1, and sample_order/gene_order both TRUE.
# lib.size holds each sample's total counts; norm.factors initially equals 1.
# DGEList retains all rows by default. This organizes existing data; expression
# filtering and TMM are later steps. This assembly is now confirmed complete.
# Constructor options were checked in the installed official edgeR help.
# The assistant did not load counts or execute this block; no checkpoint update.

# STEP 10: Filter genes with insufficient expression for exploratory QC.
# Use the same documented count settings as RP/glaucoma. Group sizes are
# three libraries each; qualifying samples need not all be in one genotype.
lhon_keep_genes <- edgeR::filterByExpr(
  lhon_dge,
  group = lhon_qc_group,
  min.count = 10,
  min.total.count = 15
)

# TRUE = retain the row; FALSE = exclude it from the filtered QC copy.
table(Retain_for_QC = lhon_keep_genes)

# Create a separate filtered object, retaining every sample column.
# Recalculate library totals over retained rows; retained counts stay unchanged.
lhon_dge_qc <- lhon_dge[
  lhon_keep_genes,
  ,
  keep.lib.sizes = FALSE
]

# Check the retained gene count and confirm all 12 samples remain.
dim(lhon_dge_qc)

# Confirm metadata and gene annotations remain correctly aligned.
c(
  sample_order = identical(
    lhon_dge_qc$samples$sample_id,
    lhon_metadata_matched$sample_id
  ),
  gene_order = identical(
    lhon_dge_qc$genes$gene_name,
    rownames(lhon_dge_qc)
  )
)

# User-confirmed: 17,860 rows retained and 31,854 excluded from this QC copy;
# dimensions 17,860 x 12, sample_order TRUE and gene_order TRUE.
# Preserve lhon_dge and lhon_counts_all as the complete inputs.

# STEP 11: TMM composition scaling and normalized log-CPM, as in RP.
# Estimate one library-scaling factor per sample from the filtered counts.
# TMM requires a sufficiently stable/balanced reference set of genes.
lhon_dge_qc <- edgeR::normLibSizes(lhon_dge_qc, method = "TMM")

# Display sample IDs, groups, retained library totals and estimated factors.
lhon_dge_qc$samples[, c("sample_id", "group", "lib.size", "norm.factors")]

# Require finite, positive factors before calculating normalized log-CPM.
stopifnot(
  all(is.finite(lhon_dge_qc$samples$norm.factors)),
  all(lhon_dge_qc$samples$norm.factors > 0)
)

# Create a separate filtered/TMM log-expression matrix for exploration.
# Keep the original all-gene lhon_log_cpm used for boxplot/density inspection.
lhon_log_cpm_qc <- edgeR::cpm(
  lhon_dge_qc,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

# Verify the normalized matrix size, finite values and sample order.
dim(lhon_log_cpm_qc)
c(
  finite_logCPM = all(is.finite(lhon_log_cpm_qc)),
  sample_order = identical(
    colnames(lhon_log_cpm_qc), lhon_dge_qc$samples$sample_id
  )
)
# Final-save verification confirms 17,860 x 12 and valid alignment/finite values.
# Saved TMM factors range from 0.8932679114 (YC-2) to 1.0394518134 (Y-1).
# Effective library sizes are lib.size * norm.factors. Retained count entries
# remain unchanged. No user data were normalized by the assistant.

# STEP 12: Create both MDS and PCA before interpreting either figure.
# Load the small RP-style formatting helper saved in this LHON folder.
# It defines a plot function with points, separated labels and outside legends.
source(file.path("scripts", "07_LHON_projection_style.R"))

# MDS: each sample pair uses its own 500 most divergent retained genes.
lhon_mds <- limma::plotMDS(
  lhon_log_cpm_qc,
  top = 500,
  gene.selection = "pairwise",
  plot = FALSE
)

# PCA: one common set of 500 most variable genes across all samples.
# This is limma's documented PCA mode, despite the plotMDS function name.
lhon_pca <- limma::plotMDS(
  lhon_log_cpm_qc,
  top = 500,
  gene.selection = "common",
  plot = FALSE
)

# Apply the same style and genotype colours to both stored projections.
# Axis percentages come from each calculation, not manually chosen values.
lhon_mds_plot <- lhon_projection_plot(
  lhon_mds, lhon_dge_qc$samples,
  title = "LHON retinal libraries: MDS after TMM",
  axis_names = c("Leading logFC dim 1", "Leading logFC dim 2")
)
lhon_pca_plot <- lhon_projection_plot(
  lhon_pca, lhon_dge_qc$samples,
  title = "LHON retinal libraries: PCA after TMM",
  axis_names = c("PC1", "PC2")
)

# Define easy redraw commands and display both figures.
lhon_draw_mds <- function() print(lhon_mds_plot)
lhon_draw_pca <- function() print(lhon_pca_plot)
lhon_draw_mds()
lhon_draw_pca()

# Set the two output filenames in the LHON folder.
lhon_mds_file <- file.path(lhon_output_dir, "LHON_GSE304707_TMM_MDS.png")
lhon_pca_file <- file.path(lhon_output_dir, "LHON_GSE304707_TMM_PCA_top500.png")

# Save MDS and close its image file.
png(lhon_mds_file, width = 2400, height = 1800, res = 300)
lhon_draw_mds()
dev.off()

# Save PCA and close its image file.
png(lhon_pca_file, width = 2400, height = 1800, res = 300)
lhon_draw_pca()
dev.off()

# Confirm both PNG files exist before reviewing the plots together.
c(MDS = file.exists(lhon_mds_file), PCA = file.exists(lhon_pca_file))
# User supplied both PNGs: PCA axes are approximately 26% and 13%;
# MDS axes approximately 18% and 11%. C-2/C-3 and WT-1/WT-3 are nearby pairs;
# WT-2 and Y-1 lie lower than most libraries. Genotype clustering is incomplete.
# Numerical factors were exported and checked during final checkpoint verification.
# M labels denote male. All 12 libraries remain retained.
# The common-gene PCA and pairwise MDS settings match RP/glaucoma; group labels
# colour the results without selecting genes or determining coordinates.
# The assistant prepared these instructions only; no real-data normalization,
# projection, figure export or updated workspace checkpoint was executed.


# CLOSURE — 30 October 2026, at the user's requested endpoint.
# Retain all 12 libraries; no failed library was established by this precheck.
# The QC copy contains 17,860 genes. No DE, pathway or flight-similarity tests
# have been performed. Save the final active-session objects with:
# source(file.path("scripts", "08_LHON_GSE304707_save_final.R"))
# Use 06_LHON_GSE304707_reload.R afterward, not this whole teaching script.
# The comprehensive handoff is saved in Downloads and archived in this folder.
# Final checkpoint is saved and restore-verified: 47 LHON objects,
# LHON_GSE304707_precheck_checkpoint_2026-10-30.RData (7,258,085 bytes).
