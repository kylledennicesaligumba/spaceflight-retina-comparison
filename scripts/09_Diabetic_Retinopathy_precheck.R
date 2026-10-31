# Portable repository copy; start R in the repository root.
# This is a processed MICROARRAY precheck: no additional log2 transformation,
# CPM calculation or TMM normalization is applied to the supplied values.
# DR row labels are Entrez IDs; supplied platform annotation links them to gene symbols.
# Original working files are preserved; GitHub sync did not execute this code.
# Capture the project folder before the original workflow changes working directory.
dr_project_dir <- normalizePath(getwd(), mustWork = TRUE)
# Set the dataset and output folders.
dr_dir <- Sys.getenv("DR_DATA_DIR", unset = file.path(dr_project_dir, "data", "raw", "Diabetic_Retinopathy_GSE87433"))
dr_output_dir <- file.path(dr_project_dir, "results", "diabetic_retinopathy")
dir.create(dr_output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the output folder as the default save location.
setwd(dr_output_dir)
getwd()

# List the downloaded files so we can choose what to unpack together.
list.files(dr_dir)

# List the family archive contents before extracting them.
dr_metadata_archive <- file.path(dr_dir, "GSE87433_family")
untar(dr_metadata_archive, list = TRUE)

# Locate the downloaded expression matrix.
dr_expression_file <- file.path(dr_dir, "GSE87433_series_matrix")

# Load the tab-separated table, skipping GEO description lines starting with !.
dr_raw <- read.delim(
  dr_expression_file,
  comment.char = "!",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the table size, column names and first five rows.
dim(dr_raw)
colnames(dr_raw)
head(dr_raw, 5)

# User-confirmed: 24,973 rows; ID_REF plus 17 GSM sample columns.
# The displayed expression values include decimals; processing is not yet confirmed.

# Check for missing entries anywhere in the imported table.
anyNA(dr_raw)

# Check for missing or blank feature identifiers.
sum(is.na(dr_raw$ID_REF) | trimws(dr_raw$ID_REF) == "")

# Check whether any feature identifier appears more than once.
anyDuplicated(dr_raw$ID_REF)

# Check that all sample columns are stored as numeric values.
sapply(dr_raw[, -1, drop = FALSE], is.numeric)

# User-confirmed: no missing entries or missing/blank IDs, no duplicate IDs,
# and all 17 sample columns are numeric.

# Build the expression matrix, keeping feature IDs as row labels.
dr_expr_all <- as.matrix(dr_raw[, -1, drop = FALSE])
rownames(dr_expr_all) <- dr_raw$ID_REF

# Check the number of features and samples in the matrix.
dim(dr_expr_all)

# Check that every expression value is finite.
all(is.finite(dr_expr_all))

# Inspect the smallest and largest expression values.
range(dr_expr_all)

# User-confirmed: expression matrix is 24,973 x 17; all values are finite.
# Observed expression range: 5.066406 to 14.095703.

# Extract the XML metadata into its own output subfolder.
dr_metadata_dir <- file.path(dr_output_dir, "metadata")
dir.create(dr_metadata_dir, recursive = TRUE, showWarnings = FALSE)
untar(dr_metadata_archive, files = "GSE87433_family.xml", exdir = dr_metadata_dir)

# Confirm that R can find the extracted metadata file.
dr_metadata_file <- file.path(dr_metadata_dir, "GSE87433_family.xml")
file.exists(dr_metadata_file)

# User-confirmed: the extracted XML metadata file exists.

# Read the extracted XML metadata and display its structure.
dr_metadata_xml <- xml2::read_xml(dr_metadata_file)
dr_metadata_xml

# Locate sample records using the XML namespace shown in the overview.
dr_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")
dr_sample_nodes <- xml2::xml_find_all(
  dr_metadata_xml, "/geo:MINiML/geo:Sample", ns = dr_ns
)

# Create a table linking each GSM identifier to its sample title.
dr_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(dr_sample_nodes, "iid"),
  title = xml2::xml_text(
    xml2::xml_find_first(dr_sample_nodes, "geo:Title", ns = dr_ns)
  ),
  stringsAsFactors = FALSE
)
dr_metadata_all

# User-confirmed sample titles: 6 normal control, 6 diabetic control,
# and 5 diabetic with transplantation. These are the submitted labels.

# Read the processing descriptions stored for the first sample.
dr_processing_first <- xml2::xml_text(
  xml2::xml_find_all(
    dr_sample_nodes[[1]], ".//geo:Data-Processing", ns = dr_ns
  ),
  trim = TRUE
)

# Print the descriptions as readable text.
cat(dr_processing_first, sep = "\n\n")

# User-reviewed first-sample processing: JMP Genomics 7; Affymetrix human
# Mogene-2.0-type arrays; custom CDF v19 with Entrez-based gene definitions;
# fluorescence intensities processed using quantile normalization,
# RMA background correction and median-polish probeset summarization.
# The excerpt does not explicitly state the log base or the sample organism.

# Add the recorded organism for each sample.
dr_metadata_all$organism <- xml2::xml_text(
  xml2::xml_find_first(dr_sample_nodes, ".//geo:Organism", ns = dr_ns),
  trim = TRUE
)

# Add the recorded sample source to help identify the tissue or material.
dr_metadata_all$source <- xml2::xml_text(
  xml2::xml_find_first(dr_sample_nodes, ".//geo:Source", ns = dr_ns),
  trim = TRUE
)

# Display the organism and source beside each sample ID.
dr_metadata_all[, c("geo_accession", "organism", "source")]

# User-confirmed: all 17 sample organism fields are Mus musculus;
# sources describe whole retina: 6 normoglycemic, 6 diabetic and 5 treated
# with transplantation and insulin-releasing pellets. Sources mention @12 weeks;
# the timing reference has not yet been established.
# Check the platform record because first-sample processing says human array.

# Locate the array platform records in the XML metadata.
dr_platform_nodes <- xml2::xml_find_all(
  dr_metadata_xml, "/geo:MINiML/geo:Platform", ns = dr_ns
)

# Display the platform ID, title and organism to check the array description.
dr_platform_info <- data.frame(
  platform_id = xml2::xml_attr(dr_platform_nodes, "iid"),
  title = xml2::xml_text(
    xml2::xml_find_first(dr_platform_nodes, "geo:Title", ns = dr_ns)
  ),
  organism = xml2::xml_text(
    xml2::xml_find_first(dr_platform_nodes, ".//geo:Organism", ns = dr_ns)
  ),
  stringsAsFactors = FALSE
)
dr_platform_info

# User-confirmed platform: GPL20710, Affymetrix Mouse Gene 2.0 ST Array
# (mogene20st_Mm_ENTREZG_19), organism Mus musculus.
# This agrees with sample organism fields. The earlier processing description's
# reference to a human array appears to be a wording error; retain this discrepancy.

# Find the biological characteristic entries recorded across the samples.
dr_characteristic_nodes <- xml2::xml_find_all(
  dr_sample_nodes, ".//geo:Characteristics", ns = dr_ns
)

# List the characteristic field names and their number of occurrences.
dr_characteristic_tags <- xml2::xml_attr(dr_characteristic_nodes, "tag")
table(dr_characteristic_tags, useNA = "ifany")

# User-confirmed characteristic tags: background strain, beta cells
# transplantation, disease, streptozotocin and tissue; 17 occurrences each.

# Read a named characteristic from each sample record.
dr_get_characteristic <- function(tag) {
  dr_path <- paste0(".//geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(dr_sample_nodes, dr_path, ns = dr_ns),
    trim = TRUE
  )
}

# Add the five recorded characteristics to the sample table.
dr_metadata_all$strain <- dr_get_characteristic("background strain")
dr_metadata_all$transplantation <- dr_get_characteristic("beta cells transplantation")
dr_metadata_all$disease <- dr_get_characteristic("disease")
dr_metadata_all$streptozotocin <- dr_get_characteristic("streptozotocin")
dr_metadata_all$tissue <- dr_get_characteristic("tissue")

# Display the extracted values beside each GSM sample ID.
dr_metadata_all[, c(
  "geo_accession", "strain", "transplantation", "disease", "streptozotocin", "tissue"
)]

# User-confirmed characteristics: strain is recorded as C57Bl6/J and tissue
# as Whole retina for all 17 samples. Six normoglycemic samples have
# transplantation=no and streptozotocin=Ctrl; six diabetic samples have
# transplantation=no and streptozotocin=Yes; five diabetic samples have
# transplantation=yes and streptozotocin=Yes.
# Sex, age and collection/extraction timing are not established by these fields.

# Read the RNA-extraction description recorded for each sample.
dr_extraction_protocol <- xml2::xml_text(
  xml2::xml_find_first(
    dr_sample_nodes, ".//geo:Extract-Protocol", ns = dr_ns
  ),
  trim = TRUE
)
names(dr_extraction_protocol) <- dr_metadata_all$geo_accession

# Print each distinct description once so we can review collection details.
cat(unique(dr_extraction_protocol), sep = "\n\n")
cat("\n")

# User-confirmed extraction protocol: total RNA isolated using Trizol
# according to manufacturer instructions; no timing details in that text.

# Paper review: Friedrichs et al. (2017), DOI 10.1007/s00125-017-4254-y.
# Source: https://link.springer.com/content/pdf/10.1007/s00125-017-4254-y.pdf
# Methods, journal page 1355: male C57BL/6J mice; STZ induction at 8 weeks.
# Expression comparisons use the 12-week endpoint after induction.
# Approximate endpoint age of 20 weeks is inferred (8 + 12), not measured per sample.
# Islet transplantation occurs after 6 weeks of diabetes; endpoint follows 6 weeks later.
# Calendar collection/RNA-extraction dates are not reported in these Methods.
# These are paper-derived annotations; no per-sample fields have been overwritten.

# Read the extracted XML as text so we can search all its fields.
dr_metadata_lines <- readLines(dr_metadata_file, warn = FALSE)

# Search for sex, age and collection or extraction details, including protocols.
dr_metadata_hits <- grep(
  "sex|gender|male|female|age|old|week|month|collect|harvest|extract|date",
  dr_metadata_lines,
  value = TRUE,
  ignore.case = TRUE
)

# Print each matching line once to avoid repeated sample descriptions.
cat(unique(trimws(dr_metadata_hits)), sep = "\n\n")
cat("\n")

# Save the current DR tables, metadata and supporting text for resuming.
dr_checkpoint_file <- file.path(
  dr_output_dir, "DR_GSE87433_precheck_checkpoint_2026-10-30.RData"
)
save(
  dr_raw, dr_expr_all, dr_metadata_all, dr_platform_info,
  dr_processing_first, dr_extraction_protocol,
  dr_metadata_lines, dr_metadata_hits,
  dr_dir, dr_output_dir, dr_expression_file, dr_metadata_archive,
  dr_metadata_dir, dr_metadata_file, dr_ns,
  file = dr_checkpoint_file
)

# Confirm that the checkpoint file was created.
file.exists(dr_checkpoint_file)

# Check that metadata and expression columns contain the same unique GSM IDs.
dr_sample_ids <- colnames(dr_expr_all)
stopifnot(
  !anyNA(dr_metadata_all$geo_accession),
  !anyNA(dr_sample_ids),
  anyDuplicated(dr_metadata_all$geo_accession) == 0L,
  anyDuplicated(dr_sample_ids) == 0L,
  setequal(dr_metadata_all$geo_accession, dr_sample_ids)
)

# Arrange a metadata copy in the same order as the expression-matrix columns.
dr_metadata_matched <- dr_metadata_all[
  match(dr_sample_ids, dr_metadata_all$geo_accession),
  , drop = FALSE
]
rownames(dr_metadata_matched) <- NULL

# Confirm the metadata dimensions and exact sample-ID order.
dim(dr_metadata_matched)
identical(dr_metadata_matched$geo_accession, colnames(dr_expr_all))

# User-confirmed: matched metadata has 17 rows and 9 columns; exact GSM order is TRUE.

# Display the names of all nine matched metadata columns.
colnames(dr_metadata_matched)

# Review all metadata values in the matched sample order.
dr_metadata_matched

# DR is a processed microarray matrix: summed expression signals are not
# sequencing-depth/library-size measurements. Proceed to distribution checks.
# Component reference: Bioconductor Affymetrix microarray workflow (intensity QC).
# https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
# Use supplied values directly for this check; no further transformation is applied.

# Define a boxplot comparing the supplied expression values across samples.
dr_draw_boxplot <- function() {
  dr_previous_par <- par(mar = c(8, 4.5, 3, 1))
  on.exit(par(dr_previous_par))
  boxplot(
    dr_expr_all,
    las = 2,
    cex.axis = 0.7,
    col = "lightblue",
    outline = FALSE,
    main = "DR retinal arrays: expression distributions",
    ylab = "Processed expression value"
  )
}

# Display the plot in R; extreme-value dots are hidden for readability.
dr_draw_boxplot()

# Save the same plot as a high-resolution PNG in the DR output folder.
dr_boxplot_file <- file.path(dr_output_dir, "DR_GSE87433_expression_boxplot.png")
png(dr_boxplot_file, width = 2400, height = 1800, res = 300)
dr_draw_boxplot()
dev.off()
file.exists(dr_boxplot_file)

# User supplied the saved boxplot: all 17 samples have visually very similar
# medians, box heights and whiskers. This is consistent with the reported
# quantile normalization; it does not establish identical gene-wise profiles
# or overall sample quality. No density plot or PCA has yet been reviewed.

# Estimate one smooth expression distribution per sample.
dr_density_curves <- lapply(
  seq_len(ncol(dr_expr_all)),
  function(i) density(dr_expr_all[, i])
)
names(dr_density_curves) <- colnames(dr_expr_all)
dr_density_colors <- hcl.colors(ncol(dr_expr_all), palette = "Dark 3")

# Define the figure once for consistent display and saving, as in LHON.
dr_draw_density <- function() {
  dr_xlim <- range(sapply(dr_density_curves, function(d) range(d$x)))
  dr_ylim <- c(0, max(sapply(dr_density_curves, function(d) max(d$y))))

  plot(
    dr_density_curves[[1]],
    xlim = dr_xlim, ylim = dr_ylim,
    col = dr_density_colors[1], lwd = 2,
    main = "DR retinal arrays: expression density",
    xlab = "Processed expression value", ylab = "Density"
  )

  # Add the remaining samples and identify their colours in the legend.
  for (i in 2:length(dr_density_curves)) {
    lines(dr_density_curves[[i]], col = dr_density_colors[i], lwd = 2)
  }
  legend(
    "topright", legend = names(dr_density_curves),
    col = dr_density_colors, lty = 1, lwd = 2,
    bty = "n", cex = 0.65
  )
}

# Display the density plot in R.
dr_draw_density()

# Save the same figure in the DR output folder.
dr_density_file <- file.path(dr_output_dir, "DR_GSE87433_expression_density.png")
png(dr_density_file, width = 2400, height = 1800, res = 300)
dr_draw_density()
dev.off()
file.exists(dr_density_file)

# User supplied the full density plot: curves visually overlap almost completely.
# This close-up changes the displayed axis window and line width only.

# Magnify the peak using the existing full-data density estimates.
dr_draw_density_zoom <- function() {
  dr_previous_par <- par(mar = c(5, 4.5, 3, 9))
  on.exit(par(dr_previous_par))
  plot(
    dr_density_curves[[1]],
    xlim = c(5.2, 5.6), ylim = c(0.75, 0.90),
    col = dr_density_colors[1], lwd = 1,
    main = "DR expression density: peak close-up",
    xlab = "Processed expression value", ylab = "Density"
  )
  for (i in 2:length(dr_density_curves)) {
    lines(dr_density_curves[[i]], col = dr_density_colors[i], lwd = 1)
  }

  # Place the legend in the right margin to keep the curves visible.
  legend(
    x = par("usr")[2], y = par("usr")[4],
    legend = names(dr_density_curves),
    col = dr_density_colors, lty = 1, lwd = 1,
    bty = "n", cex = 0.65, xpd = NA
  )
}

# Display the close-up in R.
dr_draw_density_zoom()

# Save the close-up as a separate figure in the DR output folder.
dr_density_zoom_file <- file.path(
  dr_output_dir, "DR_GSE87433_expression_density_zoom.png"
)
png(dr_density_zoom_file, width = 2400, height = 1800, res = 300)
dr_draw_density_zoom()
dev.off()
file.exists(dr_density_zoom_file)

# User supplied the zoomed density figure; curves still visually overlap.
# Saved code computes density from column i and draws curves 1 through 17.
# Verify the actual R objects before concluding the curves coincide numerically.

# Confirm that there is one correctly named curve for each sample column.
c(samples = ncol(dr_expr_all), curves = length(dr_density_curves))
identical(names(dr_density_curves), colnames(dr_expr_all))

# Compare each stored curve with a fresh calculation from its own sample.
dr_density_verified <- sapply(seq_len(ncol(dr_expr_all)), function(i) {
  dr_check_curve <- density(dr_expr_all[, i])
  isTRUE(all.equal(
    dr_density_curves[[i]][c("x", "y")],
    dr_check_curve[c("x", "y")]
  ))
})

# Display the verification result beside each GSM sample ID.
names(dr_density_verified) <- colnames(dr_expr_all)
dr_density_verified

# User's density-verification output shows TRUE for all 17 sample curves.
# These checks compare each curve with its own input, not with other samples.
# Begin identifier review before selecting an annotation mapping.
# Platform/processing descriptions report a custom CDF v19 with Entrez-based
# gene definitions; a digits-only pattern does not establish database validity.

# Retrieve the supplied feature IDs in expression-matrix row order.
dr_gene_ids <- rownames(dr_expr_all)

# Inspect ten examples and check whether the IDs contain digits only.
head(dr_gene_ids, 10)
table(Digits_only = grepl("^[0-9]+$", dr_gene_ids))

# Confirm that matrix row IDs match the original ID_REF column as text.
identical(dr_gene_ids, as.character(dr_raw$ID_REF))

# User-confirmed: all 24,973 supplied IDs contain digits only and the matrix
# row labels match as.character(dr_raw$ID_REF) in order (TRUE).
# This establishes format and alignment, not database validity or gene symbols.
# Inspect the supplied custom-platform table before choosing a gene mapping.

# Extract the supplied platform table into an annotation output subfolder.
dr_annotation_dir <- file.path(dr_output_dir, "annotation")
dir.create(dr_annotation_dir, recursive = TRUE, showWarnings = FALSE)
untar(dr_metadata_archive, files = "GPL20710-tbl-1.txt", exdir = dr_annotation_dir)

# Confirm that R can find the extracted platform table.
dr_platform_annotation_file <- file.path(dr_annotation_dir, "GPL20710-tbl-1.txt")
file.exists(dr_platform_annotation_file)

# User-confirmed: the extracted GPL20710 platform table exists (TRUE).

# Import corrected after the user's preview showed 24,972 rows and the first
# gene's annotation values as column names. This file has no header row.
# R's temporary V1...V15 labels will be resolved using column descriptions.
# Read the headerless platform table, keeping the first gene as a data row.
dr_platform_annotation <- read.delim(
  dr_platform_annotation_file,
  header = FALSE,
  colClasses = "character",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check the full table size and preview the first three columns.
dim(dr_platform_annotation)
head(dr_platform_annotation[, 1:3], 5)

# User-confirmed corrected platform import: 24,973 x 15; first row starts
# V1=11287, V2=11287, V3=Pzp. V3 previews include Pzp, Aanat, Aatk, Abca1, Abca4.
# Read the supplied column definitions before assigning biological field names.

# Find the column definitions for the platform table in the XML.
dr_platform_column_nodes <- xml2::xml_find_all(
  dr_platform_nodes, ".//geo:Data-Table/geo:Column", ns = dr_ns
)

# Display the supplied column names beside their table positions.
dr_platform_columns <- data.frame(
  position = xml2::xml_attr(dr_platform_column_nodes, "position"),
  name = xml2::xml_text(
    xml2::xml_find_first(dr_platform_column_nodes, "geo:Name", ns = dr_ns)
  ),
  stringsAsFactors = FALSE
)
dr_platform_columns

# User-confirmed XML definitions list all 15 columns in positions 1 through 15.
# ID is position 1, ENTREZ_GENE_ID position 2, Symbol position 3,
# dbXrefs position 6 and description position 9.

# Apply the supplied XML column names to the platform table.
colnames(dr_platform_annotation) <- dr_platform_columns$name

# Preview the original ID, Entrez ID, gene symbol and description together.
head(dr_platform_annotation[, c(
  "ID", "ENTREZ_GENE_ID", "Symbol", "description"
)], 5)

# User's named-table preview confirms ID, ENTREZ_GENE_ID, Symbol and description
# are accessible; the first records include 11287/Pzp and 11298/Aanat.
# Match the supplied platform ID to the expression row ID before projections.

# Check that the platform table has unique IDs matching the expression rows.
stopifnot(
  !anyNA(dr_platform_annotation$ID),
  anyDuplicated(dr_platform_annotation$ID) == 0L,
  setequal(rownames(dr_expr_all), dr_platform_annotation$ID)
)

# Arrange the supplied gene annotations in expression-matrix row order.
dr_gene_annotation <- dr_platform_annotation[
  match(rownames(dr_expr_all), dr_platform_annotation$ID),
  , drop = FALSE
]
rownames(dr_gene_annotation) <- NULL

# Confirm the annotation dimensions and exact gene-ID order.
dim(dr_gene_annotation)
identical(dr_gene_annotation$ID, rownames(dr_expr_all))

# User-confirmed annotation alignment: 24,973 x 15 and exact gene-ID order TRUE.
# Apply the RP/LHON pairwise-MDS/common-gene-PCA settings to supplied DR values.
# There is no additional transformation, normalization or row filtering here.
# The exact log base is not explicitly recorded in the deposited description;
# MDS axes are therefore labelled generically rather than as measured log2 fold changes.
# Method: limma plotMDS documentation, https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

# Load the familiar projection style and confirm sample and gene alignment.
source(file.path(dr_project_dir, "scripts", "11_DR_projection_style.R"))
stopifnot(
  identical(dr_metadata_matched$geo_accession, colnames(dr_expr_all)),
  identical(dr_gene_annotation$ID, rownames(dr_expr_all))
)

# Label the three groups using the disease and transplantation metadata.
dr_group_names <- c(
  "normoglycemic_no" = "Normal control",
  "diabetic_no" = "Diabetic control",
  "diabetic_yes" = "Diabetic + transplantation"
)
dr_group_key <- paste(
  dr_metadata_matched$disease, dr_metadata_matched$transplantation, sep = "_"
)
dr_metadata_matched$group <- factor(
  unname(dr_group_names[dr_group_key]), levels = unname(dr_group_names)
)
stopifnot(!anyNA(dr_metadata_matched$group))
table(dr_metadata_matched$group)

# Calculate MDS using 500 divergent features for each sample pair.
dr_mds <- limma::plotMDS(
  dr_expr_all, top = 500, gene.selection = "pairwise", plot = FALSE
)

# Calculate PCA using one common set of 500 variable features, as in LHON.
dr_pca <- limma::plotMDS(
  dr_expr_all, top = 500, gene.selection = "common", plot = FALSE
)

# Apply consistent points, separated GSM labels and group colours.
dr_mds_plot <- dr_projection_plot(
  dr_mds, dr_metadata_matched,
  title = "DR retinal arrays: MDS", axis_names = c("MDS1", "MDS2")
)
dr_pca_plot <- dr_projection_plot(
  dr_pca, dr_metadata_matched,
  title = "DR retinal arrays: PCA", axis_names = c("PC1", "PC2")
)
print(dr_mds_plot)
print(dr_pca_plot)

# Save both figures in the DR output folder.
dr_mds_file <- file.path(dr_output_dir, "DR_GSE87433_MDS.png")
dr_pca_file <- file.path(dr_output_dir, "DR_GSE87433_PCA_top500.png")
dr_save_projection(dr_mds_plot, dr_mds_file)
dr_save_projection(dr_pca_plot, dr_pca_file)
c(MDS = file.exists(dr_mds_file), PCA = file.exists(dr_pca_file))

# User supplied both saved projections: MDS axes display 34% and 9%;
# PCA axes display 44% and 11%. The treated samples are separated from the
# normal and diabetic controls in both views (axis signs have no biological meaning).
# MDS tends to place diabetic controls above normal controls; some are close.
# Common-top-500 PCA shows substantial normal/diabetic overlap on the left,
# with the treated group on the right. These projections use different feature
# selection strategies, so their two-dimensional arrangements need not match.
# GSM2331277 is peripheral on MDS2; GSM2331273 is peripheral on PCA PC2.
# These images alone do not justify sample exclusion or establish a technical fault.
# Treated-group separation is not proof of recovery, worse disease, or memory.
# Primary disease-signature comparison recommended from experimental design:
# diabetic controls versus normal controls (6 samples each), with the five
# transplantation/insulin samples treated as a separate condition. No model fitted.

# Save the DR data, matched tables and completed plot results for resuming.
dr_checkpoint_file <- file.path(
  dr_output_dir, "DR_GSE87433_precheck_checkpoint_2026-10-31.RData"
)
save(
  dr_raw, dr_expr_all, dr_metadata_all, dr_metadata_matched,
  dr_gene_annotation, dr_platform_annotation, dr_platform_columns,
  dr_platform_info, dr_density_curves, dr_mds, dr_pca,
  dr_mds_plot, dr_pca_plot,
  dr_processing_first, dr_extraction_protocol,
  dr_metadata_lines, dr_metadata_hits,
  dr_dir, dr_output_dir, dr_expression_file, dr_metadata_archive,
  dr_metadata_dir, dr_metadata_file, dr_ns,
  dr_annotation_dir, dr_platform_annotation_file,
  file = dr_checkpoint_file
)

# Confirm that the new checkpoint exists; TRUE means the file was created.
file.exists(dr_checkpoint_file)
