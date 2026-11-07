# Set the downloaded-data and prequality-check output folders.
amd_dir <- "/Volumes/MaxineDrive/INTERNSHIP/datasets/AMD_GSE38671"
amd_output_dir <- "/Volumes/MaxineDrive/INTERNSHIP/prequalitychecks/amd"
dir.create(amd_output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the output folder as R's default save location.
setwd(amd_output_dir)
getwd()

# List the downloaded files before choosing what to open together.
list.files(amd_dir)

# Import the microarray expression table, skipping GEO description lines.
amd_expression_file <- file.path(amd_dir, "GSE38671_series_matrix")
amd_raw <- read.delim(
  amd_expression_file,
  comment.char = "!",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the table size, column names and first five rows.
dim(amd_raw)
colnames(amd_raw)
head(amd_raw, 5)

# User-confirmed: 35,556 feature rows, ID_REF and 32 GSM sample columns.
# The displayed values include decimals; processing remains to be reviewed.

# Check for missing entries anywhere in the imported table.
anyNA(amd_raw)

# Check for missing or blank feature identifiers.
sum(is.na(amd_raw$ID_REF) | trimws(amd_raw$ID_REF) == "")

# Check whether any feature identifier appears more than once.
anyDuplicated(amd_raw$ID_REF)

# Check that all sample columns are stored as numeric values.
sapply(amd_raw[, -1, drop = FALSE], is.numeric)

# User-confirmed: no missing entries or blank/repeated feature IDs;
# all 32 sample columns are numeric.

# Build the expression matrix, keeping feature IDs as row labels.
amd_expr_all <- as.matrix(amd_raw[, -1, drop = FALSE])
rownames(amd_expr_all) <- amd_raw$ID_REF

# Check the number of features and samples in the matrix.
dim(amd_expr_all)

# Check that every expression value is finite.
all(is.finite(amd_expr_all))

# Inspect the smallest and largest expression values.
range(amd_expr_all)

# User-confirmed: expression matrix is 35,556 x 32; all values are finite.
# Observed expression range: 5.8 to 16,719.9.

# List the metadata archive contents before extracting them.
amd_metadata_archive <- file.path(amd_dir, "GSE38671_family")
untar(amd_metadata_archive, list = TRUE)

# User-confirmed: the archive includes GSE38671_family.xml.

# Extract the XML metadata into its own output subfolder.
amd_metadata_dir <- file.path(amd_output_dir, "metadata")
dir.create(amd_metadata_dir, recursive = TRUE, showWarnings = FALSE)
untar(
  amd_metadata_archive,
  files = "GSE38671_family.xml",
  exdir = amd_metadata_dir
)

# Confirm that R can find the extracted metadata file.
amd_metadata_file <- file.path(amd_metadata_dir, "GSE38671_family.xml")
file.exists(amd_metadata_file)

# User-confirmed: the extracted XML metadata file exists.

# Read the extracted XML metadata and display its structure.
amd_metadata_xml <- xml2::read_xml(amd_metadata_file)
amd_metadata_xml

# User-confirmed: the XML loaded and shows the GEO MINiML namespace.

# Use the XML namespace to locate the sample records.
amd_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")
amd_sample_nodes <- xml2::xml_find_all(
  amd_metadata_xml, "/geo:MINiML/geo:Sample", ns = amd_ns
)
length(amd_sample_nodes)

# Create a table linking each GSM identifier to its sample title.
amd_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(amd_sample_nodes, "iid"),
  title = xml2::xml_text(
    xml2::xml_find_first(amd_sample_nodes, "geo:Title", ns = amd_ns)
  ),
  stringsAsFactors = FALSE
)
amd_metadata_all

# User-confirmed titles: 16 NR and 16 RPE/choroid samples, with young/aged
# and WT/Cfh-/- labels; four labelled replicates per tissue/age/genotype group.

# Read the processing descriptions stored for the first sample.
amd_processing_first <- xml2::xml_text(
  xml2::xml_find_all(
    amd_sample_nodes[[1]], ".//geo:Data-Processing", ns = amd_ns
  ),
  trim = TRUE
)

# Print the descriptions as readable text.
cat(amd_processing_first, sep = "\n\n")

# User-reviewed first-sample processing: Expression Console 1.1 (Affymetrix).
# This field does not specify normalization or a log transformation.

# Add the recorded organism and sample source to the metadata table.
amd_metadata_all$organism <- xml2::xml_text(
  xml2::xml_find_first(amd_sample_nodes, ".//geo:Organism", ns = amd_ns),
  trim = TRUE
)
amd_metadata_all$source <- xml2::xml_text(
  xml2::xml_find_first(amd_sample_nodes, ".//geo:Source", ns = amd_ns),
  trim = TRUE
)

# Review the organism and source beside each sample ID.
amd_metadata_all[, c("geo_accession", "organism", "source")]

# User-confirmed: all 32 samples are Mus musculus; source fields identify
# 16 neuroretina and 16 RPE/choroid samples, consistent with the titles.

# Find the biological-information fields recorded across the samples.
amd_characteristic_nodes <- xml2::xml_find_all(
  amd_sample_nodes, ".//geo:Characteristics", ns = amd_ns
)

# List the characteristic field names and their number of occurrences.
amd_characteristic_tags <- xml2::xml_attr(amd_characteristic_nodes, "tag")
table(amd_characteristic_tags, useNA = "ifany")

# User-confirmed characteristic tags: age group, genotype/variation,
# strain background and tissue; 32 entries each.

# Read a named characteristic from each sample record.
amd_get_characteristic <- function(tag) {
  amd_path <- paste0(".//geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(amd_sample_nodes, amd_path, ns = amd_ns),
    trim = TRUE
  )
}

# Add the four recorded characteristics to the sample table.
amd_metadata_all$age_group <- amd_get_characteristic("age group")
amd_metadata_all$genotype <- amd_get_characteristic("genotype/variation")
amd_metadata_all$strain <- amd_get_characteristic("strain background")
amd_metadata_all$tissue <- amd_get_characteristic("tissue")

# Display the extracted values beside each GSM sample ID.
amd_metadata_all[, c("geo_accession", "age_group", "genotype", "strain", "tissue")]

# User-confirmed characteristics: all samples are C57BL/6J, with young/aged
# and WT/Cfh-/- labels; tissue fields identify 16 neuroretina and 16 RPE/choroid.

# Read the RNA-extraction description recorded for each sample.
amd_extraction_protocol <- xml2::xml_text(
  xml2::xml_find_first(
    amd_sample_nodes, ".//geo:Extract-Protocol", ns = amd_ns
  ),
  trim = TRUE
)
names(amd_extraction_protocol) <- amd_metadata_all$geo_accession

# Print each distinct description once to review tissue preparation details.
cat(unique(amd_extraction_protocol), sep = "\n\n")
cat("\n")

# User-reviewed extraction: TRIzol reagent followed by RNeasy mini spin
# columns and DNase I treatment. Sample pooling is not stated in this text.

# Check that metadata and expression columns contain the same unique GSM IDs.
amd_sample_ids <- colnames(amd_expr_all)
stopifnot(
  !anyNA(amd_metadata_all$geo_accession),
  !anyNA(amd_sample_ids),
  anyDuplicated(amd_metadata_all$geo_accession) == 0L,
  anyDuplicated(amd_sample_ids) == 0L,
  setequal(amd_metadata_all$geo_accession, amd_sample_ids)
)

# Arrange a metadata copy in expression-column order.
amd_metadata_matched <- amd_metadata_all[
  match(amd_sample_ids, amd_metadata_all$geo_accession),
  , drop = FALSE
]
rownames(amd_metadata_matched) <- NULL

# Confirm the metadata dimensions and exact sample-ID order.
dim(amd_metadata_matched)
identical(amd_metadata_matched$geo_accession, colnames(amd_expr_all))

# User-confirmed: matched metadata has 32 rows; exact sample-ID order is TRUE.

# Select the neuroretina metadata and its corresponding expression columns.
amd_metadata_selected <- amd_metadata_matched[
  amd_metadata_matched$tissue == "neuroretina",
  , drop = FALSE
]
rownames(amd_metadata_selected) <- NULL
amd_expr_selected <- amd_expr_all[
  , amd_metadata_selected$geo_accession, drop = FALSE
]

# Check the selected dimensions and exact sample-ID order.
dim(amd_metadata_selected)
dim(amd_expr_selected)
identical(amd_metadata_selected$geo_accession, colnames(amd_expr_selected))

# Count the selected samples in each age and genotype group.
with(amd_metadata_selected, table(age_group, genotype))

# User-confirmed: selected metadata is 16 x 8 and expression is 35,556 x 16;
# exact sample-ID order is TRUE; young/aged WT and Cfh-/- each have four samples.

# Paper reports RMA normalization and log2-scaled Figure 1 expression values.
# This does not by itself establish the scale of the deposited expression table.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC3465015/
# GEO guidance distinguishes logged and unlogged deposited microarray values.
# https://www.ncbi.nlm.nih.gov/geo/info/geo2r.html

# Inspect the selected expression scale before choosing a transformation.
amd_expression_quantiles <- quantile(
  amd_expr_selected, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1)
)
amd_expression_quantiles

# User-reviewed selected quantiles: 6, 34.5, 92.1, 298.2, 4315.91, 16719.9
# at 0%, 25%, 50%, 75%, 99% and 100%; consistent with unlogged intensities.

# Create a log2-scale copy for QC; all selected values are positive.
amd_expr_log2 <- log2(amd_expr_selected)

# Check the dimensions, finite values and transformed range.
dim(amd_expr_log2)
all(is.finite(amd_expr_log2))
range(amd_expr_log2)

# User-confirmed: log2 QC copy is 35,556 x 16; all values are finite.
# Observed log2 expression range: 2.584963 to 14.029279.

# Define a boxplot comparing the log2 expression distributions across samples.
amd_draw_boxplot <- function() {
  amd_previous_par <- par(mar = c(8, 4.5, 3, 1))
  on.exit(par(amd_previous_par))
  boxplot(
    amd_expr_log2,
    las = 2,
    cex.axis = 0.7,
    col = "lightblue",
    outline = FALSE,
    main = "CFH-deficiency model: neuroretina distributions",
    ylab = "Log2 expression value"
  )
}

# Display the plot in R; extreme-value dots are hidden for readability.
amd_draw_boxplot()

# Save the same plot as a high-resolution PNG in the AMD output folder.
amd_boxplot_file <- file.path(amd_output_dir, "AMD_GSE38671_log2_expression_boxplot.png")
png(amd_boxplot_file, width = 2400, height = 1800, res = 300)
amd_draw_boxplot()
dev.off()
file.exists(amd_boxplot_file)

# User-reviewed boxplot: medians are close (~6.5) and overall spreads are similar.
# GSM947360 has a slightly narrower box; no obvious global distribution outlier.
# Retain all 16 for further QC; this plot does not assess gene-specific differences.

# Estimate one smooth expression distribution per sample.
amd_density_curves <- lapply(
  seq_len(ncol(amd_expr_log2)),
  function(i) density(amd_expr_log2[, i])
)
names(amd_density_curves) <- colnames(amd_expr_log2)
amd_density_colors <- hcl.colors(ncol(amd_expr_log2), palette = "Dark 3")

# Define the figure once for consistent display and saving.
amd_draw_density <- function() {
  amd_xlim <- range(sapply(amd_density_curves, function(d) range(d$x)))
  amd_ylim <- c(0, max(sapply(amd_density_curves, function(d) max(d$y))))

  plot(
    amd_density_curves[[1]],
    xlim = amd_xlim, ylim = amd_ylim,
    col = amd_density_colors[1], lwd = 2,
    main = "CFH-deficiency model: neuroretina density",
    xlab = "Log2 expression value", ylab = "Density"
  )
  for (i in 2:length(amd_density_curves)) {
    lines(amd_density_curves[[i]], col = amd_density_colors[i], lwd = 2)
  }
  legend(
    "topright", legend = names(amd_density_curves),
    col = amd_density_colors, lty = 1, lwd = 2,
    bty = "n", cex = 0.65
  )
}

# Display the density plot in R.
amd_draw_density()

# Save the same figure in the AMD output folder.
amd_density_file <- file.path(amd_output_dir, "AMD_GSE38671_log2_expression_density.png")
png(amd_density_file, width = 2400, height = 1800, res = 300)
amd_draw_density()
dev.off()
file.exists(amd_density_file)

# User-reviewed density plot: curves broadly overlap, with similar peak locations
# and tails. GSM947360 has the tallest, more concentrated central peak.
# Retain all 16 for further QC; density differences alone do not justify exclusion.

# Retrieve the supplied feature IDs in expression-matrix row order.
amd_gene_ids <- rownames(amd_expr_log2)

# Inspect ten examples and check whether the IDs contain digits only.
head(amd_gene_ids, 10)
table(Digits_only = grepl("^[0-9]+$", amd_gene_ids))

# Confirm that matrix row IDs match the original ID_REF column as text.
identical(amd_gene_ids, as.character(amd_raw$ID_REF))

# User-confirmed: all 35,556 feature IDs contain digits only and match ID_REF
# in order (TRUE); examples begin 10338001 through 10338010.

# Extract the supplied platform table into an annotation output subfolder.
amd_annotation_dir <- file.path(amd_output_dir, "annotation")
dir.create(amd_annotation_dir, recursive = TRUE, showWarnings = FALSE)
untar(amd_metadata_archive, files = "GPL6246-tbl-1.txt", exdir = amd_annotation_dir)

# Confirm that R can find the extracted platform table.
amd_platform_annotation_file <- file.path(amd_annotation_dir, "GPL6246-tbl-1.txt")
file.exists(amd_platform_annotation_file)

# User-confirmed: the extracted GPL6246 platform table exists (TRUE).

# Preview the first lines to check for column headers before importing.
amd_annotation_preview <- readLines(
  amd_platform_annotation_file, n = 5, warn = FALSE
)
cat(amd_annotation_preview, sep = "\n\n")

# User-reviewed preview starts with data (10344614), not column headers.
# Annotation fields include long descriptions and multiple assignments.

# Read the headerless platform table, keeping all fields as text.
amd_platform_annotation <- read.delim(
  amd_platform_annotation_file,
  header = FALSE,
  colClasses = "character",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check the table size and preview only the first three columns.
dim(amd_platform_annotation)
head(amd_platform_annotation[, 1:3], 5)

# User-confirmed: platform annotation is 35,557 x 12; temporary names are V1...V12.
# First-row ID is 10344614. Expression rows remain 35,556; coverage is unchecked.

# Locate GPL6246 and its supplied annotation-column definitions.
amd_platform_nodes <- xml2::xml_find_all(
  amd_metadata_xml, "/geo:MINiML/geo:Platform[@iid='GPL6246']", ns = amd_ns
)
amd_platform_column_nodes <- xml2::xml_find_all(
  amd_platform_nodes, ".//geo:Data-Table/geo:Column", ns = amd_ns
)

# Display the supplied column names beside their table positions.
amd_platform_columns <- data.frame(
  position = xml2::xml_attr(amd_platform_column_nodes, "position"),
  name = xml2::xml_text(
    xml2::xml_find_first(amd_platform_column_nodes, "geo:Name", ns = amd_ns)
  ),
  stringsAsFactors = FALSE
)
amd_platform_columns

# User-reviewed definitions list positions 1–12 in order; ID is column 1,
# gene_assignment column 10 and mrna_assignment column 11.

# Apply the supplied XML column names to the platform table.
colnames(amd_platform_annotation) <- amd_platform_columns$name

# Preview feature IDs beside the supplied gene-assignment text.
head(amd_platform_annotation[, c("ID", "gene_assignment")], 5)

# User-reviewed named preview: 10344614 has several assignment records;
# 10344618 and 10344622 have --- in the gene_assignment field.
# Long annotation strings caused the console to print columns in separate blocks.

# The optional View call needs XQuartz in this R session; use console previews.
# View(amd_platform_annotation[, c("ID", "gene_assignment")])

# Read the supplied description before splitting the combined annotation field.
amd_assignment_description <- xml2::xml_text(
  xml2::xml_find_first(
    amd_platform_nodes,
    ".//geo:Column[geo:Name='gene_assignment']/geo:Description",
    ns = amd_ns
  ),
  trim = TRUE
)
cat(amd_assignment_description, sep = "\n")

# User-reviewed field description: gene information for assigned mRNAs (multipart).
# It does not spell out field order; use manufacturer layout and reviewed examples.

# Load the annotation-formatting helper.
source(file.path(amd_output_dir, "03_AMD_annotation_format.R"))

# Preview clearer columns using the five platform records already reviewed.
amd_annotation_preview_clean <- amd_format_annotation(
  head(amd_platform_annotation, 5)
)
amd_annotation_preview_clean[, c("ID", "reference_id", "gene_symbol", "entrez_id")]

# User-reviewed formatted preview: five platform records became eight rows;
# multiple assignments and unassigned features were retained.

# Format all platform records into clearer annotation columns.
amd_annotation_clean <- amd_format_annotation(amd_platform_annotation)

# Check the table size and confirm that every platform row was retained.
dim(amd_annotation_clean)
identical(
  sort(unique(amd_annotation_clean$platform_row)),
  seq_len(nrow(amd_platform_annotation))
)

# Save the cleaner annotation table in the AMD annotation folder.
amd_annotation_clean_file <- file.path(
  amd_annotation_dir, "AMD_GSE38671_platform_annotation_clean_2026-11-03.csv"
)
write.csv(amd_annotation_clean, amd_annotation_clean_file, row.names = FALSE, na = "")
file.exists(amd_annotation_clean_file)

# User-reported full-format error: feature 10412503 has a three-field assignment.
# The cleaner table was not created and the CSV existence check returned FALSE.

# Inspect the original annotation for the feature that stopped formatting.
amd_problem_assignment <- amd_platform_annotation$gene_assignment[
  amd_platform_annotation$ID == "10412503"
]
cat(amd_problem_assignment, sep = "\n\n")
cat("\n")

# User supplied the long annotation for 10412503, showing many assignments.
# The exact three-field entry must be isolated before deciding how to handle it.

# Separate the assignments and their individual fields.
amd_problem_records <- unlist(
  strsplit(amd_problem_assignment, "///", fixed = TRUE), use.names = FALSE
)
amd_problem_parts <- lapply(
  strsplit(amd_problem_records, "//", fixed = TRUE), trimws
)

# Display only the three-field entries reported by the formatter.
which(lengths(amd_problem_parts) == 3L)
amd_problem_parts[lengths(amd_problem_parts) == 3L]

# User-reviewed assignment 420: ENSMUST00000178728, Gm10340, and the literal
# [WARNING: THIS FIELD TRUNCATED] marker; remaining details are unavailable.

# Reload the formatter with explicit handling of the supplied truncation warning.
source(file.path(amd_output_dir, "03_AMD_annotation_format.R"))

# Rebuild the cleaner table and confirm that every platform row was retained.
amd_annotation_clean <- amd_format_annotation(amd_platform_annotation)
dim(amd_annotation_clean)
identical(
  sort(unique(amd_annotation_clean$platform_row)),
  seq_len(nrow(amd_platform_annotation))
)

# Preview the flagged assignment that previously stopped formatting.
amd_annotation_clean[
  amd_annotation_clean$ID == "10412503" & amd_annotation_clean$assignment_truncated,
  c("ID", "reference_id", "gene_symbol", "annotation_truncated")
]

# Save the cleaner table, including original text and truncation flags.
amd_annotation_clean_file <- file.path(
  amd_annotation_dir, "AMD_GSE38671_platform_annotation_clean_2026-11-03.csv"
)
write.csv(amd_annotation_clean, amd_annotation_clean_file, row.names = FALSE, na = "")
file.exists(amd_annotation_clean_file)

# User-confirmed: cleaner annotation is 298,709 x 12; all platform rows retained.
# The flagged 10412503 assignment is preserved and the CSV existence check is TRUE.

# Preview readable annotation columns directly in the R console.
head(
  amd_annotation_clean[, c("ID", "gene_symbol", "description", "entrez_id")],
  10
)

# User-reviewed full-table preview: 10344614 has four distinct supplied labels;
# 10344624 repeats Lypla1/18777. Assignment rows are not independent measurements.

# Check missing/blank IDs and duplicates in the original, unexpanded platform table.
amd_platform_ids <- amd_platform_annotation$ID
sum(is.na(amd_platform_ids) | trimws(amd_platform_ids) == "")
anyDuplicated(amd_platform_ids)

# Count expression IDs without a platform record and unused platform IDs.
amd_ids_without_annotation <- setdiff(amd_gene_ids, amd_platform_ids)
amd_extra_platform_ids <- setdiff(amd_platform_ids, amd_gene_ids)
c(
  missing_from_platform = length(amd_ids_without_annotation),
  extra_in_platform = length(amd_extra_platform_ids)
)

# Display up to ten examples from each group.
head(amd_ids_without_annotation, 10)
head(amd_extra_platform_ids, 10)

# User-confirmed: no missing/blank or duplicate platform IDs; all 35,556
# expression IDs are covered. Extra platform-only ID: 10338063.

# Match one platform record to each expression feature in matrix order.
amd_feature_annotation <- amd_platform_annotation[
  match(rownames(amd_expr_log2), amd_platform_annotation$ID),
  , drop = FALSE
]
rownames(amd_feature_annotation) <- NULL

# Check annotation dimensions and exact feature-ID order.
dim(amd_feature_annotation)
identical(amd_feature_annotation$ID, rownames(amd_expr_log2))

# User-confirmed: feature annotation is 35,556 x 12 and exact feature-ID order TRUE.
# Follow the existing DR/RP/LHON pairwise top-500 MDS setting on the log2 QC matrix.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

# Load the familiar projection style and label the age/genotype groups.
source(file.path(amd_output_dir, "04_AMD_projection_style.R"))
amd_metadata_selected$group <- factor(
  paste(amd_metadata_selected$age_group, amd_metadata_selected$genotype),
  levels = c("young WT", "young Cfh-/-", "aged WT", "aged Cfh-/-")
)
stopifnot(
  identical(amd_metadata_selected$geo_accession, colnames(amd_expr_log2)),
  !anyNA(amd_metadata_selected$group)
)

# Calculate MDS using the 500 most different features for each sample pair.
amd_mds <- limma::plotMDS(
  amd_expr_log2, top = 500, gene.selection = "pairwise", plot = FALSE
)

# Display the sample map with group colours and GSM labels.
amd_mds_plot <- amd_projection_plot(
  amd_mds, amd_metadata_selected,
  title = "CFH-deficiency model: neuroretina MDS",
  axis_names = c("MDS1", "MDS2")
)
print(amd_mds_plot)

# Save the figure in the AMD output folder.
amd_mds_file <- file.path(amd_output_dir, "AMD_GSE38671_MDS.png")
amd_save_projection(amd_mds_plot, amd_mds_file)
file.exists(amd_mds_file)

# Calculate PCA using one shared set of the 500 most variable features.
amd_pca <- limma::plotMDS(
  amd_expr_log2, top = 500, gene.selection = "common", plot = FALSE
)

# Display PCA with the same group colours and GSM labels as MDS.
amd_pca_plot <- amd_projection_plot(
  amd_pca, amd_metadata_selected,
  title = "CFH-deficiency model: neuroretina PCA",
  axis_names = c("PC1", "PC2")
)
print(amd_pca_plot)

# Save the figure in the AMD output folder.
amd_pca_file <- file.path(amd_output_dir, "AMD_GSE38671_PCA_top500.png")
amd_save_projection(amd_pca_plot, amd_pca_file)
file.exists(amd_pca_file)

# User supplied both figures: MDS axes 24% and 12%; PCA axes 33% and 13%.
# GSM947360 (young WT) is visually separated in both views and needs QC review.
# The other three young WT samples are closer; the four groups partly overlap.
# Aged Cfh-/- samples are dispersed; aged WT GSM947365 is apart from its other three.
# Earlier GSM947360 distribution differences support reviewing it, not a proven cause.
# Retain all 16 for now; flag GSM947360 before a final sample-inclusion decision.
# These views do not establish differential expression or resemblance to spaceflight.

# Save a dated script copy, AMD checkpoint, package versions and console history.
source(file.path(amd_output_dir, "05_AMD_save_session.R"))
