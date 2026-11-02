# Completed user-run precheck: 2 November 2026; see docs/autoimmune-uveitis-precheck.md.
# Preserve original mixed labels; the separate annotation export records lookup gaps.
# Live Ensembl lookups may change; the dated exported table preserves this run.

# Capture the repository root before changing R's working directory.
precheck_project_dir <- normalizePath(getwd(), mustWork = TRUE)

# Use existing downloads or the default input folder; keep new results local.
uveitis_dir <- Sys.getenv(
  "UVEITIS_DATA_DIR",
  unset = file.path(precheck_project_dir, "data", "raw", "Autoimmune_Uveitis_GSE284505")
)
uveitis_output_dir <- file.path(precheck_project_dir, "results", "autoimmune_uveitis")
dir.create(uveitis_output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the output folder as R's default save location.
setwd(uveitis_output_dir)
getwd()

# List the downloaded files before choosing what to open together.
list.files(uveitis_dir)

# Import the tab-separated count table, preserving its original column names.
uveitis_count_file <- file.path(uveitis_dir, "GSE284505_R.HTSeq.counts.txt")
uveitis_raw <- read.delim(
  uveitis_count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect the table size, column names and first five rows before further checks.
dim(uveitis_raw)
colnames(uveitis_raw)
head(uveitis_raw, 5)

# User-confirmed: 56,859 rows, AccID and 12 sample columns.
# CON, EAU, CSD-EAU and CSD-CON each have three labelled columns.
# Biological group meanings and sample independence await metadata review.

# Check for missing entries anywhere in the imported table.
anyNA(uveitis_raw)

# Check for missing or blank gene labels.
sum(is.na(uveitis_raw$AccID) | trimws(uveitis_raw$AccID) == "")

# Check whether any gene label appears more than once.
anyDuplicated(uveitis_raw$AccID)

# Check that all sample columns are stored as numeric values.
sapply(uveitis_raw[, -1, drop = FALSE], is.numeric)

# User-confirmed: no missing entries, blank/repeated labels; all samples numeric.

# Build the numeric count matrix, keeping gene labels as row names.
uveitis_counts_all <- as.matrix(uveitis_raw[, -1, drop = FALSE])
rownames(uveitis_counts_all) <- uveitis_raw$AccID

# Confirm the number of gene rows and sample columns.
dim(uveitis_counts_all)

# Check for finite values and inspect the smallest and largest counts.
all(is.finite(uveitis_counts_all))
range(uveitis_counts_all)

# Check that counts are whole numbers without rounding or changing them.
all(uveitis_counts_all == round(uveitis_counts_all))

# User-confirmed: 56,859 x 12; all finite, nonnegative and whole-valued.
# Observed count range: 0 to 616,271.

# Locate the downloaded metadata archive.
uveitis_metadata_archive <- file.path(uveitis_dir, "GSE284505_family")

# List its contents before extracting or reading the metadata.
untar(uveitis_metadata_archive, list = TRUE)

# User-confirmed archive member: GSE284505_family.xml.

# Create a subfolder for extracted metadata within the output folder.
uveitis_metadata_dir <- file.path(uveitis_output_dir, "metadata")
dir.create(uveitis_metadata_dir, recursive = TRUE, showWarnings = FALSE)

# Extract the XML metadata file into that subfolder.
untar(
  uveitis_metadata_archive,
  files = "GSE284505_family.xml",
  exdir = uveitis_metadata_dir
)

# Confirm that R can find the extracted file.
uveitis_metadata_file <- file.path(uveitis_metadata_dir, "GSE284505_family.xml")
file.exists(uveitis_metadata_file)

# User-confirmed: the extracted XML metadata file exists.

# Read the XML metadata and display its structure before selecting sample records.
uveitis_metadata_xml <- xml2::read_xml(uveitis_metadata_file)
uveitis_metadata_xml

# User-confirmed: GEO MINiML XML overview shows 12 Sample records.

# Use the XML namespace shown in the overview to locate sample records.
uveitis_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")
uveitis_sample_nodes <- xml2::xml_find_all(
  uveitis_metadata_xml, "/geo:MINiML/geo:Sample", ns = uveitis_ns
)
length(uveitis_sample_nodes)

# Create one row per sample with its GEO identifier and descriptive title.
uveitis_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(uveitis_sample_nodes, "iid"),
  title = xml2::xml_text(
    xml2::xml_find_first(uveitis_sample_nodes, "geo:Title", ns = uveitis_ns)
  ),
  stringsAsFactors = FALSE
)
uveitis_metadata_all

# User-confirmed: 12 rows, with three titles for each of four conditions:
# healthy control, chronic sleep deprivation, IRBP/CFA and combined treatment.
# Metadata/count matching and independent replication remain to be checked.

# Find the biological-information fields recorded in the sample metadata.
uveitis_characteristic_nodes <- xml2::xml_find_all(
  uveitis_sample_nodes, ".//geo:Characteristics", ns = uveitis_ns
)
uveitis_characteristic_tags <- xml2::xml_attr(uveitis_characteristic_nodes, "tag")

# Show the available field names and how often each appears.
table(uveitis_characteristic_tags, useNA = "ifany")

# User-confirmed fields: genotype, tissue and treatment, each appearing 12 times.
# Sex and collection time are not separate fields in this inventory.

# Read growth, treatment and extraction descriptions for sex and collection timing.
uveitis_protocol_nodes <- xml2::xml_find_all(
  uveitis_sample_nodes,
  ".//geo:Growth-Protocol | .//geo:Treatment-Protocol | .//geo:Extract-Protocol",
  ns = uveitis_ns
)
uveitis_protocol_text <- xml2::xml_text(uveitis_protocol_nodes, trim = TRUE)

# Print each distinct description once so the full text is easy to review.
cat(unique(uveitis_protocol_text), sep = "\n\n")
cat("\n")

# Reviewed protocols describe RNA processing and library preparation.
# Sex and tissue-collection timing are unreported in the inspected fields.
# Verification against the full paper remains pending.

# Read one recorded characteristic per sample, following the RP/LHON pattern.
uveitis_get_characteristic <- function(tag) {
  uveitis_path <- paste0(".//geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(uveitis_sample_nodes, uveitis_path, ns = uveitis_ns),
    trim = TRUE
  )
}

# Add the available biological information while preserving metadata row order.
uveitis_metadata_all$genotype <- uveitis_get_characteristic("genotype")
uveitis_metadata_all$tissue <- uveitis_get_characteristic("tissue")
uveitis_metadata_all$treatment <- uveitis_get_characteristic("treatment")
uveitis_metadata_all

# User-confirmed: all 12 samples are labelled Retina.
# The genotype field contains condition codes, with three samples per condition.

# Build sample labels from the recorded condition code and title's rep number.
uveitis_metadata_all$sample_id <- paste0(
  uveitis_metadata_all$genotype,
  sub(".*rep", "", uveitis_metadata_all$title)
)

# Display the link between GEO accessions, sample labels and condition codes.
uveitis_metadata_all[, c("geo_accession", "sample_id", "genotype")]

# Check that labels are unique and cover the same samples as the count columns.
anyDuplicated(uveitis_metadata_all$sample_id)
setequal(uveitis_metadata_all$sample_id, colnames(uveitis_counts_all))

# User-confirmed: sample labels are unique and the sample-name sets match.

# Arrange a metadata copy in the same order as the count-matrix columns.
uveitis_match <- match(colnames(uveitis_counts_all), uveitis_metadata_all$sample_id)
uveitis_metadata_matched <- uveitis_metadata_all[uveitis_match, , drop = FALSE]
rownames(uveitis_metadata_matched) <- NULL

# Check the metadata size and exact agreement of sample names and order.
dim(uveitis_metadata_matched)
identical(colnames(uveitis_counts_all), uveitis_metadata_matched$sample_id)

# Display the sample-to-metadata links in their aligned order.
uveitis_metadata_matched[, c("sample_id", "geo_accession", "genotype")]

# User-confirmed: matched metadata is 12 x 6 and exact sample-order check is TRUE.

# Select EAU and healthy-control metadata for the autoimmune-uveitis comparison.
uveitis_metadata_selected <- uveitis_metadata_matched[
  uveitis_metadata_matched$genotype %in% c("CON", "EAU"),
  , drop = FALSE
]
rownames(uveitis_metadata_selected) <- NULL

# Check the subset size and number of measured samples per condition.
dim(uveitis_metadata_selected)
table(uveitis_metadata_selected$genotype)

# Display the sample identities retained for this comparison.
uveitis_metadata_selected[, c("sample_id", "geo_accession", "genotype")]

# User-confirmed: selected metadata is 6 x 6, with three CON and three EAU samples.
# The user agreed to retain EAU versus CON for the primary comparison.

# Select count columns using the selected metadata's sample labels.
uveitis_counts_selected <- uveitis_counts_all[
  , uveitis_metadata_selected$sample_id, drop = FALSE
]

# Confirm count dimensions and exact sample alignment with selected metadata.
dim(uveitis_counts_selected)
identical(colnames(uveitis_counts_selected), uveitis_metadata_selected$sample_id)

# Preview the first five gene rows across the selected samples.
head(uveitis_counts_selected, 5)

# User-confirmed: selected counts are 56,859 x 6 and metadata alignment is TRUE.

# Sum gene counts within each selected sample.
uveitis_total_counts <- colSums(uveitis_counts_selected)
uveitis_total_counts

# Display totals from smallest to largest for comparison.
sort(uveitis_total_counts)

# User-reported totals span 11,745,016 (CON1) to 14,955,036 (EAU2).

# Count genes with a positive count in each selected sample.
uveitis_detected_genes <- colSums(uveitis_counts_selected > 0)
uveitis_detected_genes

# Display detected-gene totals from smallest to largest for comparison.
sort(uveitis_detected_genes)

# User-reported detected-gene totals span 27,183 (CON1) to 28,509 (EAU2).

# Calculate all-gene log-CPM for the initial distribution plot.
# This uses sample count totals; expression filtering and TMM are later steps.
uveitis_log_cpm <- edgeR::cpm(
  uveitis_counts_selected,
  log = TRUE,
  prior.count = 2
)

# Define the plot once for consistent on-screen display and PNG export.
uveitis_draw_boxplot <- function() {
  boxplot(
    uveitis_log_cpm,
    las = 2,
    col = "lightblue",
    outline = FALSE,
    main = "Autoimmune uveitis: expression distributions",
    xlab = "Sample",
    ylab = "Log2 counts per million (CPM)"
  )
}

# Display the boxplot in R; outline = FALSE hides dots without removing genes.
uveitis_draw_boxplot()

# Save the same plot in the Autoimmune Uveitis output folder.
uveitis_boxplot_file <- file.path(
  uveitis_output_dir,
  "Autoimmune_Uveitis_GSE284505_logCPM_boxplot.png"
)
png(uveitis_boxplot_file, width = 2400, height = 1800, res = 300)
uveitis_draw_boxplot()
dev.off()
file.exists(uveitis_boxplot_file)

# User-reviewed boxplot: distributions broadly overlap; retain all six samples.
# EAU2's median is above the zero-count floor, consistent with >50% positive rows.

# Estimate one density curve per sample using the existing all-gene log-CPM.
uveitis_density_curves <- lapply(
  seq_len(ncol(uveitis_log_cpm)),
  function(i) density(uveitis_log_cpm[, i])
)
names(uveitis_density_curves) <- colnames(uveitis_log_cpm)
uveitis_density_colors <- hcl.colors(ncol(uveitis_log_cpm), palette = "Dark 3")

# Define the plot once, using axis limits that include every sample curve.
uveitis_draw_density <- function() {
  uveitis_xlim <- range(sapply(uveitis_density_curves, function(d) range(d$x)))
  uveitis_ylim <- c(0, max(sapply(uveitis_density_curves, function(d) max(d$y))))
  plot(
    uveitis_density_curves[[1]],
    xlim = uveitis_xlim, ylim = uveitis_ylim,
    col = uveitis_density_colors[1], lwd = 2,
    main = "Autoimmune uveitis: expression density",
    xlab = "Log2 counts per million (CPM)", ylab = "Density"
  )
  for (i in 2:length(uveitis_density_curves)) {
    lines(uveitis_density_curves[[i]], col = uveitis_density_colors[i], lwd = 2)
  }
  legend(
    "topright", legend = names(uveitis_density_curves),
    col = uveitis_density_colors, lty = 1, lwd = 2, bty = "n", cex = 0.75
  )
}

# Display the density plot in R.
uveitis_draw_density()

# Save the same plot, close the PNG device and confirm the file exists.
uveitis_density_file <- file.path(
  uveitis_output_dir,
  "Autoimmune_Uveitis_GSE284505_logCPM_density.png"
)
png(uveitis_density_file, width = 2400, height = 1800, res = 300)
uveitis_draw_density()
dev.off()
file.exists(uveitis_density_file)

# User-reviewed density curves broadly overlap, with a large zero/low-count peak.
# Continue all six samples; these initial plots were generated before TMM.

# Retrieve gene labels in the selected count matrix's row order.
uveitis_gene_labels <- rownames(uveitis_counts_selected)
head(uveitis_gene_labels, n = 10)

# Check which labels start with the mouse Ensembl gene-ID prefix.
uveitis_ensembl_prefix <- grepl("^ENSMUSG", uveitis_gene_labels)
table(Ensembl_prefix = uveitis_ensembl_prefix)

# Show examples that use another label format.
head(uveitis_gene_labels[!uveitis_ensembl_prefix], n = 10)

# Confirm the gene labels still match the original AccID column in order.
identical(uveitis_gene_labels, uveitis_raw$AccID)

# User-confirmed: 1,548 labels have the ENSMUSG prefix; 55,311 do not.
# Gene-label order matches the original AccID column exactly.

# Preserve the supplied mixed-format labels under their original field name.
uveitis_gene_annotation <- data.frame(
  AccID = uveitis_gene_labels,
  stringsAsFactors = FALSE
)

# Check the annotation size and preview the supplied labels.
dim(uveitis_gene_annotation)
head(uveitis_gene_annotation, n = 10)

# Inspect examples from the Ensembl-prefixed subset found in this dataset.
head(uveitis_gene_labels[uveitis_ensembl_prefix], n = 10)

# Confirm annotations follow the count matrix's gene-row order exactly.
identical(uveitis_gene_annotation$AccID, rownames(uveitis_counts_selected))

# User-confirmed: annotation is 56,859 x 1 and gene-row alignment is TRUE.
# Preserve full identifiers, including ENSMUSG; common-ID mapping is a later task.

# Check the platform section for a supplied gene-annotation table.
xml2::xml_find_all(
  uveitis_metadata_xml,
  "/geo:MINiML/geo:Platform/geo:Data-Table/*",
  ns = uveitis_ns
)

# Inspect processing notes for the reference genome or annotation used for gene IDs.
uveitis_processing_text <- xml2::xml_text(
  xml2::xml_find_all(
    uveitis_sample_nodes, ".//geo:Data-Processing", ns = uveitis_ns
  ),
  trim = TRUE
)

# Display each distinct description once before choosing an annotation lookup.
cat(unique(uveitis_processing_text), sep = "\n\n")
cat("\n")

# User-reviewed processing notes report GRCm38; no annotation release is stated.
# User-reported platform table query returned zero nodes: no table content there.

# BioMart was unavailable on 2 November; its failed connections created no objects.

# Read the annotation release from Ensembl's responding API.
uveitis_ensembl_release <- jsonlite::fromJSON(
  "https://rest.ensembl.org/info/data?content-type=application/json"
)
uveitis_ensembl_release$releases

# Select ten supplied Ensembl IDs for the first annotation lookup.
uveitis_lookup_ids <- head(uveitis_gene_labels[uveitis_ensembl_prefix], n = 10)

# Request their gene records from Ensembl and check that the request succeeded.
uveitis_lookup_response <- httr::POST(
  "https://rest.ensembl.org/lookup/id",
  httr::accept_json(),
  body = list(ids = unname(uveitis_lookup_ids)),
  encode = "json",
  httr::timeout(30)
)
httr::stop_for_status(uveitis_lookup_response)
uveitis_lookup_result <- httr::content(
  uveitis_lookup_response, as = "parsed",
  type = "application/json", encoding = "UTF-8"
)

# Display the first two records before building the annotation table.
uveitis_lookup_result[uveitis_lookup_ids[1:2]]

# User-reviewed: first two mouse gene records are Gm57191 (lncRNA) and
# Gm56487 (unprocessed_pseudogene), annotated on GRCm39 in release 116.

# Extract a field for each requested ID, preserving absent values as NA.
uveitis_lookup_field <- function(field) {
  vapply(uveitis_lookup_ids, function(id) {
    value <- uveitis_lookup_result[[id]][[field]]
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1), USE.NAMES = FALSE)
}

# Build a separate lookup table, keeping source IDs and annotation provenance.
uveitis_lookup_preview <- data.frame(
  AccID = uveitis_lookup_ids,
  ensembl_id = uveitis_lookup_field("id"),
  gene_name = uveitis_lookup_field("display_name"),
  description = uveitis_lookup_field("description"),
  biotype = uveitis_lookup_field("biotype"),
  species = uveitis_lookup_field("species"),
  annotation_assembly = uveitis_lookup_field("assembly_name"),
  annotation_release = uveitis_ensembl_release$releases[1],
  stringsAsFactors = FALSE
)

# Review all ten mappings and count IDs with a returned Ensembl record.
uveitis_lookup_preview[, c("AccID", "gene_name", "biotype", "species")]
table(Found_in_Ensembl = !is.na(uveitis_lookup_preview$ensembl_id))

# User-confirmed: all ten IDs were found in mouse; two lacked display_name.

# Split the remaining Ensembl IDs into batches of at most 1,000.
uveitis_ensembl_ids <- uveitis_gene_labels[uveitis_ensembl_prefix]
uveitis_remaining_ids <- setdiff(uveitis_ensembl_ids, uveitis_lookup_ids)
uveitis_lookup_batches <- split(
  uveitis_remaining_ids,
  ceiling(seq_along(uveitis_remaining_ids) / 1000)
)
uveitis_all_lookup_result <- uveitis_lookup_result

# Retrieve each batch and retain the records alongside the first ten results.
for (i in seq_along(uveitis_lookup_batches)) {
  uveitis_lookup_response <- httr::POST(
    "https://rest.ensembl.org/lookup/id",
    httr::accept_json(),
    body = list(ids = unname(uveitis_lookup_batches[[i]])),
    encode = "json",
    httr::timeout(60)
  )
  httr::stop_for_status(uveitis_lookup_response)
  uveitis_batch_result <- httr::content(
    uveitis_lookup_response, as = "parsed",
    type = "application/json", encoding = "UTF-8"
  )
  uveitis_all_lookup_result <- c(uveitis_all_lookup_result, uveitis_batch_result)
  cat("Completed batch", i, "of", length(uveitis_lookup_batches), "\n")
}

# Count returned records across all 1,548 requested Ensembl IDs.
uveitis_all_lookup_found <- vapply(
  uveitis_ensembl_ids,
  function(id) !is.null(uveitis_all_lookup_result[[id]]),
  logical(1)
)
table(Found_in_Ensembl = uveitis_all_lookup_found)

# User-confirmed: both batches completed; 1,464 records found and 84 not found.

# Extract fields for all requested IDs, preserving missing records/fields as NA.
uveitis_all_lookup_field <- function(field) {
  vapply(uveitis_ensembl_ids, function(id) {
    value <- uveitis_all_lookup_result[[id]][[field]]
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1), USE.NAMES = FALSE)
}

# Store the Ensembl annotations separately, retaining every original requested ID.
uveitis_ensembl_annotation <- data.frame(
  AccID = uveitis_ensembl_ids,
  ensembl_id = uveitis_all_lookup_field("id"),
  gene_name = uveitis_all_lookup_field("display_name"),
  description = uveitis_all_lookup_field("description"),
  biotype = uveitis_all_lookup_field("biotype"),
  species = uveitis_all_lookup_field("species"),
  annotation_assembly = uveitis_all_lookup_field("assembly_name"),
  annotation_release = uveitis_ensembl_release$releases[1],
  stringsAsFactors = FALSE
)

# Check table size, distinguish returned records from available names, and review species.
dim(uveitis_ensembl_annotation)
table(
  Record_found = !is.na(uveitis_ensembl_annotation$ensembl_id),
  Name_available = !is.na(uveitis_ensembl_annotation$gene_name) &
    trimws(uveitis_ensembl_annotation$gene_name) != ""
)
table(uveitis_ensembl_annotation$species, useNA = "ifany")

# User-confirmed: 1,548 x 8; 1,440 named, 24 found without names, 84 not found.
# All 1,464 returned gene records are mus_musculus.

# Match Ensembl annotations to the full gene table using the original AccID labels.
uveitis_annotation_match <- match(
  uveitis_gene_labels, uveitis_ensembl_annotation$AccID
)
uveitis_gene_annotation <- data.frame(
  AccID = uveitis_gene_labels,
  uveitis_ensembl_annotation[uveitis_annotation_match, -1, drop = FALSE],
  row.names = NULL,
  stringsAsFactors = FALSE
)

# Distinguish unqueried labels, unmatched IDs and records with or without names.
uveitis_gene_annotation$lookup_status <- "not_queried"
uveitis_gene_annotation$lookup_status[uveitis_ensembl_prefix] <- "not_found"
uveitis_gene_annotation$lookup_status[
  !is.na(uveitis_gene_annotation$ensembl_id)
] <- "found_without_name"
uveitis_gene_annotation$lookup_status[
  !is.na(uveitis_gene_annotation$gene_name) &
    trimws(uveitis_gene_annotation$gene_name) != ""
] <- "found_with_name"

# Confirm the full table size, gene-row alignment and annotation-status totals.
dim(uveitis_gene_annotation)
identical(uveitis_gene_annotation$AccID, rownames(uveitis_counts_selected))
table(uveitis_gene_annotation$lookup_status)

# User-confirmed: full annotation is 56,859 x 9 and gene-row alignment is TRUE.
# Status totals: 1,440 named, 24 unnamed, 84 not found and 55,311 not queried.

# Confirm sample and gene annotations match the selected count matrix.
stopifnot(
  identical(uveitis_metadata_selected$sample_id, colnames(uveitis_counts_selected)),
  identical(uveitis_gene_annotation$AccID, rownames(uveitis_counts_selected))
)

# Define the two experimental groups for the QC workflow.
uveitis_qc_group <- factor(
  uveitis_metadata_selected$genotype,
  levels = c("CON", "EAU")
)
stopifnot(!anyNA(uveitis_qc_group))
table(uveitis_qc_group)

# Assemble raw counts, sample metadata and gene annotations in one edgeR object.
uveitis_dge <- edgeR::DGEList(
  counts = uveitis_counts_selected,
  samples = uveitis_metadata_selected,
  group = uveitis_qc_group,
  genes = uveitis_gene_annotation
)

# Check dimensions, sample information and both alignments in the assembled object.
dim(uveitis_dge)
uveitis_dge$samples[, c("sample_id", "group", "lib.size", "norm.factors")]
c(
  sample_order = identical(uveitis_dge$samples$sample_id, colnames(uveitis_dge)),
  gene_order = identical(uveitis_dge$genes$AccID, rownames(uveitis_dge))
)

# User-confirmed: 56,859 x 6, three samples per group, factors 1, both orders TRUE.

# Apply the same expression-based count settings used for RP and LHON.
uveitis_keep_genes <- edgeR::filterByExpr(
  uveitis_dge,
  group = uveitis_qc_group,
  min.count = 10,
  min.total.count = 15
)
table(Retain_for_QC = uveitis_keep_genes)

# Create a QC copy and recalculate library totals over the retained gene rows.
uveitis_dge_qc <- uveitis_dge[
  uveitis_keep_genes,
  , keep.lib.sizes = FALSE
]

# Check retained dimensions and sample/gene alignment after filtering.
dim(uveitis_dge_qc)
c(
  sample_order = identical(
    uveitis_dge_qc$samples$sample_id, uveitis_metadata_selected$sample_id
  ),
  gene_order = identical(
    uveitis_dge_qc$genes$AccID, rownames(uveitis_dge_qc)
  )
)

# User-confirmed: 18,187 retained, 38,672 excluded; six samples and both orders TRUE.

# Estimate TMM library-scaling factors from the filtered counts.
uveitis_dge_qc <- edgeR::normLibSizes(uveitis_dge_qc, method = "TMM")

# Display retained library totals and the additional TMM scaling factors.
uveitis_dge_qc$samples[, c("sample_id", "group", "lib.size", "norm.factors")]

# Check that every estimated factor is finite and positive.
c(
  finite_factors = all(is.finite(uveitis_dge_qc$samples$norm.factors)),
  positive_factors = all(uveitis_dge_qc$samples$norm.factors > 0)
)

# User-confirmed: factors are finite/positive, ranging from 0.9857363 to 1.0340989.

# Calculate filtered log-CPM using TMM-adjusted library sizes.
uveitis_log_cpm_qc <- edgeR::cpm(
  uveitis_dge_qc,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

# Check the matrix size, finite values and sample order.
dim(uveitis_log_cpm_qc)
c(
  finite_logCPM = all(is.finite(uveitis_log_cpm_qc)),
  sample_order = identical(
    colnames(uveitis_log_cpm_qc), uveitis_dge_qc$samples$sample_id
  )
)

# User-confirmed: 18,187 x 6; finite_logCPM and sample_order are both TRUE.

# Load the familiar point/label style with an outside condition legend.
source(file.path(precheck_project_dir, "scripts", "14_Uveitis_projection_style.R"))

# Calculate pairwise top-500 MDS and common top-500 PCA from the same matrix.
uveitis_mds <- limma::plotMDS(
  uveitis_log_cpm_qc, top = 500, gene.selection = "pairwise", plot = FALSE
)
uveitis_pca <- limma::plotMDS(
  uveitis_log_cpm_qc, top = 500, gene.selection = "common", plot = FALSE
)

# Format both projections using their calculated coordinates and axis percentages.
uveitis_mds_plot <- uveitis_projection_plot(
  uveitis_mds, uveitis_dge_qc$samples,
  title = "Autoimmune uveitis: MDS after TMM",
  axis_names = c("Leading logFC dim 1", "Leading logFC dim 2")
)
uveitis_pca_plot <- uveitis_projection_plot(
  uveitis_pca, uveitis_dge_qc$samples,
  title = "Autoimmune uveitis: PCA after TMM",
  axis_names = c("PC1", "PC2")
)

# Display both figures and keep simple redraw commands.
uveitis_draw_mds <- function() print(uveitis_mds_plot)
uveitis_draw_pca <- function() print(uveitis_pca_plot)
uveitis_draw_mds()
uveitis_draw_pca()

# Set filenames in the current output folder.
uveitis_mds_file <- file.path(uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_TMM_MDS.png")
uveitis_pca_file <- file.path(uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_TMM_PCA_top500.png")

# Save each figure and close its PNG device.
png(uveitis_mds_file, width = 2400, height = 1800, res = 300)
uveitis_draw_mds()
dev.off()
png(uveitis_pca_file, width = 2400, height = 1800, res = 300)
uveitis_draw_pca()
dev.off()

# Confirm both files exist before reviewing the plots together.
c(MDS = file.exists(uveitis_mds_file), PCA = file.exists(uveitis_pca_file))

# User-reviewed: MDS axes 78%/7%; PCA axes 91%/4%; both separate CON from EAU.
# CON3 differs on the second axis; retain all six samples at this stage.

# Save the current reviewed stage to a new 2 November checkpoint and tables.
source(file.path(precheck_project_dir, "scripts", "15_Autoimmune_Uveitis_GSE284505_save.R"))
