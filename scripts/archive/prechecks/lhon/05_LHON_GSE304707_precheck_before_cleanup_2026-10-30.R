# LHON-linked mitochondrial optic neuropathy model: GSE304707
# Manual processed-data precheck, started 29 October 2026 (Asia/Manila).
# Run one block at a time and review its output before proceeding.
# Use lhon_ names so previously saved RP/glaucoma objects are preserved.
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304707
# Local filesystem checks confirmed the expression file is uncompressed CSV.
# GSE304707_family is a tar archive, not a plain SOFT text file.
# No expression import or analysis has been executed by the assistant.

# STEP 0: Set the destination for all LHON precheck outputs.
# The existing folder is listed on disk as lowercase "lhon".
# No packages or previous dataset checkpoints are required for Steps 0-1.
# Every later figure, table and checkpoint save must use lhon_output_dir.
# Setting this variable alone does not write output or change getwd().
# Store the folder path so later save commands can reuse it.
lhon_output_dir <- "/Volumes/MaxineDrive/INTERNSHIP/prequalitychecks/lhon"

# Check that R can find the output folder.
dir.exists(lhon_output_dir)
# Expected: TRUE. Review this before proceeding to Step 1.

# STEP 1: Import the deposited table and inspect its structure.
# Store the location of the downloaded dataset files.
lhon_dir <- "/Volumes/MaxineDrive/INTERNSHIP/datasets/LHON_GSE304707"

# Combine the folder and filename to identify the count file.
lhon_count_file <- file.path(lhon_dir, "GSE304707_gene_count_matrix.csv")

# Read the CSV into a table named lhon_raw.
# check.names = FALSE keeps the original column labels, including hyphens.
# stringsAsFactors = FALSE keeps text as text rather than categorical factors.
lhon_raw <- read.csv(
  lhon_count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check the number of rows and columns, in that order.
dim(lhon_raw)

# Display the column labels to distinguish gene information from samples.
colnames(lhon_raw)

# Preview the first five rows to see how the table is arranged.
# This preview does not establish that every row contains valid counts.
head(lhon_raw, 5)

# check.names = FALSE preserves the deposited column labels.
# dim() reports rows and columns; columns may include gene identifiers or
# annotation, so the total column count is not automatically sample count.
# Keep the original table intact. Identify gene/sample columns and inspect
# metadata before creating a numeric count matrix or selecting groups.
# Confirmed by user screenshot on 29 October 2026:
# lhon_raw has 49,714 rows and 13 columns: gene_name plus 12 sample columns.
# Sample order: WT-1, WT-2, WT-3, Y-1, Y-2, Y-3,
#               C-1, C-2, C-3, YC-1, YC-2, YC-3.
# The displayed first five gene labels are Gm10568, Gm7357, Gm6085,
# Sox17 and Gm37587. Full gene-ID uniqueness/count validity is still unchecked.
#
# Metadata interpretation verified by reading the downloaded family XML
# directly inside its tar archive (no extraction or expression analysis):
# WT = wild type; Y = nuclear YARS2 G186V mutant;
# C = mitochondrial COI V421A mutant; YC = both mutations.
# GSM9155578-580 map to WT-1..3; GSM9155581-583 to Y-1..3;
# GSM9155584-586 to C-1..3; GSM9155587-589 to YC-1..3.
# All 12 records list retina, male, C57BL/6J and age 12 months.
# Each group has three labelled replicate libraries. Animal independence,
# pooling and exact allelic state are not established by these labels alone.
# Keep all four genotype categories distinct in the sample metadata.
# Pending: manual metadata import/matching and subsequent count checks.

# TEACHING RESET, requested 29 October 2026:
# The source details above were read ahead by the assistant. They are provenance,
# not metadata checks completed by the user. From here, inspect files together
# one user-run block at a time and review output before advancing.

# STEP 2A: List what is inside the downloaded metadata archive.
# Metadata describes the samples, such as their experimental groups.
# First inspect the filenames inside the archive to decide how to open them.

# Identify the downloaded metadata archive using its folder and filename.
# This stores a path; it does not read the metadata yet.
lhon_metadata_archive <- file.path(lhon_dir, "GSE304707_family")

# List the filenames inside the archive.
# list = TRUE displays the contents without extracting files or loading them.
untar(lhon_metadata_archive, list = TRUE)

# Review the returned filenames together before extracting or importing them.
# Confirmed by user output on 29 October 2026:
# [1] "GSE304707_family.xml"
# XML stores information in labelled fields. Its contents have not yet been
# imported or inspected together in the user's R session.

# STEP 2B: Extract the metadata XML into the LHON output folder.
# Create a metadata subfolder to keep extracted files organized.
lhon_metadata_dir <- file.path(lhon_output_dir, "metadata")
dir.create(lhon_metadata_dir, recursive = TRUE, showWarnings = FALSE)

# Extract only the XML file identified in Step 2A.
# exdir specifies the destination; the downloaded archive stays intact.
untar(
  lhon_metadata_archive,
  files = "GSE304707_family.xml",
  exdir = lhon_metadata_dir
)

# Store the extracted file's full path for the next import step.
lhon_metadata_file <- file.path(
  lhon_metadata_dir,
  "GSE304707_family.xml"
)

# Confirm that the extracted file exists. Expected: TRUE.
# This checks file presence, not whether the XML content is valid.
file.exists(lhon_metadata_file)
# Confirmed by user on 29 October 2026: TRUE.
# The extracted file exists; its contents still need to be loaded and reviewed.

# STEP 2C: Load and preview the XML metadata, as we did for RP.
# Read the extracted file using the xml2 package.
# The :: calls its function directly, so library(xml2) is not required.
lhon_metadata_xml <- xml2::read_xml(lhon_metadata_file)

# Display an overview of the document and its labelled sections.
lhon_metadata_xml

# Confirmed by user screenshot on 29 October 2026:
# xml_document with MINiML root, version 0.5.0, and the displayed GEO namespace.
# The 20 top-level entries comprise 5 Contributor, 1 Database, 1 Platform,
# 12 Sample (GSM9155578 through GSM9155589), and 1 Series.
# These 20 entries are document sections, not 20 biological samples.

# STEP 2D: Create a table of sample identifiers and titles, as we did for RP.
# Give the XML namespace shown in the overview the short label "geo".
# This identifies the XML format for searches; it does not download anything.
lhon_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")

# Select only the Sample sections directly inside MINiML.
# Read /geo:MINiML/geo:Sample as: start at MINiML, then select its Sample sections.
# xml_find_all() returns every matching section, stored in lhon_sample_nodes.
lhon_sample_nodes <- xml2::xml_find_all(
  lhon_metadata_xml,
  "/geo:MINiML/geo:Sample",
  ns = lhon_ns
)

# Check how many sample records were selected.
length(lhon_sample_nodes)

# Build one row per sample with its GEO identifier and recorded title.
# iid holds the GSM identifier; Title contains the sample's descriptive name.
# Keep the title as "title" until we inspect how it relates to count labels.
# xml_attr(..., "iid") reads the ID attached to each Sample section.
# xml_find_first(..., "geo:Title") selects the Title within each sample.
# xml_text() reads that Title's text; data.frame() combines the two columns.
lhon_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(lhon_sample_nodes, "iid"),
  title = xml2::xml_text(
    xml2::xml_find_first(lhon_sample_nodes, "geo:Title", ns = lhon_ns)
  ),
  stringsAsFactors = FALSE
)

# Display the table for us to review together.
lhon_metadata_all

# Confirmed by user screenshot on 29 October 2026:
# 12 rows and two columns, geo_accession and title.
# GSM9155578-580: WT, replicate1 through replicate3.
# GSM9155581-583: YARS2, replicate1 through replicate3.
# GSM9155584-586: COI, replicate1 through replicate3.
# GSM9155587-589: YARS2/COI, replicate1 through replicate3.
# These are deposited sample titles, not an independent verification of
# animal replication or a completed match to the count-column labels.

# STEP 2E: Check which biological-information fields the samples contain.
# Find Characteristics fields anywhere inside the selected Sample sections.
# .// means search within those sections, including their nested sections.
lhon_characteristic_nodes <- xml2::xml_find_all(
  lhon_sample_nodes,
  ".//geo:Characteristics",
  ns = lhon_ns
)

# Read each field's label from its "tag" attribute.
# This identifies available categories before we choose fields to extract.
lhon_characteristic_tags <- xml2::xml_attr(lhon_characteristic_nodes, "tag")

# Display the field labels and count how often each occurs across the samples.
# Include missing labels, if any, rather than silently hiding them.
table(lhon_characteristic_tags, useNA = "ifany")

# Confirmed by user screenshot on 29 October 2026:
# age = 12, genotype = 12, Sex = 12, strain = 12, tissue = 12.
# No missing-tag category is displayed. Preserve the capital S in "Sex"
# when selecting this label from the XML.
# These are occurrence totals across the 12 sample records; individual
# per-sample values and coverage are still to be reviewed during extraction.
# Next: add the observed characteristic fields to the sample metadata table.

# STEP 2F: Add the recorded biological characteristics, as we did for RP.
# Keep lhon_metadata_all in its existing XML order so its rows correspond
# to lhon_sample_nodes. Use a separate matched table later for count columns.

# Define a reusable function that reads one named characteristic per sample.
# tag is the field label we request, such as "age".
# The search selects Characteristics with that exact tag inside each sample.
lhon_get_characteristic <- function(tag) {
  lhon_path <- paste0(".//geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(lhon_sample_nodes, lhon_path, ns = lhon_ns),
    trim = TRUE
  )
}
# xml_find_first() returns the first matching field for each sample;
# a missing field becomes NA. trim = TRUE removes surrounding whitespace.

# Read the five labels observed in Step 2E and add their values as columns.
# Use the exact XML spelling "Sex" while naming our table column "sex".
lhon_metadata_all$genotype <- lhon_get_characteristic("genotype")
lhon_metadata_all$sex <- lhon_get_characteristic("Sex")
lhon_metadata_all$age <- lhon_get_characteristic("age")
lhon_metadata_all$tissue <- lhon_get_characteristic("tissue")
lhon_metadata_all$strain <- lhon_get_characteristic("strain")

# Display the expanded table to inspect each sample's recorded values.
lhon_metadata_all

# Confirmed by user screenshot on 29 October 2026:
# The expanded table displays 12 rows. Genotype labels are WT, YARS2, COI,
# and YARS2/COI; all displayed sex/age/tissue values are male/12 months/retina.
# The strain column is clipped at the right edge of the supplied screenshot;
# its full values have not been visibly confirmed in this joint walkthrough.
# This updates the R object, not a saved CSV. Preserve recorded labels;
# animal independence and count-column matching remain to be established.

# STEP 2G: Check the extraction-protocol descriptions for sample timing.
# Requested by the user: determine whether extraction timing is reported.
# Age is separate from a collection/extraction date or time of day.
# Start with the extraction description; absence here would not establish
# absence from every part of the metadata or from the original paper.

# Read the recorded extraction protocol for each sample, as we did for RP.
# A missing Extract-Protocol field returns NA for that sample.
lhon_extraction_protocol <- xml2::xml_text(
  xml2::xml_find_first(
    lhon_sample_nodes, ".//geo:Extract-Protocol", ns = lhon_ns
  ),
  trim = TRUE
)

# Keep each sample's GEO identifier attached to its protocol description.
names(lhon_extraction_protocol) <- xml2::xml_attr(lhon_sample_nodes, "iid")

# Display each different description once, with readable line breaks.
# unique() simplifies only the display; the per-sample entries remain stored.
cat(unique(lhon_extraction_protocol), sep = "\n\n")
cat("\n")

# Confirmed by user screenshot on 29 October 2026:
# The displayed extraction description reports retina dissection, RNA
# isolation with the RNeasy mini kit, 2 micrograms RNA input per sample,
# and sequencing-library preparation steps. It provides no explicit
# collection/extraction date, clock time or collection-to-extraction interval.
# The user explicitly requests documenting this missing timing information.
#
# DOCUMENTED METADATA LIMITATION:
# Retinal collection and RNA-extraction timing are not reported in the
# inspected Characteristics and Extract-Protocol fields. Age is recorded
# separately as 12 months. Timing therefore remains unspecified for this
# precheck; do not infer shared extraction time or processing batch.
# This statement covers the fields reviewed together, not every other XML
# section or an unreviewed publication. No sample exclusion follows from it.

# STEP 3A: Check how the sample descriptions identify the count columns.
# Resume point confirmed on 30 October 2026: the restored metadata has
# 12 rows and seven columns, without description. Its GEO IDs match the
# reopened XML sample records in the same order (TRUE).
# The observed metadata title "WT, replicate1" differs from the count label
# "WT-1". Inspect the deposited descriptions for an explicit connection.
# lhon_metadata_all and lhon_sample_nodes remain in their original XML order.

# Read each sample's Description field and add it to the metadata table.
# A missing description returns NA in the corresponding sample's row.
lhon_metadata_all$description <- xml2::xml_text(
  xml2::xml_find_first(
    lhon_sample_nodes,
    ".//geo:Description",
    ns = lhon_ns
  ),
  trim = TRUE
)

# Display the title and description together so we can inspect the names.
lhon_metadata_all[, c("title", "description")]

# Confirmed by user screenshot on 30 October 2026:
# Each Description explicitly identifies its column in gene_count_matrix.csv.
# WT replicate1-3 -> WT-1..3; YARS2 replicate1-3 -> Y-1..3;
# COI replicate1-3 -> C-1..3; YARS2/COI replicate1-3 -> YC-1..3.
# The description has a library-name line followed by a count-column line.
# The screenshot prints its line break as \n and retains the deposited
# spelling "colunm". Neither prevents extraction of the final column label.
# The displayed two-column view does not remove other metadata columns.
# The user's second screenshot also shows the full strain value C57BL/6J
# for all 12 records, resolving the earlier cropped-display limitation.

# PAUSE CHECKPOINT — requested 29 October 2026.
# Run this block in the current R session before ejecting MaxineDrive.
# Step 3A has not yet been confirmed; save the table in its current state.

# Set the filenames for the checkpoint and a readable metadata copy.
lhon_checkpoint_file <- file.path(
  lhon_output_dir,
  "LHON_GSE304707_precheck_checkpoint_2026-10-29.RData"
)
lhon_metadata_csv <- file.path(
  lhon_output_dir,
  "LHON_GSE304707_metadata_2026-10-29.csv"
)

# Save the imported table, current metadata, protocol text and file paths.
# XML document/node objects contain pointers; reopen the XML when resuming.
save(
  lhon_raw, lhon_metadata_all, lhon_extraction_protocol,
  lhon_characteristic_tags, lhon_ns,
  lhon_dir, lhon_count_file, lhon_output_dir,
  lhon_metadata_archive, lhon_metadata_dir, lhon_metadata_file,
  lhon_checkpoint_file, lhon_metadata_csv,
  file = lhon_checkpoint_file
)

# Also save the current metadata as a readable CSV table.
write.csv(lhon_metadata_all, lhon_metadata_csv, row.names = FALSE)

# Confirm both files exist before ejecting the drive.
c(
  checkpoint = file.exists(lhon_checkpoint_file),
  metadata = file.exists(lhon_metadata_csv)
)
# Expected: TRUE for both, after the save commands finish without errors.
# The dated pause checkpoint was subsequently confirmed saved and restored
# on 30 October. It predates the newly added description/matching work.

# STEP 3B: Match metadata rows to the count-column names.
# Take the label after the last ": " in each observed description.
# trimws() removes surrounding whitespace from the resulting sample ID.
lhon_metadata_all$sample_id <- trimws(
  sub(".*: ", "", lhon_metadata_all$description)
)

# Read the count-column names; the first column contains gene_name.
lhon_sample_ids <- colnames(lhon_raw)[-1]

# Check for missing/duplicated metadata IDs or duplicated count names,
# and require the same set of IDs in both sources before arranging rows.
# stopifnot() stops here if any of these checks fails.
stopifnot(
  !anyNA(lhon_metadata_all$sample_id),
  anyDuplicated(lhon_metadata_all$sample_id) == 0L,
  anyDuplicated(lhon_sample_ids) == 0L,
  setequal(lhon_metadata_all$sample_id, lhon_sample_ids)
)

# Create a separate table with metadata rows in count-column order.
# Keep lhon_metadata_all in its original XML order.
lhon_metadata_matched <- lhon_metadata_all[
  match(lhon_sample_ids, lhon_metadata_all$sample_id),
]
rownames(lhon_metadata_matched) <- NULL

# Confirm that the matched metadata IDs follow the count columns exactly.
identical(lhon_metadata_matched$sample_id, lhon_sample_ids)

# Display a short view of the result; all metadata columns remain stored.
lhon_metadata_matched[, c("geo_accession", "sample_id", "genotype")]

# Expected: TRUE and 12 matched rows. Actual matching results are pending
# user execution/review. The assistant checked label parsing only on a small
# constructed example using the layout shown in the user's screenshot.
# It did not run the matching or read underlying dataset values.
# Next after successful matching: prepare/check the numeric count matrix,
# then calculate total counts. The pause checkpoint has not been refreshed.
