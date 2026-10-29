# Portable repository copy. Start R in the repository root.
# This file records the interactive workflow; later result/closure notes
# supersede earlier pending notes. See the execution guides under docs/.
# Original working files are preserved. Repository sync did not run analysis.
# Retinitis pigmentosa GSE156533: manual precheck
# Started: 28 October 2026, after the user closed the glaucoma precheck.
# Work one small block at a time and review the author's output before advancing.
# All new objects use rp_ names. Do not clear the workspace or alter glaucoma
# objects. Save the latest glaucoma checkpoint before starting this script.
# Preserve downloaded filenames; each future figure must include save code.

# STEP 0: RP-specific location for this script and all subsequent outputs.
# The author created this folder inside prequalitychecks on 28 October 2026.
# Use file.path(rp_output_dir, <filename>) for RP figures, tables and checkpoints.
# Input downloads remain separate. Start R in the repository root.
rp_output_dir <- file.path("results", "retinitis_pigmentosa")
dir.create(rp_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.exists(rp_output_dir)
# Expected: TRUE.

# STEP 1: Identify the actual expression file and metadata archive.
rp_dir <- Sys.getenv("RP_DATA_DIR", unset = file.path("data", "raw", "Retinitis_Pigmentosa_GSE156533"))

rp_count_file <- file.path(
  rp_dir,
  "GSE156533_readcounts_organized_by_largest_to_smallest_count.xlsx"
)

# This existing download is a tar archive despite having no filename suffix.
# File inspection confirmed it contains GSE156533_family.xml (GEO MINiML).
# It is not a plain SOFT file; do not readLines() the archive or rename it.
rp_metadata_archive <- file.path(rp_dir, "GSE156533_family")

# Check the folder and both input files. Expected: all three TRUE.
c(
  rp_folder = dir.exists(rp_dir),
  rp_counts = file.exists(rp_count_file),
  rp_metadata = file.exists(rp_metadata_archive)
)

# List the metadata archive's contents without extracting/modifying files.
# Expected member: GSE156533_family.xml.
utils::untar(rp_metadata_archive, list = TRUE)

# STEP 0/1 RESULT — 28 October 2026:
# User reports all four combined input/output checks TRUE and archive listing
# [1] "GSE156533_family.xml". File locations and archive inventory are confirmed.
# These results do not yet assess count values or biological sample quality.

# STEP 2: Inspect sheet names and a small first-sheet preview.
# readxl 1.5.0 is installed in the inspected local R environment.
# Run manually, then review the output before importing a final count matrix.
rp_sheets <- readxl::excel_sheets(rp_count_file)
rp_sheets

rp_preview <- readxl::read_excel(
  rp_count_file,
  sheet = rp_sheets[1],
  col_names = FALSE,
  col_types = "text",
  n_max = 6
)

print(rp_preview, width = Inf)

# excel_sheets() lists the workbook tabs. The first sheet is being previewed,
# not yet selected as the final expression input.
# col_names = FALSE keeps the potential header row among the displayed values;
# automatic column labels such as ...1 and ...2 are temporary preview labels.
# col_types = "text" prevents type guessing across headers and data in this
# inspection object. The final expression counts will later be numeric.
# n_max = 6 limits the preview to six rows, with leading empty rows skipped by
# readxl. width = Inf asks the printed tibble to show every imported column.
# Expected: sheet name(s) and a small table exposing labels and example values.
# Await the author's outputs; do not assume which columns are sample counts.
# No user workbook data were loaded by the assistant. No archive extraction,
# full matrix import, expression QC, figure export or checkpoint save performed.

# STEP 2 RESULT / EXPLANATION — 28 October 2026:
# User screenshot shows a 6 x 17 tibble: the header row plus five gene rows.
# This is the limited preview, not the full workbook dimensions. Sheet-name
# output is not visible in this screenshot, so no sheet name is recorded here.
# ...1 through ...17 are temporary labels; the real headers are in preview row 1.
# All <chr> tags reflect col_types = "text", including the count values.
# The wrapped blocks are continuations of the same 17-column table.
# Column 1: gene_id (mouse Ensembl gene IDs).
# Columns 2:8, in workbook order:
# P23H_6, WT_7, P23H_8, WT_9, P23H_11, WT_12, P23H_13.
# These provide four P23H-labelled and three WT-labelled sample count columns.
# Genotype meaning is consistent with previously retrieved GEO study metadata;
# exact sample matching, sex and experimental-unit verification remain pending.
# Columns 9:17 are annotations: gene_name, gene_chr, gene_start, gene_end,
# gene_strand, gene_length, gene_biotype, gene_description, tf_family.
# Numeric annotation fields such as genomic coordinates and gene length must
# not enter the expression count matrix. Retain annotations separately.
# Visible genes: Rho, Gnat1, mt-Co1, Gnb1, Prph2.
# Example: Rho has reported count 224353 in P23H_6 and 555127 in WT_7.
# These unnormalised single-library values do not establish differential
# expression, biological causes or a quality verdict.
# In tf_family, '-' is a placeholder for no listed annotation, not a count zero.
# gene_length is the supplied length annotation; its calculation is not yet
# verified and must not be assumed equal to gene_end - gene_start + 1.
# Interpretation only: no new executable analysis block or workbook import
# was run by the assistant. Next: finish sheet-name context, then import with
# the observed header and appropriate types before full-table integrity checks
# and metadata matching, continuing one explained block at a time.

# STEP 3: Import the complete first-sheet table using its actual header.
# The author requests loading the count file into an R object.
# rp_count_file is already the filepath; retain it as a filepath.
# rp_raw will contain all gene rows, seven sample columns and annotations.
# The preview established a header in the first non-empty row of this sheet.
# rp_sheets[1] selects the same sheet that produced the author's preview; its
# literal name has not yet been reported. Display it here for provenance.
rp_sheets[1]

rp_raw <- readxl::read_excel(
  rp_count_file,
  sheet = rp_sheets[1],
  col_names = TRUE
)

dim(rp_raw)
colnames(rp_raw)
print(head(rp_raw[, 1:8], n = 5), width = Inf)

# col_names = TRUE uses the header as column names, not as a gene row.
# Omitting n_max removes the six-row preview limit and imports all data rows.
# Omitting col_types = "text" restores readxl's default type detection.
# Confirm actual column types and import warnings before creating a numeric
# count matrix; type detection is not itself a passed numeric integrity check.
# Expected: 17 columns and more than five gene rows; review the actual total.
# First eight columns are the observed gene_id plus seven sample count fields.
# This display does not remove annotation columns from rp_raw.
# The <- assignment stores the table in the active R session. It does not
# create a new disk file or save a persistent RData checkpoint.
# Await user execution and results. The assistant has not imported the author's
# workbook or created a counts-only matrix or matched sample metadata.

# STEP 3 RESULT — 28 October 2026:
# User screenshot confirms rp_raw has 54532 rows and 17 columns.
# Actual column names match the preview: gene_id, seven sample columns,
# and nine gene-annotation columns. The first-five-row display is 5 x 8
# only because head() and [, 1:8] limit the display; rp_raw remains 54532 x 17.
# gene_id is <chr>; all seven sample columns are <dbl> (numeric storage).
# <dbl> does not by itself establish whole-valued, finite, nonnegative counts.
# Sheet-name output and the import call/warnings are not shown in this image.
# Do not claim unique genes or completed numeric integrity checks yet.

# STEP 4: Check missing identifiers/counts, duplicate IDs and numeric types.
# The observed columns 1:8 are gene_id plus the seven sample-count columns.
# Restrict this missing-value check to these essential fields; optional gene
# annotation gaps are a separate issue and should not be mistaken for missing
# expression counts. This indexing does not remove anything from rp_raw.
anyNA(rp_raw[, 1:8])

# Return 0 if all gene IDs are distinct; otherwise the first duplicate position.
anyDuplicated(rp_raw$gene_id)

# Apply is.numeric to each of the seven sample-count columns.
sapply(rp_raw[, 2:8], is.numeric)

# Expected: FALSE, then 0, then TRUE for all seven sample columns.
# Report the actual outputs; do not silently discard rows or coerce values.
# These are structural input checks. Finite/nonnegative/whole-valued checks
# follow before sample-level QC. The assistant has not executed this block,
# loaded user expression data, matched XML metadata or saved a checkpoint.

# STEP 4 RESULT — 28 October 2026:
# User reports FALSE, 0, and TRUE for all seven sample columns.
# Thus there are no NA values in gene_id/sample counts, no duplicate gene_id
# entries, and all seven count columns are numeric. Annotation completeness
# was not tested by this block. No samples or genes have been filtered.

# STEP 5: Build the numeric count matrix and check the count values.
# Columns 2:8 are the seven verified sample columns in their original order.
# The original table, including every gene annotation, stays in rp_raw.
rp_counts_all <- as.matrix(rp_raw[, 2:8])
rownames(rp_counts_all) <- rp_raw$gene_id

# Expected: 54532 gene rows x 7 sample columns.
dim(rp_counts_all)

# Check that every entry is finite (no NA, NaN, Inf or -Inf).
all(is.finite(rp_counts_all))

# Display the minimum and maximum. The minimum should be nonnegative;
# zero is a valid count. Review the actual range rather than assuming it.
range(rp_counts_all)

# Check whether every entry is a whole number. This compares to rounded
# values without assigning them, so it does not change or round the matrix.
all(rp_counts_all == round(rp_counts_all))

# Expected: dimensions 54532 x 7, finite TRUE, nonnegative minimum, whole TRUE.
# A different output should be reviewed, not repaired by rounding/deleting.
# Gene IDs are row labels, not an extra numeric column. Sample column names
# are retained by as.matrix(). No normalization or sample filtering is done.
# Await user output. Matrix construction and checks have not been run by the
# assistant; XML sample metadata and saved workspace remain pending.

# STEP 5 RESULT — 28 October 2026:
# User screenshot confirms dimensions 54532 x 7; finite TRUE;
# count range 0 to 670529; whole-valued TRUE.
# The numeric matrix therefore contains finite, nonnegative whole counts.
# 670529 is the largest individual gene/sample entry, not a library total.
# These input-format checks do not establish sample quality or independence.

# STEP 6A: Extract and open the XML metadata for manual inspection.
# xml2 1.5.2 is installed in the inspected local R environment.
# Keep the downloaded archive intact; extract only its known XML member into
# an RP-specific metadata subfolder. Rerunning extracts the same member again.
rp_metadata_dir <- file.path(rp_output_dir, "metadata")
dir.create(rp_metadata_dir, showWarnings = FALSE)

utils::untar(
  rp_metadata_archive,
  files = "GSE156533_family.xml",
  exdir = rp_metadata_dir
)

rp_metadata_file <- file.path(rp_metadata_dir, "GSE156533_family.xml")
file.exists(rp_metadata_file)

rp_xml <- xml2::read_xml(rp_metadata_file)
rp_xml

# Expected: file.exists TRUE, followed by an XML document overview with a
# MINiML root. XML stores labelled sections, not a flat sample table yet.
# The file is written to:
# results/retinitis_pigmentosa/metadata/GSE156533_family.xml
# rp_metadata_file remains a path; rp_xml is the parsed XML object in R.
# Await the overview before extracting sample IDs/titles/characteristics.
# No metadata extraction, XML parsing, new user-data calculation or checkpoint
# save has been performed by the assistant. User execution is pending.

# STEP 6A RESULT — 28 October 2026:
# User screenshot confirms the extracted file exists (TRUE) and rp_xml is
# an xml_document with MINiML root, version 0.5.0, using the namespace
# http://www.ncbi.nlm.nih.gov/geo/info/MINiML.
# There are 13 displayed top-level sections: 3 Contributor, 1 Database,
# 1 Platform, 7 Sample (GSM4733112 through GSM4733118), and 1 Series.
# These are not 13 biological samples. The XML is loaded but not yet tabulated.
# Raw XML line inspection confirms Sample identifiers use attribute "iid"
# and sample names are in each Sample's direct Title child.

# STEP 6B: Extract sample accessions and titles into a small metadata table.
# A namespace identifies the XML vocabulary. Give the observed GEO namespace
# the short prefix geo for searches; this URL is an identifier, not a download.
rp_ns <- c(geo = "http://www.ncbi.nlm.nih.gov/geo/info/MINiML")

# Select only Sample sections directly inside the MINiML root.
rp_sample_nodes <- xml2::xml_find_all(
  rp_xml,
  "/geo:MINiML/geo:Sample",
  ns = rp_ns
)

length(rp_sample_nodes)

# Preserve the exact GEO title as sample_id for comparison to count columns.
# xml_find_first searches separately within each Sample and preserves a missing
# entry as NA rather than shortening/reordering the result if a Title is absent.
rp_metadata_all <- data.frame(
  geo_accession = xml2::xml_attr(rp_sample_nodes, "iid"),
  sample_id = xml2::xml_text(
    xml2::xml_find_first(rp_sample_nodes, "geo:Title", ns = rp_ns)
  ),
  stringsAsFactors = FALSE
)

rp_metadata_all

# Expected: 7 sample sections and a seven-row, two-column table containing
# GSM accessions and the P23H/WT sample titles. This is the beginning of the
# metadata table; genotype, sex, age and tissue have not been extracted yet.
# Retain source XML order at this step; do not assume it matches count order.
# After reviewing the table, match explicitly by names and verify uniqueness
# and exact ordering before associating biological metadata with the counts.
# The assistant has inspected XML text structure only, not executed this
# extraction code or generated the table. Await user output.

# STEP 6B RESULT — 28 October 2026:
# User screenshot confirms 7 sample nodes and the seven-row identity table:
# GSM4733112 P23H_6; GSM4733113 P23H_8; GSM4733114 P23H_11;
# GSM4733115 P23H_13; GSM4733116 WT_7; GSM4733117 WT_9;
# GSM4733118 WT_12.
# The author asks where the biological metadata is. Explain that the last code
# copied only identity fields. Biological details remain in rp_xml beneath
# each Sample's Channel/Characteristics elements until explicitly extracted.
# The table is currently in source XML order, different from count order.

# STEP 6C: Add the recorded biological sample characteristics.
# Run after Step 6B, while rp_metadata_all and rp_sample_nodes share XML order.
# Keep rp_metadata_all in source order; use a separate matched object later.
# The exact tag spellings below were verified in all seven Sample records.
# This helper returns one field per sample, retaining source labels while
# trimming only surrounding whitespace/newlines. Missing tags remain NA.
rp_get_characteristic <- function(tag) {
  rp_path <- paste0("geo:Channel/geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(rp_sample_nodes, rp_path, ns = rp_ns),
    trim = TRUE
  )
}

rp_metadata_all$genotype <- rp_get_characteristic("genotype/variation")
rp_metadata_all$sex <- rp_get_characteristic("Sex")
rp_metadata_all$age <- rp_get_characteristic("age")
rp_metadata_all$tissue <- rp_get_characteristic("tissue")
rp_metadata_all$strain <- rp_get_characteristic("strain background")

rp_metadata_all

# Expected: 7 rows and 7 fields: geo_accession, sample_id, genotype, sex, age,
# tissue, strain. These fields are read from the source, not inferred from IDs.
# Genotype describes the genetic state; strain is the mouse genetic background.
# This covers the requested biological characteristics, not every XML field:
# extraction protocols and processing/design statements remain in rp_xml.
# Inspect this output before matching to count columns or summarising design.
# The assistant inspected raw XML structure only and did not execute this
# extraction code, create the expanded table, alter counts or save a checkpoint.

# STEP 6C RESULT — 28 October 2026:
# User screenshot confirms the expanded seven-row metadata table.
# P23H_6/P23H_8 are female; P23H_11/P23H_13 are male.
# WT_7/WT_9 are female; WT_12 is male.
# All samples: Postnatal day 90; Whole retina; C57BL/6J.
# Mutant label: RhoP23H/WT (mutant); control: RhoWT (wild-type).
# These are deposited labels. Exact count order and biological independence
# have not yet been verified. Source order remains four P23H followed by WT.

# STEP 6D: Inventory ALL characteristic labels, without preselecting field names.
# The author explicitly asks to check other characteristics because they may be
# relevant. Read-only audit of the complete downloaded XML found exactly the
# same five characteristic tags once per sample, across all seven samples.
# This audit read metadata only; it did not execute user R code or analyse counts.
rp_characteristic_nodes <- xml2::xml_find_all(
  rp_sample_nodes,
  ".//geo:Characteristics",
  ns = rp_ns
)

rp_characteristic_tags <- xml2::xml_attr(rp_characteristic_nodes, "tag")
table(rp_characteristic_tags, useNA = "ifany")

# Expected: Sex, age, genotype/variation, strain background, tissue;
# each should occur 7 times. This aggregates entries; the assistant's separate
# source inspection confirmed once per sample, not just 35 entries overall.
# No extra named Characteristics field is present in this XML. This does NOT
# imply these five fields exhaust the biologically relevant study metadata.
# Other source sections contain:
# - Source and Organism: genotype-specific whole-retina source; Mus musculus.
# - Molecule: total RNA; Treatment-Protocol: None, animals never anesthetized.
# - Extract-Protocol: collection/freezing, RNA extraction, library preparation.
# - Library/Instrument: RNA-Seq, transcriptomic, cDNA; Illumina NovaSeq 6000.
# - Data-Processing: read filtering, mm10/TopHat2 and FPKM reporting.
# - Overall-Design: heterozygous P23H n=4 vs WT littermates n=3 at P90.
# Treatment, extraction, processing and instrument texts are identical across
# the seven XML records. Equal descriptions do not prove one processing batch.
# No explicit batch, litter ID, eye/laterality, collection-time-of-day or pooling
# field is recorded here. Missing reporting must not be translated into "none".
# Review original-paper methods for biological units/pooling and design details.
# The processing paragraph mentions FPKM, while Series separately links count
# and FPKM workbooks. The current input is the deposited readcounts workbook;
# this paragraph alone does not identify the count-generation algorithm.
# Await the author's inventory result; then inspect non-Characteristics protocol
# and study-design statements in R. Do not finalize the model or invent covariates.

# STEP 6D RESULT / FILE-FORMAT CLARIFICATION — 28 October 2026:
# User screenshot confirms age, genotype/variation, Sex, strain background,
# and tissue each appear 7 times. No additional Characteristics tag is shown.
# The earlier full-source audit independently confirmed one of each per sample.
# User asks why another XML appeared and whether a SOFT file contains many.
# Clarify: RP's GSE156533_family download is a tar archive, not a SOFT file.
# Its inventory contains one member: GSE156533_family.xml. Extraction unpacked
# that existing member into the RP metadata subfolder; it did not create a new
# or different study metadata source. rp_xml is the in-memory parsed document.
# The same XML contains 7 Sample records plus Platform/Series/etc. sections.
# Characteristics, Treatment-Protocol and Extract-Protocol are different
# sections within each sample record, not separate XML files. Overall-Design
# belongs to the Series section in this same document. Only Characteristics
# labels were inventoried by Step 6D, so protocol fields are not in that table.
# Glaucoma used a plain-text .soft family file: a different GEO format, also
# capable of holding multiple related records within one file. Do not treat
# RP's archive as SOFT or describe multiple records as multiple physical files.
# Conceptual explanation only; no new executable code, extraction or analysis.

# STEP 6E: Read and display the treatment and extraction protocols.
# User explicitly requests viewing these two sections of the existing XML.
# Each vector holds one recorded protocol per sample, in rp_sample_nodes order.
rp_treatment_protocol <- xml2::xml_text(
  xml2::xml_find_first(
    rp_sample_nodes, "geo:Channel/geo:Treatment-Protocol", ns = rp_ns
  ),
  trim = TRUE
)

rp_extraction_protocol <- xml2::xml_text(
  xml2::xml_find_first(
    rp_sample_nodes, "geo:Channel/geo:Extract-Protocol", ns = rp_ns
  ),
  trim = TRUE
)

# Retain the source accession with each entry for subsequent reference.
names(rp_treatment_protocol) <- xml2::xml_attr(rp_sample_nodes, "iid")
names(rp_extraction_protocol) <- xml2::xml_attr(rp_sample_nodes, "iid")

# Display every distinct description, preserving line breaks for readability.
# unique() affects only the displayed result; the seven-entry vectors remain.
cat("TREATMENT PROTOCOL\n")
cat(unique(rp_treatment_protocol), sep = "\n\n")

cat("\n\nEXTRACTION PROTOCOL\n")
cat(unique(rp_extraction_protocol), sep = "\n\n")
cat("\n")

# Expected from the earlier complete XML audit: one distinct treatment text
# and one distinct extraction text, each shared by all seven records.
# Treatment: None. Animals never anesthetized.
# Extraction: cervical dislocation; retinas removed and snap-frozen in liquid
# nitrogen; Qiagen RNAeasy Mini-Kit; genomic DNA removed by on-column digestion;
# standard Illumina library preparation, with RNA-seq supplied by Novogene.
# Preserve the original text in the vectors. These descriptions do not specify
# an animal/eye-to-library mapping or rule out unreported pooling/batch effects.
# The assistant has not run this R block. Await the author's protocol output
# before interpreting it together; no count analysis or checkpoint save occurs.

# PROTOCOL REVIEW / METADATA STORAGE CLARIFICATION — 28 October 2026:
# User says the protocol descriptions reveal nothing further relevant and
# explicitly requests proceeding to metadata/count matching. No protocol
# output screenshot was supplied; record this as the author's review, not a new
# independently reproduced R result. Source audit found no added grouping field.
# The five $<- assignments in Step 6C persist in rp_metadata_all in the current
# R session. This is a 7 x 7 table, not a display-only transformation. It has
# not yet been exported to CSV or saved in a confirmed RP checkpoint.
# No additional Characteristics tags were found; original-paper biological-unit
# verification remains part of later design review, not a block on ID matching.

# STEP 7: Match metadata rows to count-matrix columns by exact sample names.
# Preserve rp_metadata_all in XML order so earlier extraction code remains safe.
# The verified XML and workbook orders differ, so do not join by row position.
# Stop on missing/duplicated names or different sample membership.
stopifnot(
  !anyNA(rp_metadata_all$sample_id),
  anyDuplicated(rp_metadata_all$sample_id) == 0L,
  anyDuplicated(colnames(rp_counts_all)) == 0L,
  setequal(rp_metadata_all$sample_id, colnames(rp_counts_all))
)

rp_metadata_matched <- rp_metadata_all[
  match(colnames(rp_counts_all), rp_metadata_all$sample_id),
  , drop = FALSE
]

# Reset displayed row numbers; retain every metadata field.
rownames(rp_metadata_matched) <- NULL

dim(rp_metadata_matched)
identical(rp_metadata_matched$sample_id, colnames(rp_counts_all))
rp_metadata_matched

# stopifnot is silent when checks pass. If it errors, inspect before continuing.
# Expected dimensions: 7 x 7. Exact name-and-order check: TRUE.
# Expected row order: P23H_6, WT_7, P23H_8, WT_9, P23H_11, WT_12, P23H_13.
# This reorders a copy of the metadata only; count values/column order and
# original source-order metadata are preserved. Counts remain a numeric matrix
# separate from the descriptive table. No flat merge or sample exclusion occurs.
# Await user output. The assistant has not executed matching, changed the live
# R workspace, run expression calculations or saved a checkpoint.

# STEP 7 RESULT — 28 October 2026:
# User screenshot shows rp_metadata_matched with seven rows/seven fields in
# order P23H_6, WT_7, P23H_8, WT_9, P23H_11, WT_12, P23H_13. This visually
# matches the count-column order previously supplied. GEO accessions and all
# five characteristics remain correctly attached to each sample.
# The separate dim()/identical() outputs are not visible in this screenshot;
# do not claim a user-reported TRUE yet. Include the exact check below.

# STEP 8: Confirm alignment and summarise the experimental groups.
identical(rp_metadata_matched$sample_id, colnames(rp_counts_all))

with(rp_metadata_matched, table(genotype, sex))

unique(rp_metadata_matched[, c("age", "tissue", "strain")])

# Expected: TRUE; mutant 2 female/2 male, WT 2 female/1 male; one common
# age/tissue/strain row (Postnatal day 90, Whole retina, C57BL/6J).
# These counts describe deposited sample libraries; retain sex for subsequent
# design/QC interpretation. Age is a common collection stage in this dataset,
# not multiple time-point groups. Do not drop sex or fit a model in this step.
# Source review: original eLife paper (doi:10.7554/eLife.59422), RNA isolation,
# quantitative PCR and RNA-sequencing methods, confirms the 3-month cohort
# has 3 WT and 2 female + 2 male P23H samples. Its sexes-pooled wording follows
# sample QC and describes combining sexes in subsequent genotype analyses;
# this is not evidence of physical RNA pooling before library preparation.
# Exact eyes-per-library mapping is not explicit in that methods paragraph.
# The main methods also mention HiSeq/RNAlater, while this cohort's GEO record
# lists NovaSeq/liquid-nitrogen freezing; preserve both sources without guessing
# whether differing cohorts explain this. No input/model change follows here.
# Await user results. The assistant has not executed these R checks, calculated
# sample totals, normalised counts or saved a checkpoint.

# STEP 8 RESULT — 28 October 2026:
# User screenshot explicitly confirms identical(...) TRUE.
# Genotype/sex table confirms mutant 2 female + 2 male and WT 2 female + 1 male.
# Shared age/tissue/strain output is one row: Postnatal day 90, Whole retina,
# C57BL/6J. Metadata/count name-and-order alignment is now confirmed in R.
# Keep all seven libraries for the next within-study precheck summaries.

# STEP 9A: Total gene-matrix counts per sample.
# colSums adds every gene-row count within each sample column.
# These are library totals from the supplied count matrix, not necessarily
# total sequencer reads, uniquely mapped reads or a complete raw-read QC metric.
rp_total_counts <- colSums(rp_counts_all)

# Display in the original matched sample order.
rp_total_counts

# Also display from smallest to largest, retaining sample labels.
sort(rp_total_counts)

# Expected: seven named totals. Actual values await the author's output.
# sort() here changes only the display; it does not reorder rp_total_counts,
# rp_counts_all or rp_metadata_matched. No genes/samples are removed.
# Compare the spread before proceeding to detected-gene summaries. A lower
# total alone is not a sample failure; interpret with later profile checks.
# The assistant has not computed the author's totals or saved a checkpoint.

# STEP 9A RESULT — 28 October 2026:
# User screenshot reports:
# P23H_6  15469799; WT_7 16695962; P23H_8 18022977; WT_9 17758147;
# P23H_11 14709319; WT_12 18315284; P23H_13 17401717.
# Lowest: P23H_11 (mutant male), 14709319. Highest: WT_12 (WT male), 18315284.
# Ratio from the user-reported summaries: 1.2451483 (highest about 24.5% larger).
# No extreme low-total sample is apparent relative to these seven totals;
# this is descriptive, not a universal pass threshold or a quality certificate.
# Retain all seven while reviewing detected genes and later sample profiles.
# The assistant computed only the ratio of supplied summary values, not totals
# from the expression matrix. No checkpoint was saved.

# STEP 9B: Count gene rows with a nonzero count in each sample.
# For this initial descriptive summary only, detected means count > 0.
# This is not the later statistical expressed-gene filtering criterion.
rp_detected_genes <- colSums(rp_counts_all > 0)

rp_detected_genes
sort(rp_detected_genes)

# The > 0 comparison produces TRUE/FALSE entries. colSums treats TRUE as 1
# and FALSE as 0, counting positive-count gene rows in each sample column.
# Example: counts 0, 5, 100 have total count 105 but 2 detected gene rows.
# Expected: seven named gene-row totals, each between 0 and 54532.
# Actual outputs await user execution. No filtering, normalization, plotting,
# sample exclusion or checkpoint save occurs in this block.
# A zero count is no observed assigned count here, not proof that a gene can
# never be expressed biologically. Very small positive counts may be weak.
# The assistant has not executed this calculation on the author's data.

# STEP 9B RESULT — 28 October 2026:
# User screenshot reports detected gene rows (count > 0):
# P23H_6 25694; WT_7 24518; P23H_8 25094; WT_9 24550;
# P23H_11 24360; WT_12 24771; P23H_13 25317.
# Minimum P23H_11 = 24360; maximum P23H_6 = 25694. Difference 1334,
# highest about 5.48% larger (arithmetic from supplied summaries only).
# P23H_11 has the lowest count and detection totals, but no extreme deficit
# appears from these summaries alone. P23H_6 has the most detected gene rows
# despite the second-lowest count total; these summaries measure different things.
# Every detected total is below 54532/2, so over half the rows have zero counts
# in every sample. Initial all-row boxplot medians will reflect zero-count rows.
# Keep all seven samples and continue planned distribution/profile checks.

# STEP 9C: Initial all-gene log-CPM boxplot with display and complete PNG export.
# edgeR 4.8.2 and the cpm.default signature were confirmed locally.
# CPM scales counts by each library's count total; log2 compresses the range.
# prior.count = 2 is an average prior count scaled with library size to keep
# zero values finite on the log scale. It is not +2 CPM added to every entry.
# This uses every gene row before expression filtering and before TMM scaling.
rp_log_cpm <- edgeR::cpm(
  rp_counts_all,
  log = TRUE,
  prior.count = 2
)

# Define the figure once so screen and saved versions use identical settings.
rp_draw_boxplot <- function() {
  boxplot(
    rp_log_cpm,
    las = 2,
    col = "lightblue",
    outline = FALSE,
    main = "RP retinal libraries: expression distributions",
    xlab = "Sample",
    ylab = "Log2 counts per million (CPM)"
  )
}

rp_draw_boxplot()

rp_boxplot_file <- file.path(
  rp_output_dir,
  "RP_GSE156533_logCPM_boxplot.png"
)

png(rp_boxplot_file, width = 2400, height = 1800, res = 300)
rp_draw_boxplot()
dev.off()
file.exists(rp_boxplot_file)

# Full target:
# results/retinitis_pigmentosa/RP_GSE156533_logCPM_boxplot.png
# Expected: on-screen boxplot, saved PNG, and TRUE from file.exists.
# Rerunning replaces this PNG. A saved figure is not an R workspace checkpoint.
# Each box summarises all gene-row values within one sample. The box shows the
# middle 50%, its line the median; outline = FALSE hides individual outlier dots
# without removing gene rows or altering the box/whisker calculations.
# Medians at the low zero-count floor are expected here, not automatic failures.
# Initial CPM scaling is not TMM or differential-expression analysis. Counts
# and metadata stay unchanged. No figure or log-CPM values were generated by
# the assistant; await user execution and interpretation of the actual plot.

# STEP 9C CREATION / SAVE STATUS — 28 October 2026:
# User supplies a screenshot of the RP log-CPM boxplot and asks for saving code.
# Filesystem inspection confirms RP_GSE156533_logCPM_boxplot.png already exists
# in the RP output folder (~162 KB, timestamp 15:51). Creation/export should
# no longer be described as pending. Detailed plot interpretation is pending.
# Reuse the PNG export block above with the full RP output directory; running
# it again replaces that PNG. The assistant has not redrawn or recomputed it,
# and no workspace checkpoint is confirmed by the presence of this image.

# STEP 9C INTERPRETATION — 29 October 2026:
# The author supplied the saved PNG for interpretation. Seven overall
# distributions are broadly similar, with no isolated sample showing a large
# distribution shift. P23H samples have modestly higher upper box edges than
# WT samples. This does not establish global upregulation or differential
# expression; distributions do not preserve individual gene identities.
# Each sample has fewer than 27266 positive-count rows out of 54532, so more
# than half of its rows are zero. The median and lower box edge coincide at
# the transformed zero-count floor (approximately -3.1). The finite floor
# arises from the library-size-scaled prior count used by edgeR::cpm.
# Negative log-CPM values are not negative counts or evidence of downregulation.
# Whiskers use the standard 1.5-IQR rule and do not show the highest expression
# value when more extreme gene values exist. outline=FALSE hides these dots.
# Retain all seven libraries at this stage. The plot supports continuing QC,
# not a declaration that all samples pass or that genotype has no effect.
# Initial CPM uses column totals; no gene filtering or TMM has been applied.
# Next distribution view: density curves, followed by the planned preparation
# for gene-profile comparisons. No additional analysis was executed here.
# The author has also restored the objects and saved the explicit RP checkpoint
# on 29 October; this supersedes the earlier pending checkpoint status.

# STEP 9D: Initial density plot of the same all-gene log-CPM values.
# Uses base R density(), following the density-plot approach illustrated in
# the official RNAseq123 workflow. No filtering or TMM is performed here.
# Each curve summarises gene values within one sample, not gene identities.
rp_density_curves <- lapply(
  seq_len(ncol(rp_log_cpm)),
  function(i) density(rp_log_cpm[, i])
)
names(rp_density_curves) <- colnames(rp_log_cpm)

rp_density_colors <- hcl.colors(ncol(rp_log_cpm), palette = "Dark 3")

rp_draw_density <- function() {
  # Include the full plotted range of all seven curves.
  rp_xlim <- range(sapply(rp_density_curves, function(d) range(d$x)))
  rp_ylim <- c(0, max(sapply(rp_density_curves, function(d) max(d$y))))

  plot(
    rp_density_curves[[1]],
    xlim = rp_xlim, ylim = rp_ylim,
    col = rp_density_colors[1], lwd = 2,
    main = "RP retinal libraries: expression density",
    xlab = "Log2 counts per million (CPM)",
    ylab = "Density"
  )

  for (i in 2:length(rp_density_curves)) {
    lines(rp_density_curves[[i]], col = rp_density_colors[i], lwd = 2)
  }

  legend(
    "topright",
    legend = names(rp_density_curves),
    col = rp_density_colors, lty = 1, lwd = 2,
    bty = "n", cex = 0.85
  )
}

# Display the plot.
rp_draw_density()

# Save the same plot.
rp_density_file <- file.path(
  rp_output_dir,
  "RP_GSE156533_logCPM_density.png"
)
png(rp_density_file, width = 2400, height = 1800, res = 300)
rp_draw_density()
dev.off()

file.exists(rp_density_file)
# Expected: TRUE
# Full target:
# results/retinitis_pigmentosa/RP_GSE156533_logCPM_density.png
# Each curve integrates to approximately one. Density is relative concentration
# along the log-CPM axis, not a raw number of genes or an expression level.
# Expect a prominent low-value peak from the many zero-count gene rows.
# Kernel smoothing can extend the curve below the smallest observed value;
# this extension is not additional observed expression.
# Use all curves to set bounds rather than clipping them to the first sample.
# Defaults use automatic per-sample smoothing, as in the referenced workflow.
# Creation/save and interpretation await user execution; assistant only
# prepared the code. The existing workspace checkpoint remains at boxplot stage.

# STEP 9D RESULT / INTERPRETATION — 29 October 2026:
# User supplies the density plot and questions whether lack of obvious sample
# outliers signals an error, including in the earlier glaucoma work.
# Density PNG confirmed on disk: 238680 bytes, 29 October 10:57 timestamp.
# All seven curves overlap broadly, with a dominant peak around -3.1 and
# a broader smaller hump around log-CPM 4-6, followed by a right tail.
# The low peak mainly reflects the many zero-count rows; differences in peak
# heights also depend on the automatic smoothing bandwidth and are not direct
# gene percentages or sample-quality scores. No isolated large distribution
# shift is visible. The curves are similar, not identical.
# Reviewed code confirms separate density(rp_log_cpm[, i]) calculations for
# all seven columns. CPM scaling does not force identical distributions.
# Earlier outline=FALSE only hides individual gene dots in the boxplot; it
# does not remove sample columns or exclude genes from these density curves.
# No clear marginal-distribution outlier is identified; this does not imply
# no sample-profile outlier, no technical errors, or no disease response.
# Glaucoma R-08 was distinctive in MDS/PCA and had lower totals/detection;
# it was retained because failure was not established, per the prior decision.
# Keep all seven RP libraries and continue the agreed annotation/filtering,
# TMM and sample-profile stages. No new expression analysis, sample exclusion,
# reanalysis of glaucoma, or workspace checkpoint update occurred here.

# STEP 10: Preserve and check gene identifiers and supplied annotations.
# RP already has separate gene_id and gene_name fields, unlike glaucoma's
# combined labels. Preserve all nine annotation columns without splitting,
# renaming, updating from a database, or collapsing by gene symbol.
rp_gene_annotation <- as.data.frame(
  rp_raw[, c(
    "gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
    "gene_strand", "gene_length", "gene_biotype",
    "gene_description", "tf_family"
  )]
)

# Check dimensions and inspect five examples.
dim(rp_gene_annotation)
head(rp_gene_annotation[, c("gene_id", "gene_name", "gene_biotype")], 5)

# Count missing values and empty/whitespace-only entries in each field.
colSums(is.na(rp_gene_annotation))
sapply(rp_gene_annotation, function(x) {
  sum(trimws(as.character(x)) == "", na.rm = TRUE)
})

# Count duplicated gene IDs (not duplicated gene symbols).
sum(duplicated(rp_gene_annotation$gene_id))

# Check complete mouse Ensembl gene-ID structure, allowing a version suffix.
rp_gene_id_format_ok <- !is.na(rp_gene_annotation$gene_id) & grepl(
  "^ENSMUSG[0-9]{11}(\\.[0-9]+)?$",
  rp_gene_annotation$gene_id
)
table(Complete_ID_format = rp_gene_id_format_ok)

# Confirm annotations follow the count matrix's gene-row order exactly.
identical(
  rp_gene_annotation$gene_id,
  rownames(rp_counts_all)
)

# Expected dimensions: 54532 x 10. Earlier input checks found no missing or
# duplicated IDs; newly reported annotation completeness/format await review.
# Desired ID duplicate count: 0; format TRUE for all; row alignment TRUE.
# A missing gene_name or other annotation does not automatically invalidate
# a count row. Symbols need not be unique; gene_id remains the row key.
# A literal "-" is neither NA nor an empty string (e.g. tf_family).
# This is a structural check of supplied labels, not verification that every
# ID is present/current in a particular Ensembl release or cross-study mapping.
# Existing distribution plots use numeric counts and do not require this
# annotation table, so their validity does not depend on performing this step
# first. Review output before bundling counts/metadata/annotations for QC.
# No data were loaded or analysed by the assistant. No filtering, TMM, MDS,
# PCA, sample exclusion, checkpoint save or changes to existing inputs occur.

# STEP 10 RESULT — 29 October 2026:
# User screenshot confirms 54532 x 10 gene annotation table; NA counts are
# zero in every field except gene_description (43). Blank/whitespace-only
# counts are zero in all fields. Duplicate gene IDs: 0. Complete ID format:
# TRUE for all 54532 rows. Annotation/count-row identity: TRUE.
# rp_gene_id_format_ok is a logical result per gene, not the ID vector.
# Retain the 43 rows with missing descriptions; free-text descriptions are
# not used to calculate filtering, normalization or MDS/PCA coordinates.
# No database membership or gene-symbol uniqueness conclusion is implied.

# STEP 11: Assemble counts, matched sample information and gene annotations.
# A DGEList stores these as separate components in one coordinated object.
# Raw counts are supplied, not the initial log-CPM matrix.
stopifnot(
  identical(rp_metadata_matched$sample_id, colnames(rp_counts_all)),
  identical(rp_gene_annotation$gene_id, rownames(rp_counts_all))
)

# Record genotype categories; preserve the full genotype text in the metadata.
rp_qc_group <- factor(
  rp_metadata_matched$genotype,
  levels = c("RhoWT (wild-type)", "RhoP23H/WT (mutant)"),
  labels = c("WT", "P23H")
)
stopifnot(!anyNA(rp_qc_group))
table(rp_qc_group)

rp_dge <- edgeR::DGEList(
  counts = rp_counts_all,
  samples = rp_metadata_matched,
  group = rp_qc_group,
  genes = rp_gene_annotation
)

# Inspect the size and sample information.
dim(rp_dge)
rp_dge$samples[, c("sample_id", "group", "sex", "lib.size", "norm.factors")]

# Confirm both sample and gene alignment inside the new object.
c(
  sample_order = identical(rp_dge$samples$sample_id, colnames(rp_dge)),
  gene_order = identical(rp_dge$genes$gene_id, rownames(rp_dge))
)

# Expected: WT = 3, P23H = 4; dim = 54532 x 7; both order checks TRUE.
# lib.size contains each sample's raw count total; norm.factors start at 1.
# These initial ones mean no additional composition scaling has been applied.
# DGEList default remove.zeros=FALSE retains all rows. Missing descriptive
# gene text does not remove rows. Samples and annotations remain separate
# from the numerical count matrix; sex and all original metadata are retained.
# group labels identify categories; this does not fit a genotype model,
# average replicates, adjust for sex, or force subsequent plot separation.
# Review the author's output before expression filtering and TMM normalization.
# Assistant has inspected constructor source/documentation but has not
# executed this block, read expression data, or updated the checkpoint.

# STEP 11 RESULT — 29 October 2026:
# User screenshot confirms WT = 3 and P23H = 4; rp_dge is 54532 x 7.
# Sample IDs, genotype groups, sex and library totals match the supplied
# metadata and previous summaries. All norm.factors are initially 1.
# Both sample_order and gene_order checks return TRUE.

# STEP 12: Filter genes with insufficient expression for exploratory QC.
# Genotype group sizes define the required sample support (minimum size 3).
# Use standard filterByExpr count thresholds explicitly for reproducibility.
rp_keep_genes <- edgeR::filterByExpr(
  rp_dge,
  group = rp_qc_group,
  min.count = 10,
  min.total.count = 15
)

# TRUE = keep this gene row; FALSE = leave it out of the QC copy.
table(Retain_for_QC = rp_keep_genes)

# Keep the original rp_dge and create a separate filtered QC object.
# Blank column selection retains all seven samples.
rp_dge_qc <- rp_dge[
  rp_keep_genes,
  ,
  keep.lib.sizes = FALSE
]

dim(rp_dge_qc)

# Confirm sample order and annotation alignment after filtering.
c(
  sample_order = identical(
    rp_dge_qc$samples$sample_id,
    rp_metadata_matched$sample_id
  ),
  gene_order = identical(
    rp_dge_qc$genes$gene_id,
    rownames(rp_dge_qc)
  )
)

# Expected: an as-yet-unknown number of retained gene rows x 7 columns;
# both alignment checks TRUE. Review the author's actual retained/excluded
# counts before TMM. No retained-gene total has been computed by the assistant.
# From reported totals, the median original library size is 17401717.
# CPM cutoff = 10 / 17401717 * 1e6, approximately 0.574656, in at least
# 3 of the 7 libraries; at least 15 counts in total are also required.
# Qualifying libraries need not be the same genotype. Group sizes determine
# support, not a requirement to have differential expression between groups.
# Gene descriptions and symbols do not determine retention.
# keep.lib.sizes=FALSE recalculates totals over retained rows; retained count
# entries are unchanged. All seven samples and corresponding metadata remain.
# This genotype-group filter is for the current exploratory QC preparation;
# it does not fit or fix a later sex-adjusted differential-expression model.
# Source and defaults checked in installed edgeR and official documentation.
# The assistant has not executed this block or changed the saved checkpoint.

# STEP 12 RESULT / THRESHOLD EXPLANATION — 29 October 2026:
# User screenshot confirms FALSE = 37181, TRUE = 17351; rp_dge_qc is
# 17351 x 7, with sample_order and gene_order both TRUE. Retained proportion
# is about 31.8% (arithmetic from the supplied row counts). All seven samples
# remain; excluded gene rows remain in the original rp_dge/rp_counts_all.
# The author asks where the filtering thresholds originate. min.count=10 and
# min.total.count=15 are documented edgeR defaults, supplied explicitly here.
# The approximately 0.575 CPM cutoff is calculated from min.count and the
# original median library size; it is not a fixed CPM default for all studies.
# Required support of 3 comes from min(WT=3, P23H=4). The counts of retained
# or excluded genes did not determine these pre-specified settings.
# large.n=10 and min.prop=0.7 are further defaults affecting large groups;
# they do not reduce the sample requirement for these 3/4 genotype groups.
# Documentation can be opened with help("filterByExpr", package="edgeR").
# Passing the filter means sufficient support under this rule, not proof of
# differential expression; failing does not prove a gene is biologically absent.
# Group labels supply group sizes, not preferential retention of genes that
# separate WT and P23H. No threshold or workflow change is requested by this
# methodological question. Next remains TMM, not executed or supplied yet.

# STEP 13: TMM composition scaling and normalized log-CPM for exploration.
# The author requests MDS; filtering is complete but TMM is not yet confirmed.
# Complete this step before calculating MDS coordinates.
rp_dge_qc <- edgeR::normLibSizes(rp_dge_qc, method = "TMM")

rp_dge_qc$samples[
  , c("sample_id", "group", "lib.size", "norm.factors")
]
stopifnot(
  all(is.finite(rp_dge_qc$samples$norm.factors)),
  all(rp_dge_qc$samples$norm.factors > 0)
)

rp_log_cpm_qc <- edgeR::cpm(
  rp_dge_qc,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)
stopifnot(
  all(is.finite(rp_log_cpm_qc)),
  identical(colnames(rp_log_cpm_qc), rp_dge_qc$samples$sample_id)
)

# TMM stores one factor per sample; effective library size is the retained
# library total times that factor. Counts in rp_dge_qc remain raw counts.
# A neutral factor is 1, not 0. Factor magnitudes are not quality scores.
# Initial all-gene rp_log_cpm remains separate from this filtered/TMM matrix.
# TMM assumes a sufficient reference set without systematic directional change;
# it is not a method to remove sex or force genotype clusters to separate.

# STEP 14: Pairwise-gene MDS, calculated once and reused for display/export.
# Same top=500, gene.selection="pairwise" approach as the glaucoma precheck.
# Correct exported generic is limma::plotMDS. A log-CPM matrix is supplied.
rp_mds <- limma::plotMDS(
  rp_log_cpm_qc,
  top = 500,
  gene.selection = "pairwise",
  plot = FALSE
)

rp_mds_palette <- c(WT = "#0072B2", P23H = "#D55E00")
rp_mds_colours <- rp_mds_palette[as.character(rp_dge_qc$samples$group)]
rp_mds_labels <- paste(
  rp_dge_qc$samples$sample_id,
  substr(rp_dge_qc$samples$sex, 1, 1)
)

rp_draw_mds <- function() {
  limma::plotMDS(
    rp_mds,
    labels = rp_mds_labels,
    col = rp_mds_colours,
    cex = 0.8,
    main = "RP retinal libraries: MDS after TMM"
  )
  legend(
    "topright",
    legend = names(rp_mds_palette),
    text.col = rp_mds_palette,
    bty = "n", cex = 0.85
  )
}

# Display the MDS plot.
rp_draw_mds()

# Save the same plot.
rp_mds_file <- file.path(rp_output_dir, "RP_GSE156533_TMM_MDS.png")
png(rp_mds_file, width = 2400, height = 1800, res = 300)
rp_draw_mds()
dev.off()
file.exists(rp_mds_file)

# Full target:
# results/retinitis_pigmentosa/RP_GSE156533_TMM_MDS.png
# Expected: seven sample labels, F=Female/M=Male; blue WT/orange P23H;
# a saved PNG and file.exists TRUE. Actual factors, axes and layout are pending.
# Each of the 21 sample pairs uses its own top 500 most divergent retained
# genes to calculate a distance; MDS summarizes those distances in two axes.
# The 500-gene choice does not reduce rp_dge_qc below 17351 rows.
# Colors, genotype and sex labels do not determine the MDS coordinates.
# Reusing the stored MDS object avoids recomputation between screen and PNG.
# Rerunning replaces this figure. The workspace checkpoint is not updated
# by these plot-save commands; refresh it after the user confirms the result.
# Installed method signatures/source and official edgeR/limma docs checked.
# No user expression calculations, figure creation or checkpoint writes were
# executed by the assistant. PCA remains a later separate plot.

# STEP 14 RESULT / INTERPRETATION — 29 October 2026:
# The author supplies the MDS figure. The saved PNG exists (169436 bytes,
# 29 October 12:49 timestamp). MDS creation/export are confirmed.
# Dimension 1 is labelled 39%, dimension 2 22% (about 61% combined in this
# MDS representation). These are not causal percentages assigned to genotype
# and sex, nor a percentage for every gene's expression variation.
# All three WT libraries are left of the four P23H libraries on dimension 1;
# both sexes follow that genotype-associated separation. Females are higher
# than males on dimension 2 in both genotypes. WT_12 is the sole male WT;
# its lower position parallels the two male P23H libraries rather than
# establishing a technical failure. WT female variation is visible; no
# isolated unexplained sample outlier is established from this view.
# Retain all seven libraries. The pattern supports genotype comparisons with
# sex retained for later modelling, not gene-level differential-expression
# calls, a causal allocation of the axes, or complete QC certification.
# Top-right female P23H labels overlap each other/the legend. This is a
# display limitation; improve layout before a presentation/export deliverable.
# Actual TMM factors were not shown, so their magnitudes remain unreviewed.
# No coordinates, normalization or images were recomputed by the assistant.
# PCA is still pending. Existing R checkpoint has not yet been refreshed.

# STEP 15: PCA after filtering and TMM.
# The author clarified that the preceding request meant PCA, not another
# density plot. This replaces the unexecuted post-TMM density block.
# Use the same common top-500 PCA approach as the glaucoma precheck.
# limma's gene.selection="common" mode selects one common set of genes
# with greatest variance across samples and constructs PCA.
rp_pca <- limma::plotMDS(
  rp_log_cpm_qc,
  top = 500,
  gene.selection = "common",
  plot = FALSE
)

rp_pca_palette <- c(WT = "#0072B2", P23H = "#D55E00")
rp_pca_colours <- rp_pca_palette[
  as.character(rp_dge_qc$samples$group)
]
rp_pca_labels <- paste(
  rp_dge_qc$samples$sample_id,
  substr(rp_dge_qc$samples$sex, 1, 1)
)

# Redraw existing PCA scores with points, separated labels and an outside legend.
rp_pca_plot_data <- data.frame(
  PC1 = rp_pca$x,
  PC2 = rp_pca$y,
  group = rp_dge_qc$samples$group,
  label = rp_pca_labels
)

rp_pca_plot <- ggplot2::ggplot(
  rp_pca_plot_data,
  ggplot2::aes(x = PC1, y = PC2, colour = group)
) +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    size = 3.5,
    seed = 1,
    box.padding = 0.6,
    point.padding = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 2,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = rp_pca_palette) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = 0.2)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = 0.2)
  ) +
  ggplot2::labs(
    title = "RP retinal libraries: PCA after TMM",
    x = sprintf("PC1 (%.0f%%)", 100 * rp_pca$var.explained[1]),
    y = sprintf("PC2 (%.0f%%)", 100 * rp_pca$var.explained[2]),
    colour = "Genotype"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

rp_draw_pca <- function() {
  print(rp_pca_plot)
}

# Display.
rp_draw_pca()

# Save the same coordinates and plot settings.
rp_pca_file <- file.path(
  rp_output_dir,
  "RP_GSE156533_TMM_PCA_top500.png"
)
png(rp_pca_file, width = 2400, height = 1800, res = 300)
rp_draw_pca()
dev.off()
file.exists(rp_pca_file)
# Expected: TRUE.
# Full target:
# results/retinitis_pigmentosa/RP_GSE156533_TMM_PCA_top500.png
# The input uses all retained genes to select the same top 500 variable genes
# for every sample. Coordinates use gene centering without unit-variance
# gene scaling. This matches the documented common-gene PCA mode used for
# glaucoma. In contrast, the earlier pairwise MDS can use different genes for
# each sample pair. PCA and MDS are complementary views of the same libraries.
# The genotype legend is outside the panel on the right. Points mark true
# sample scores and ggrepel separates labels with connector lines. Confirm
# the exported label layout from the author's revised plot.
# Object class is MDS, with x/y coordinate vectors, despite this PCA mode.
# Do not index rp_pca$x as if it were a prcomp score matrix.
# Official limma documentation and installed source were checked. The
# assistant has not run PCA, loaded expression data, created a PNG or
# refreshed the R checkpoint. Review user output, then save a new checkpoint.

# STEP 15 PCA RECEIVED / LAYOUT REVISION — 29 October 2026:
# User's PCA image shows PC1 48%, PC2 24%, seven samples. Genotype remains
# separated horizontally; males are above females on PC2, the reverse of
# the displayed MDS vertical orientation. Full PCA interpretation is pending.
# PCA PNG exists (177279 bytes, 29 October 13:30). The author asks what the
# top-centre WT text means; clarified that WT/P23H are legend entries.
# User then requests separation of the overlapping P23H_6 F / P23H_8 F
# labels and an outside colour legend. Updated the drawing code above to
# use existing PCA scores with ggplot2/ggrepel (installed 4.0.3 / 0.9.8).
# Coordinates and normalization are not recomputed. Actual points are shown;
# labels can move with segments back to the points. Right-side outside
# legend is titled Genotype and contains point keys, not floating text.
# The helper/export block replaces the existing PCA PNG when the user runs
# it. The assistant has not executed the revised drawing or overwritten the
# PNG. Revised visual QA, detailed interpretation and checkpoint refresh
# remain pending. Official ggrepel examples support these layout settings.

# APPROVED SHARED FIGURE STYLE — 29 October 2026:
# User confirms the revised RP PCA appearance and requests the same style
# for RP MDS and glaucoma MDS/PCA. The existing PCA is the visual reference.
# Standalone code reuses stored scores, uses points/repelled sex-labelled
# sample names, and places colour legends outside on the right.
# Run manually from the current RP session:
# source(file.path("scripts", "04_Restyle_Retinal_MDS_PCA.R"))
# It exports all three requested PNGs and replaces the drawing helpers in
# the active session. Glaucoma remains closed; this is figure formatting.
# User execution, revised three-figure QA and checkpoint refresh are pending.

# PCA/MDS COMPARISON AND STYLE CONFIRMATION — 29 October 2026:
# User reports the three restyled exports are done and requests RP PCA
# interpretation, noting similar sample relationships to MDS. Opened all
# four final PNGs for visual comparison/layout QA, without reading matrices
# or rerunning projections. Points, labels and outside legends are readable.
# RP PCA: WT left/P23H right on PC1; males above/females below on PC2.
# P23H_6 F/P23H_8 F form a particularly close pair; P23H_11 M/P23H_13 M
# are also near. WT_7 F/WT_9 F remain in the same neighbourhood; WT_12 M
# is separate vertically in the same sex-associated pattern as P23H males.
# PC1 48% + PC2 24% = approximately 72% of variation for the selected
# common 500 genes' log-expression values, not all 17351 retained genes.
# The main genotype/sex associations persist under common versus pairwise
# gene selection. These are related summaries of the same data, not
# independent validation. Exact pair distances and axis percentages are
# not interchangeable between methods; do not infer numerical rankings
# from visual inspection. The opposite vertical orientation is not a
# biological reversal: eigenvector/component signs are arbitrary.
# Keep all seven libraries and retain sex for later genotype modelling.
# No gene-level significance, pathways, causal variance assignment or
# resemblance to spaceflight is inferred from these exploratory plots.
# Actual TMM factor table and refreshed R checkpoint are still outstanding.

# PRECHECK CLOSED AT USER'S CHOSEN ENDPOINT — 29 October 2026:
# Retain all seven libraries. The QC copy has 17351 retained genes.
# Genotype and sex structure recur in MDS/PCA; retain sex for later modelling.
# Completed summary:
# docs/retinitis-pigmentosa-precheck.md
# Save the current R objects below before leaving this session. This is a
# separate dated checkpoint; keep the existing initial-stage file as backup.

rp_tmm_summary <- rp_dge_qc$samples[
  , c("sample_id", "group", "sex", "lib.size", "norm.factors")
]
rp_checkpoint_file <- file.path(
  rp_output_dir,
  "RP_GSE156533_precheck_checkpoint_2026-10-29.RData"
)

# Repository checkpoint bookkeeping: the interactive recovery helper creates
# this base list. Initialize it here too for a fresh block-by-block run.
# Run the shared styling block before saving the final styled plot objects.
if (!exists("rp_checkpoint_objects", inherits = FALSE)) {
  rp_sample_ids <- colnames(rp_counts_all)
  rp_checkpoint_objects <- c(
    "rp_output_dir", "rp_dir", "rp_count_file",
    "rp_metadata_archive", "rp_metadata_dir", "rp_metadata_file",
    "rp_sheets", "rp_sample_ids", "rp_raw",
    "rp_counts_all", "rp_metadata_all", "rp_metadata_matched",
    "rp_characteristic_tags", "rp_treatment_protocol", "rp_extraction_protocol",
    "rp_total_counts", "rp_detected_genes", "rp_log_cpm",
    "rp_draw_boxplot", "rp_boxplot_file", "rp_checkpoint_file",
    "rp_checkpoint_objects"
  )
}

rp_checkpoint_objects <- unique(c(
  rp_checkpoint_objects,
  "rp_gene_annotation", "rp_gene_id_format_ok", "rp_qc_group",
  "rp_dge", "rp_keep_genes", "rp_dge_qc", "rp_log_cpm_qc", "rp_tmm_summary",
  "rp_density_curves", "rp_density_colors", "rp_draw_density", "rp_density_file",
  "rp_mds", "rp_mds_palette", "rp_mds_colours", "rp_mds_labels",
  "rp_mds_plot", "rp_draw_mds", "rp_mds_file",
  "rp_pca", "rp_pca_palette", "rp_pca_colours", "rp_pca_labels",
  "rp_pca_plot_data", "rp_pca_plot", "rp_draw_pca", "rp_pca_file"
))

save(list = rp_checkpoint_objects, file = rp_checkpoint_file)
file.exists(rp_checkpoint_file)
# Expected: TRUE. User execution and save confirmation are pending.
# XML document/node pointers and the isolated glaucoma environment are not
# in this explicit list. The current RP plotting helpers use saved plot objects.
# Numerical TMM metadata are preserved in rp_tmm_summary and rp_dge_qc.
# Future restart after successful save:
# load(file.path("results", "retinitis_pigmentosa", "RP_GSE156533_precheck_checkpoint_2026-10-29.RData"))


# FINAL CHECKPOINT / HANDOFF VERIFIED — 29 October 2026:
# Final dated checkpoint now exists: 9637922 bytes, 29 October 14:06.
# An isolated R-process read confirmed 49 saved RP objects; original count
# dimensions 54532 x 7 and QC dimensions 17351 x 7; projection modes
# pairwise (MDS) and common (PCA). The existing rp_tmm_summary was read,
# without recalculating normalization: factors range 0.9409120709-1.0847612352.
# Metadata/projection inventory was checked; no expression analysis rerun.
# Detailed handoff with full code, results, explanations and next-chat prompt:
# docs/retinitis-pigmentosa-precheck.md (portable completion record)
# Archive copy is in this RP output folder. RP remains closed with all seven
# libraries retained. Next agreed panel dataset: DR GSE87433 (microarray).
