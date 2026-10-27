# Glaucoma GSE241782: pre-quality checks
# Started: 26 October 2026
# Run one small block at a time and review its output before proceeding.
# Portable repository copy; original interactive working script is retained.
# Start R in the repository root. See docs/running-the-precheck.md.
# Set GLAUCOMA_DATA_DIR in R to use inputs stored outside this repository.
# Recorded observations below come from the author's interactive runs.
# Repository preparation checked syntax only; it did not rerun the analysis.

# STEP 1: Confirm that R can find the dataset folder and both input files.
# The dataset folder has been renamed to "Glaucoma_GSE241782".
# The downloaded files are already uncompressed: .csv and .soft.

# <- assigns the folder path to an R object named glaucoma_dir.
glaucoma_dir <- Sys.getenv(
  "GLAUCOMA_DATA_DIR",
  unset = file.path("data", "raw", "Glaucoma_GSE241782")
)

# file.path() joins the folder path and filename to locate each input.
# glaucoma_count_file points to the gene-expression measurements.
glaucoma_count_file <- file.path(glaucoma_dir, "GSE241782_gene_count_matrix.csv")

# glaucoma_metadata_file points to the GEO study and sample descriptions.
glaucoma_metadata_file <- file.path(glaucoma_dir, "GSE241782_family.soft")

# dir.exists() asks whether R can find the dataset folder.
dir.exists(glaucoma_dir)

# file.exists() asks whether R can find each input file.
file.exists(glaucoma_count_file)
file.exists(glaucoma_metadata_file)

# Expected output: TRUE for each of the three checks above.
# These checks confirm file locations; they do not establish data quality.
# The author confirmed that all three checks returned TRUE on 26 October 2026.

# STEP 2: Read the gene-count table and inspect its structure.
# Run this block only after completing Step 1 in the same R session.
# read.csv() imports a comma-separated file as an R data frame (a table).
# check.names = FALSE preserves sample names such as "R-01" exactly.
# stringsAsFactors = FALSE keeps text columns as text.
glaucoma_raw <- read.csv(
  glaucoma_count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# dim() reports the number of rows followed by the number of columns.
# At this point, the gene_id column is included in the column count.
dim(glaucoma_raw)

# colnames() displays the column labels, including gene_id and sample names.
colnames(glaucoma_raw)

# In [rows, columns], a blank row position selects all rows.
# 1:6 selects columns 1 through 6; head(..., n = 5) shows their first 5 rows.
head(glaucoma_raw[, 1:6], n = 5)

# This full table contains multiple eye-tissue regions.
# We will identify retinal samples using metadata before selecting them.
# Review these three outputs before checking numeric values or filtering.

# STEP 2 RESULT (recorded interactive R output, 26 October 2026):
# 55,487 rows and 37 columns: gene_id plus 36 sample columns.
# The labels contain 12 MON, 12 R and 12 UON sample columns.
# The preview shows identifiers containing an Ensembl ID and gene symbol,
# separated by |. Keep the original identifiers during these input checks.

# STEP 3: Check missing values, duplicate gene IDs and numeric sample columns.
# anyNA() returns TRUE if any cell contains an R missing value (NA).
# Expected: FALSE. This does not check for every possible invalid value.
anyNA(glaucoma_raw)

# $ selects a named column. anyDuplicated() returns 0 when no duplicates
# exist, or the row position of the first duplicated identifier otherwise.
# This checks the complete identifiers exactly as imported.
# Expected: 0.
anyDuplicated(glaucoma_raw$gene_id)

# [, -1] selects all rows and excludes the first column (gene_id).
# sapply() applies is.numeric() to each of the remaining sample columns.
# is.numeric() checks whether R stored each column as numeric data.
# Expected: TRUE for all 36 sample columns.
sapply(glaucoma_raw[, -1], is.numeric)

# Review the results before proceeding to count-value and sample-design checks.

# STEP 3 RESULT (recorded interactive R output, 26 October 2026):
# anyNA(): FALSE; anyDuplicated(gene_id): 0;
# all 36 sample columns returned TRUE for is.numeric().

# STEP 4: Prepare a numeric count matrix and inspect the count values.
# as.matrix() converts the numeric sample columns into a rectangular matrix.
# The suffix _all indicates that all 36 samples are still included.
# glaucoma_raw remains available as the original imported table.
glaucoma_counts_all <- as.matrix(glaucoma_raw[, -1])

# rownames() attaches each original gene identifier to its count-matrix row.
# The gene identifiers become row labels, not numeric measurements.
rownames(glaucoma_counts_all) <- glaucoma_raw$gene_id

# is.finite() checks each value for an ordinary finite number.
# all() returns TRUE only when every value passes.
# Expected: TRUE. This excludes Inf, -Inf, NA and NaN values.
all(is.finite(glaucoma_counts_all))

# range() returns the smallest and largest values in the entire matrix.
# Check that the minimum is at least zero; zeros are valid count values.
range(glaucoma_counts_all)

# round() supplies whole-number values for comparison with ==.
# TRUE means every count equals its nearest whole number.
# This comparison does not assign rounded values back to the count matrix.
# A FALSE result requires reviewing the count-generation method.
all(glaucoma_counts_all == round(glaucoma_counts_all))

# Review these outputs before examining sample groups or sample-level plots.

# STEP 4 RESULT (recorded interactive R output, 26 October 2026):
# All values finite: TRUE; count range: 0 to 6,794,430;
# all counts equal their rounded values: TRUE.
# Basic count-file integrity checks have passed. Sample-level quality and
# experimental-design checks remain pending.

# STEP 5A: Read and inspect the local GEO metadata before matching samples.
# A SOFT file is structured text rather than a CSV table.
# readLines() creates a character vector with one element per file line.
glaucoma_soft_lines <- readLines(glaucoma_metadata_file)

# grep() finds lines matching a pattern.
# ^ means the beginning of a line; value = TRUE returns the matching text.
# Sample-title lines describe tissue, time, treatment and sex.
glaucoma_title_lines <- grep(
  pattern = "^!Sample_title",
  x = glaucoma_soft_lines,
  value = TRUE
)

# Supplementary-file lines include filenames with labels such as UON-01.
# These labels provide a link to the count-matrix column names.
glaucoma_supplement_lines <- grep(
  pattern = "^!Sample_supplementary_file",
  x = glaucoma_soft_lines,
  value = TRUE
)

# length() counts the sample-title entries. Expected: 36.
length(glaucoma_title_lines)

# Inspect the first three entries of each kind before building a sample table.
head(glaucoma_title_lines, n = 3)
head(glaucoma_supplement_lines, n = 3)

# Do not assign treatments from the count-column order. We will match
# explicit sample identifiers and confirm one-to-one correspondence next.

# STEP 5A RESULT (recorded interactive R output, 26 October 2026):
# 36 sample-title lines were found. The first three titles and supplementary
# filenames were displayed successfully after correcting the grep pattern.

# STEP 5B: Clean the field labels and obtain supplementary filenames.
# sub() replaces the first occurrence of a pattern in each text entry.
# replacement = "" removes that text; fixed = TRUE treats the pattern literally.
glaucoma_titles <- sub(
  pattern = "!Sample_title = ",
  replacement = "",
  x = glaucoma_title_lines,
  fixed = TRUE
)

# This file uses the numbered field !Sample_supplementary_file_1 for every
# sample. Remove its label to retain only the file URL as text.
glaucoma_file_urls <- sub(
  pattern = "!Sample_supplementary_file_1 = ",
  replacement = "",
  x = glaucoma_supplement_lines,
  fixed = TRUE
)

# basename() returns the filename after the final slash in each URL.
# This is a text operation; it does not download the files.
glaucoma_file_names <- basename(glaucoma_file_urls)

# Inspect the cleaned entries and confirm the number of supplementary filenames.
head(glaucoma_titles, n = 3)
head(glaucoma_file_names, n = 3)
length(glaucoma_file_names)

# Expected: three clean titles, three GSM..._UON-..._genes.tsv.gz filenames,
# and a filename count of 36. Review before extracting sample codes.

# STEP 5B RESULT (recorded interactive R output, 26 October 2026):
# Cleaned sample titles were displayed successfully. The filename preview
# contains GSM7743510_UON-01_genes.tsv.gz and the following sample filenames;
# length(glaucoma_file_names) returned 36. Identifier matching is still pending.

# STEP 5C: Extract matrix sample identifiers and check identifier coverage.
# This block uses glaucoma_file_names created in Step 5B.
# Example: GSM7743510_UON-01_genes.tsv.gz -> UON-01.
# ^GSM[0-9]+_ matches GSM followed by one or more digits and an underscore
# at the beginning of each filename. The first sub() removes that prefix.
glaucoma_sample_ids <- sub(
  pattern = "^GSM[0-9]+_",
  replacement = "",
  x = glaucoma_file_names
)

# Remove the literal suffix to retain only the matrix sample identifier.
glaucoma_sample_ids <- sub(
  pattern = "_genes.tsv.gz",
  replacement = "",
  x = glaucoma_sample_ids,
  fixed = TRUE
)

# Inspect the first three extracted identifiers.
# Expected: UON-01, UON-02, UON-03.
head(glaucoma_sample_ids, n = 3)

# Confirm the metadata identifier count. Expected: 36.
length(glaucoma_sample_ids)

# Check for duplicate identifiers on both sides of the planned match.
# Expected: 0 for each check.
anyDuplicated(glaucoma_sample_ids)
anyDuplicated(colnames(glaucoma_counts_all))

# setequal() checks whether both sides contain the same identifiers,
# irrespective of order. Expected: TRUE.
setequal(glaucoma_sample_ids, colnames(glaucoma_counts_all))

# These checks establish unique identifier coverage before attaching titles.
# We will verify title-to-file pairing within SOFT sample records, build the
# metadata table, and align its rows to the count-column order next.

# STEP 5C RESULT (recorded interactive R output, 26 October 2026):
# 36 identifiers; 0 duplicate metadata identifiers; 0 duplicate count-column
# names; setequal(metadata identifiers, count-column names): TRUE.

# STEP 5D: Build the metadata table and align its rows to the count columns.
# Source-file inspection confirmed 36 SOFT sample records, each containing
# exactly one title and one supplementary-file entry. The separately extracted
# title and filename vectors therefore share the same SOFT sample-record order.
# That record order differs from the count-matrix column order.

# data.frame() combines the related entries into one row per sample.
# Keep metadata as a data frame because its columns contain descriptions/IDs.
glaucoma_metadata_all <- data.frame(
  sample_id = glaucoma_sample_ids,
  title = glaucoma_titles,
  file_name = glaucoma_file_names,
  stringsAsFactors = FALSE
)

# match(target IDs, source IDs) returns the source-row position for each target.
# Our target order is the count-matrix column order.
glaucoma_metadata_order <- match(
  colnames(glaucoma_counts_all),
  glaucoma_metadata_all$sample_id
)

# In [rows, columns], the positions before the comma reorder the rows.
# Leaving the column position blank retains all metadata columns.
glaucoma_metadata_all <- glaucoma_metadata_all[glaucoma_metadata_order, ]

# Reset the displayed row numbers after reordering.
# The actual identifiers remain in the sample_id column.
rownames(glaucoma_metadata_all) <- NULL

# Expected dimensions: 36 rows and 3 columns.
dim(glaucoma_metadata_all)

# identical() checks both identifier values and their exact order.
# Expected: TRUE.
identical(glaucoma_metadata_all$sample_id, colnames(glaucoma_counts_all))

# Preview the first three matched rows; expected sample IDs: MON-01 to MON-03.
head(glaucoma_metadata_all, n = 3)

# Review the aligned metadata before separating tissue/time/treatment/sex fields
# and choosing the retinal glaucoma/control comparison.

# STEP 5D RESULT (recorded interactive R output, 26 October 2026):
# Metadata dimensions: 36 rows, 3 columns. Exact identifier/order check: TRUE.
# The first rows are MON-01, MON-02 and MON-03 with their sample descriptions.

# STEP 5E: Separate the titles into tissue, time, treatment and sex.
# For this dataset, each of the 36 titles has four comma-separated fields.
# read.csv() can read existing text in R as well as a file on disk.
# text supplies the title strings, one sample per line.
# header = FALSE keeps the first sample as data rather than a header row.
# col.names supplies the four column labels; c() combines them into a vector.
# strip.white = TRUE removes spaces surrounding unquoted field values.
glaucoma_title_fields <- read.csv(
  text = glaucoma_metadata_all$title,
  header = FALSE,
  col.names = c("tissue", "time", "treatment", "sex"),
  strip.white = TRUE,
  stringsAsFactors = FALSE
)

# Expected: 36 rows and 4 columns, in the aligned metadata's row order.
dim(glaucoma_title_fields)

# Expected: FALSE, indicating no missing parsed field values.
anyNA(glaucoma_title_fields)

# Inspect the parsed fields before adding them to the sample metadata.
head(glaucoma_title_fields, n = 3)

# Keep the original title column for provenance. Review these outputs before
# attaching the four fields and checking the available experimental groups.

# STEP 5E RESULT (recorded interactive R output, 26 October 2026):
# Parsed fields: 36 rows, 4 columns, no missing values.
# The first six rows agree with the aligned MON sample titles.

# STEP 5F: Add the four parsed columns to the aligned metadata.
# These rows correspond because the fields were parsed directly from
# glaucoma_metadata_all$title and neither table has since been reordered.
# $ selects a named column. Each assignment copies values from the parsed
# table into that named column in the full metadata, creating it if absent.
# Rerunning these assignments updates the same columns without duplicating them.
glaucoma_metadata_all$tissue <- glaucoma_title_fields$tissue
glaucoma_metadata_all$time <- glaucoma_title_fields$time
glaucoma_metadata_all$treatment <- glaucoma_title_fields$treatment
glaucoma_metadata_all$sex <- glaucoma_title_fields$sex

# Expected: 36 rows and 7 columns (the original 3 plus the new 4).
dim(glaucoma_metadata_all)

# Confirm that the sample IDs still follow the count-matrix column order.
# Expected: TRUE.
identical(
  glaucoma_metadata_all$sample_id,
  colnames(glaucoma_counts_all)
)

# Display just the identifier and four parsed fields for the first six rows.
# This preview does not remove the original title or file_name columns.
head(
  glaucoma_metadata_all[, c("sample_id", "tissue", "time", "treatment", "sex")],
  n = 6
)

# Review these outputs before checking the available tissue/treatment groups.

# STEP 5F OBSERVATION (recorded interactive R output, 26 October 2026):
# The metadata preview contains the original three columns and all four parsed
# columns. Each row still represents one of the 36 imported sample libraries.
# This is the assembled metadata for all tissues, before selecting a retinal
# analysis subset or completing experimental-design and expression QC checks.

# STEP 6: Save the current R workspace and a readable copy of the metadata.
# Run this block in the active R session after Steps 1-5.
# Adding this code to the script does not itself export any R objects.
# These filenames identify an in-progress checkpoint; prechecks are incomplete.

# Use an explicit output folder so saving does not depend on getwd().
glaucoma_precheck_dir <- file.path("results", "glaucoma")
dir.create(glaucoma_precheck_dir, recursive = TRUE, showWarnings = FALSE)

# Confirm the output folder exists and the metadata is still aligned.
# stopifnot() raises an error if a supplied condition is not TRUE.
stopifnot(
  dir.exists(glaucoma_precheck_dir),
  identical(glaucoma_metadata_all$sample_id, colnames(glaucoma_counts_all))
)

# file.path() joins the output folder and each filename.
glaucoma_metadata_csv <- file.path(
  glaucoma_precheck_dir,
  "Glaucoma_GSE241782_metadata_all.csv"
)
glaucoma_workspace_file <- file.path(
  glaucoma_precheck_dir,
  "Glaucoma_GSE241782_precheck_checkpoint.RData"
)

# write.csv() exports the seven-column metadata table for inspection.
# row.names = FALSE omits R's displayed row numbers; sample_id is retained.
# UTF-8 preserves accented treatment labels present in the original titles.
write.csv(
  glaucoma_metadata_all,
  file = glaucoma_metadata_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# save.image() saves all objects currently in the global R workspace, including
# glaucoma_raw, glaucoma_counts_all and glaucoma_metadata_all. It preserves the
# R object types and matrix row/column names. The code itself is kept in this
# separate .R script. Rerunning the saves replaces these checkpoint files.
save.image(file = glaucoma_workspace_file)

# Expected: TRUE TRUE. These confirm that both output files exist.
file.exists(c(glaucoma_metadata_csv, glaucoma_workspace_file))

# For a later fresh R session, restore this workspace with the command below.
# Do not run it now; continue using the objects already present in this session.
# load(file.path("results", "glaucoma", "Glaucoma_GSE241782_precheck_checkpoint.RData"))

# NEXT: Review tissue/treatment/time/sex groups and independent replication,
# select the retinal glaucoma/control samples, and perform sample-level QC.
# No normalized expression, differential expression or disease ranking has
# been calculated by this script at the checkpoint above.

# STEP 7A: Experimental-group check: tissue and treatment inventory.
# Added: 27 October 2026. Run in the restored R session.
# The checkpoint export paths above now match the actual glaucoma subfolder.
# This block counts existing metadata rows; it does not remove samples or
# change the counts or metadata. Each row represents a sequenced library.
# Numbers of libraries do not establish numbers of independent animals;
# pooling and biological replication must be checked against study records.

# table() counts how often each tissue label occurs.
# useNA = "ifany" displays a missing-value category only if NA values exist.
# Expected tissue totals: MON 12, Retina 12, UON 12.
table(glaucoma_metadata_all$tissue, useNA = "ifany")

# Cross-tabulate tissue and treatment: rows are tissues, columns are treatments.
# Tissue and Treatment here label the two dimensions in the printed table.
# Each cell is the number of libraries with that tissue/treatment combination.
table(
  Tissue = glaucoma_metadata_all$tissue,
  Treatment = glaucoma_metadata_all$treatment,
  useNA = "ifany"
)

# Await the author's two outputs before interpreting the treatment counts.
# Next: inspect collection times and sex within the observed groups, then
# verify controls and independent replication/pooling in the study records.

# STEP 7A RESULT (recorded interactive R output, 27 October 2026):
# Tissue totals: MON 12, Retina 12, UON 12 (36 libraries overall).
# Within EACH tissue: Crush 4, Glaucoma 6, Naive (accented label in data) 2.
# The primary tissue match to the fixed OSD-758 retinal reference is Retina.
# The candidate retinal glaucoma comparison therefore has six glaucoma and
# two naive-control libraries; four retinal crush libraries are outside that
# comparison. No subset has yet been created by this block. Collection-time,
# sex and independent-replication checks remain pending.

# STEP 7B: Inspect collection time and sex within the retinal libraries.
# Added: 27 October 2026. Run this block in the active R session.
# Select rows whose tissue is exactly "Retina". == compares values; <- assigns.
# In [rows, columns], leaving columns blank retains every metadata column.
# This new inspection table contains all 12 retinal libraries, including crush.
# The all-tissue master metadata and count matrix are retained unchanged.
glaucoma_metadata_retina <- glaucoma_metadata_all[
  glaucoma_metadata_all$tissue == "Retina",
]

# Expected: 12 rows and 7 columns.
dim(glaucoma_metadata_retina)

# Three-variable table: treatments are rows, collection times are columns,
# and separate treatment-by-time tables are printed for each sex.
# Each cell counts libraries with that treatment, collection time and sex.
# Here zero means no library in that combination, not zero gene expression.
# useNA = "ifany" displays missing categories only if NA values are present.
table(
  Treatment = glaucoma_metadata_retina$treatment,
  Time = glaucoma_metadata_retina$time,
  Sex = glaucoma_metadata_retina$sex,
  useNA = "ifany"
)

# Await the author's dimensions and sex-specific tables before interpretation.
# No retinal count subset or glaucoma/control analysis subset is created here.
# Next: verify source-record pooling/independence and controls, then derive
# the matched count/metadata subset while preserving collection-time labels.

# STEP 7B RESULT (recorded interactive R output, 27 October 2026):
# glaucoma_metadata_retina has 12 rows and 7 columns.
# Each listed retinal treatment/time combination has one female library and
# one male library: Naive 0D; Glaucoma 3D, 2W, 6W; Crush 3D, 2W.
# Thus the candidate glaucoma/control set contains eight libraries in four
# groups: two baseline control libraries and two glaucoma libraries per time.
# The six glaucoma libraries span three collection times, not six replicates
# collected at a single time. The metadata has one shared Naive 0D group.

# STEP 7C: Display the study's pooling statements from the original SOFT text.
# Source verification: Keuthan et al. (2023), Methods 4.8 and 4.9,
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10531004/
# and the GEO design/extraction protocol in the local family SOFT file.
# The study reports four original RNA samples per tissue/treatment/time,
# pooled within sex into two RNA-seq libraries (one male, one female).
# A count-matrix column is one pooled library; the original input samples
# cannot be treated as separately observed expression profiles.
# This limits precision and must remain explicit in the eventual analysis.
# This is a study-design finding, not a decision to discard the dataset.

# grep() retrieves lines containing pool or pooled, regardless of letter case.
# unique() displays each distinct statement once instead of repeating the
# same extraction protocol for every sample. It does not modify the SOFT text.
unique(
  grep(
    pattern = "pool",
    x = glaucoma_soft_lines,
    value = TRUE,
    ignore.case = TRUE
  )
)

# Review the displayed source statements before building the eight-library
# glaucoma/control subset and its matching count matrix. No count subset,
# gene filtering, normalization or differential expression is run here.

# STEP 7C RESULT (recorded interactive R output, 27 October 2026):
# unique(grep(...)) returned three distinct pooling-related SOFT statements:
# the design reports four original samples pooled into two libraries, RNA
# from sex replicates was pooled in equal amounts, and pooled RNA was sent
# for sequencing-library preparation. The author confirmed this interpretation.
# Expression measurements are available for two pooled libraries per group.

# STEP 8A: Select retinal glaucoma and naive-control metadata for prechecks.
# Added: 27 October 2026. Run this block manually in the active R session.
# The source table already contains only retinal libraries. Keep explicitly
# named treatments, preserving time and sex rather than combining stages.
# c() creates the vector of allowed labels; %in% checks membership in it.
# The original all-tissue and all-retinal metadata objects are retained.
glaucoma_metadata_selected <- glaucoma_metadata_retina[
  glaucoma_metadata_retina$treatment %in% c("Glaucoma", "Naïve"),
]

# Reset displayed row numbers. The sample identifiers stay in sample_id.
rownames(glaucoma_metadata_selected) <- NULL

# Expected: 8 rows and 7 columns.
dim(glaucoma_metadata_selected)

# Expected: Glaucoma 6, Naive 2 (accented label displayed as in the source).
table(glaucoma_metadata_selected$treatment)

# Display all eight selected identifiers, treatments, collection times and sex.
glaucoma_metadata_selected[, c("sample_id", "treatment", "time", "sex")]

# Await these outputs before selecting count columns by the sample_id values.
# The count matrix has not yet been subset, and no expression QC is run here.

# STEP 8A OBSERVATION (recorded interactive R output, 27 October 2026):
# The selected metadata preview shows eight libraries, with treatment, time
# and sex retained: R-03/R-04 (Glaucoma 2W), R-07/R-08 (Glaucoma 3D),
# R-09/R-10 (Glaucoma 6W), R-11/R-12 (Naive 0D).
# Each glaucoma pair is Male/Female; the naive pair is Female/Male.
# This preview displays four fields; the selection code retains all seven.

# STEP 8B: Select the corresponding eight columns from the count matrix.
# Added: 27 October 2026. Run manually after reviewing Step 8A.
# In [rows, columns], blank rows retain every gene. The character sample_id
# vector selects count columns by name and places them in metadata row order.
# drop = FALSE keeps the result as a matrix, including for a one-column subset.
# Store a new matrix; retain the original all-tissue count matrix.
glaucoma_counts_selected <- glaucoma_counts_all[
  ,
  glaucoma_metadata_selected$sample_id,
  drop = FALSE
]

# Expected: 55487 rows (genes) and 8 columns (selected libraries).
# No genes have been filtered out by the sample-selection step.
dim(glaucoma_counts_selected)

# Expected: TRUE. Exact sample identity AND order must agree on both sides.
identical(
  glaucoma_metadata_selected$sample_id,
  colnames(glaucoma_counts_selected)
)

# Preview five genes and the first four selected sample columns.
# This display does not reduce the saved count matrix.
head(glaucoma_counts_selected[, 1:4, drop = FALSE], n = 5)

# Await these outputs before inspecting sample count totals and detected genes.
# No normalization, gene filtering or expression comparison is run here.

# STEP 8B RESULT (recorded interactive R output, 27 October 2026):
# glaucoma_counts_selected dimensions: 55487 genes by 8 selected libraries.
# Exact metadata sample_id / count-column identity and order check: TRUE.
# The preview shows numeric counts for R-03, R-04, R-07 and R-08, with the
# original full gene identifiers retained as row labels.
# The selected count matrix and metadata are now an aligned pair: each
# count-matrix column corresponds to the metadata row at the same position.
# Retain these separate structures; counts stay numeric and sample
# descriptions remain in the metadata. No flat merge is required.
# Next, when continuing the precheck: sample count totals and detected genes.

# STEP 9A: Calculate total gene-level counts per selected library.
# Added: 27 October 2026. Run manually in the active R session.
# Each count-matrix column contains values for the same 55,487 gene rows.
# colSums() adds the values down each column, returning one total per library.
# The resulting named numeric vector retains the column sample IDs as names.
# This is a sum of counts, not a count of nonzero/detected gene rows.
glaucoma_total_counts <- colSums(glaucoma_counts_selected)

# Expected: eight totals, named R-03, R-04, R-07, R-08, R-09, R-10, R-11, R-12.
# Numerical totals must be reviewed from the actual interactive output.
glaucoma_total_counts

# Display the same named totals from smallest to largest for comparison.
# This sorted display does not reorder the count matrix or metadata.
sort(glaucoma_total_counts)

# Await these outputs before interpreting differences or moving to detected
# genes. Library totals alone do not establish whether a sample is usable.
# No sample exclusion, normalization or new export is performed in this block.

# STEP 9A RESULT (recorded interactive R output, 27 October 2026):
# R-03: 40946004; R-04: 32514596; R-07: 46940497; R-08: 19931516;
# R-09: 26597848; R-10: 34211389; R-11: 43795226; R-12: 37153668.
# R-08 (Glaucoma 3D Female) has the lowest total; R-07 (3D Male) the highest.
# The largest total is about 2.36 times the smallest. These are matrix-based
# gene-count totals, not an independent raw-read sequencing QC report.
# Keep all eight libraries at this stage and review other QC summaries;
# having the lowest library total alone does not establish sample failure.

# STEP 9B: Count detected gene rows per selected library.
# Added: 27 October 2026, following the author's request for this next check.
# For this initial summary only, detected means the gene count is > 0.
# This is separate from later expression filtering for statistical analysis.
# > 0 produces a logical matrix: TRUE for positive counts, FALSE for zeros.
# colSums() treats TRUE as 1 and FALSE as 0, so it counts positive rows in
# each sample column, retaining the sample IDs as the vector's names.
glaucoma_detected_genes <- colSums(glaucoma_counts_selected > 0)

# Expected: eight named detected-gene counts, each between 0 and 55487.
# Review the actual interactive values before drawing sample-quality conclusions.
glaucoma_detected_genes

# Display these counts from smallest to largest, retaining sample names.
# This does not reorder or modify the count matrix or metadata.
sort(glaucoma_detected_genes)

# Await the author's results. In particular, inspect whether R-08's lower count
# total is accompanied by fewer detected genes. Differences in library totals
# alone do not establish that a library has fewer detected genes.
# No gene filtering, normalization, sample exclusion or export is run here.

# STEP 9B RESULT (recorded interactive R output, 27 October 2026):
# Detected gene rows (count > 0):
# R-03: 28154; R-04: 27497; R-07: 28191; R-08: 25647;
# R-09: 26320; R-10: 28001; R-11: 29611; R-12: 27855.
# R-08 is lowest and R-11 highest. Range: 25647 to 29611; R-08 is about
# 13.4% below the highest detected-gene count (and 8.2% below median 27928).
# Detected-gene counts vary less than library totals (max/min about 1.15
# versus 2.36). Positive counts of different magnitudes still count as one
# detected gene row, so these two summaries measure different properties.
# Similar numbers detected do not establish that gene identities or relative
# expression profiles are similar. These summaries alone do not justify
# excluding R-08; retain all eight libraries for the remaining QC checks.
# No new executable analysis block is added during this interpretation.
# Next: explain and inspect expression distributions before the later gene
# identifier/annotation and normalized sample-relationship checks.

# STEP 9C: Initial expression-distribution boxplots on the log-CPM scale.
# Added: 27 October 2026, following the author's request to inspect distributions.
# Environment-only check: edgeR 4.8.2 is available in the installed R 4.5.3.
# No expression data or workspace was loaded by that package-availability check.
# Run the following block manually in the active R session.
#
# Rationale: library totals differ by about 2.36-fold. CPM scales gene counts
# per million total counts within each sample; log2 makes their broad range
# easier to compare visually. This is an initial library-size-adjusted view,
# before low-expression filtering or composition-normalization decisions.
# All 55,487 gene rows and all eight libraries remain included.
# The original glaucoma_counts_selected matrix is retained for later analysis.
#
# Official references for this QC component (not validation of the full project):
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
# For a matrix input, cpm() uses column sums as library sizes by default.
# No TMM or other composition-normalization factors are computed by this call.
# edgeR:: calls a function from edgeR without attaching the whole package.
# log = TRUE returns log2-CPM values.
# prior.count = 2 supplies the average prior count, scaled by library size,
# so zero counts can be included in the log calculation without log(0).
glaucoma_log_cpm <- edgeR::cpm(
  glaucoma_counts_selected,
  log = TRUE,
  prior.count = 2
)

# Each matrix column becomes a box labelled with its sample ID.
# las = 2 rotates library labels; col supplies the box fill colour.
# outline = FALSE hides individual gene-outlier dots to keep boxes readable;
# no genes are removed from the calculation or from the count matrix.
# main, xlab and ylab label the plot and its axes.
boxplot(
  glaucoma_log_cpm,
  las = 2,
  col = "lightblue",
  outline = FALSE,
  main = "Retinal libraries: expression distributions",
  xlab = "Library",
  ylab = "Log2 counts per million (CPM)"
)

# Expected: eight labelled boxes in R's graphics window. The line inside a box
# is the median (middle value); the box covers the middle 50% of gene values.
# Negative log-CPM values are valid and do not mean raw counts are negative.
# All gene rows, including zeros, are present. A median at the low-expression
# floor can reflect the proportion of zeros rather than a sample failure.
# Similar boxes do not establish similar gene identities or relative profiles;
# different boxes require interpretation alongside library totals and design.
# Await the author's plot before interpreting or adding further QC steps.
# This block does not export an image file, overwrite counts, select a final
# statistical model, run differential expression or refresh the checkpoint.
