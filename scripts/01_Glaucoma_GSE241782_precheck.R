# Portable repository copy. Start R in the repository root.
# This file records the interactive workflow; later result/closure notes
# supersede earlier pending notes. See the execution guides under docs/.
# Original working files are preserved. Repository sync did not run analysis.
# Glaucoma GSE241782: pre-quality checks
# Started: 26 October 2026
# Run one small block at a time and review its output before proceeding.
# Completed glaucoma precheck; figures restyled on 29 October 2026.

# STEP 1: Confirm that R can find the dataset folder and both input files.
# The dataset folder has been renamed to "Glaucoma_GSE241782".
# The downloaded files are already uncompressed: .csv and .soft.

# <- assigns the folder path to an R object named glaucoma_dir.
glaucoma_dir <- Sys.getenv("GLAUCOMA_DATA_DIR", unset = file.path("data", "raw", "Glaucoma_GSE241782"))

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
# Added: 27 October 2026. Run in the author's restored R session.
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
# Numerical totals must be reviewed from the author's actual output.
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
# Review the author's actual values before drawing sample-quality conclusions.
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

# STEP 9D: Save the same initial boxplot as a PNG image.
# Added: 27 October 2026 at the author's explicit request to save the boxplot.
# Run after Step 9C has created glaucoma_log_cpm in the active R session.
# The existing directory is prequalitychecks/glaucoma (not prechecks/glaucoma).
# Start R in the repository root when using relative output paths.
glaucoma_precheck_dir <- file.path("results", "glaucoma")
glaucoma_boxplot_file <- file.path(
  glaucoma_precheck_dir,
  "Glaucoma_GSE241782_logCPM_boxplot.png"
)

# png() opens a graphics device that writes the following plot to this file.
# Width and height are pixels; res specifies 300 pixels per inch.
# These settings give an 8-by-6-inch image at 300 dpi.
png(
  filename = glaucoma_boxplot_file,
  width = 2400,
  height = 1800,
  res = 300
)

# Draw the same boxplot on the PNG device instead of the graphics window.
boxplot(
  glaucoma_log_cpm,
  las = 2,
  col = "lightblue",
  outline = FALSE,
  main = "Retinal libraries: expression distributions",
  xlab = "Library",
  ylab = "Log2 counts per million (CPM)"
)

# Close the PNG device to finish writing the image file.
dev.off()

# Expected: TRUE, confirming the PNG exists after executing this block.
file.exists(glaucoma_boxplot_file)

# Rerunning this block replaces the PNG with the same filename.
# This image export does not save or refresh the R workspace checkpoint.
# Await the author's saved plot for visual interpretation before further QC.

# STEP 9C OBSERVATION (user's boxplot screenshot, 27 October 2026):
# The eight initial log-CPM boxes have broadly comparable shapes and upper
# quartiles; R-08 is not grossly separated in this marginal-distribution view.
# R-04, R-08 and R-09 have median lines coinciding with the lower floor.
# This is consistent with their zero fractions, calculated from the author's
# detected-gene summaries: about 50.4%, 53.8% and 52.6%, respectively.
# Keeping all 55,487 gene rows means zeros affect these medians strongly.
# For the reported mean library total of 35,261,343 and prior.count = 2,
# zero counts map to a common log-CPM floor near -4.14. This is the plotting
# transformation, not a negative raw count or evidence of detected expression.
# Similar marginal distributions do not establish matching gene profiles.
# Retain all eight libraries for the remaining planned checks; do not treat
# this initial boxplot as final sample-quality approval.
#
# The author requested a conceptual explanation of CPM, log2 and the prior.
# CPM expresses a gene's count per million total gene counts in that library.
# log = TRUE returns a new matrix of log2-CPM values, with unchanged row/column
# identities. Negative logarithms are possible for positive values below 1.
# prior.count = 2 is the installed edgeR cpm() default (verified from local
# edgeR 4.8.2 function formals). It is an average prior, scaled by library size,
# rather than adding 2 CPM or 2 fixed counts to every library. It allows
# finite log-values for zeros and reduces instability at very low counts.
# Original counts and the raw-count-based detected-gene definition are retained.
# No new executable analysis was added during this explanation; PNG export
# success is not established from the screenshot alone.

# STEP 9E: Overlay density curves for the eight libraries' log-CPM values.
# Added: 27 October 2026 at the author's explicit request for density plots.
# Use the same glaucoma_log_cpm matrix as the initial boxplot: all 55,487
# rows, including zeros. No new scaling or expression filtering occurs here.
# Reference for this QC component:
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# density() estimates a smooth distribution curve for one numeric vector.
# apply(..., 2, ...) operates on columns, one library at a time.
# simplify = FALSE retains one density result per column in a list.
glaucoma_density_curves <- apply(
  glaucoma_log_cpm,
  2,
  density,
  simplify = FALSE
)

# Each density result contains horizontal positions ($x) and curve heights
# ($y). sapply() extracts these into two matrices with one column per library.
glaucoma_density_x <- sapply(glaucoma_density_curves, function(d) d$x)
glaucoma_density_y <- sapply(glaucoma_density_curves, function(d) d$y)

# Assign one visible colour per library.
glaucoma_density_colours <- hcl.colors(
  ncol(glaucoma_log_cpm),
  palette = "Dark 3"
)

# Bundle the drawing instructions so the same plot can be displayed in R
# and written to a PNG without copying the plotting code twice.
glaucoma_draw_density <- function() {
  # matplot() draws paired columns of x/y as separate curves on shared axes.
  # Its default axis limits include the ranges of all eight curves.
  # type = "l" draws lines; lty = 1 makes them solid; lwd = 2 sets thickness.
  matplot(
    glaucoma_density_x,
    glaucoma_density_y,
    type = "l",
    lty = 1,
    lwd = 2,
    col = glaucoma_density_colours,
    main = "Retinal libraries: log-CPM density",
    xlab = "Log2 counts per million (CPM)",
    ylab = "Density"
  )
  legend(
    "topright",
    legend = colnames(glaucoma_log_cpm),
    col = glaucoma_density_colours,
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.8
  )
}

# Draw the plot in R's graphics window.
glaucoma_draw_density()

# Each curve is a smoothed distribution, not a count of detected genes.
# Its area over the whole expression scale is approximately one. Peaks show
# where gene-expression values are concentrated; height depends on smoothing.
# The large low-expression peak near -4.14 is expected from the included zeros.
# Smoothing can extend curves slightly beyond the observed expression range.
# Similar curves do not establish matching gene-specific expression profiles.
# Review the author's plot before deciding whether any differences need follow-up.

# STEP 9F: Save the same density plot in the agreed glaucoma output folder.
glaucoma_density_file <- file.path(
  file.path("results", "glaucoma"),
  "Glaucoma_GSE241782_logCPM_density.png"
)
png(filename = glaucoma_density_file, width = 2400, height = 1800, res = 300)
glaucoma_draw_density()
dev.off()

# Expected: TRUE after the author's successful export.
file.exists(glaucoma_density_file)

# Rerunning the export replaces the PNG with the same filename. The original
# count matrix and metadata are retained. No final normalization, sample
# exclusion, gene annotation, differential expression or checkpoint refresh
# is performed in these plotting/export blocks.

# DENSITY EXPORT TROUBLESHOOTING (recorded interactive R output, 27 October 2026):
# The author's first png() call lacked a closing parenthesis; later calls fixed
# that syntax. The export path was typed as prequality/glaucoma, which does
# not exist; the confirmed folder is prequalitychecks/glaucoma as above.
# Typing dev.off without () prints the function instead of closing/writing
# the PNG device. Use dev.off(). Use file.exists(), with the final s.
# For this interactive recovery only, graphics.off() can close unfinished
# plotting devices before rerunning the correct Step 9F block. It also closes
# plot windows, while preserving counts, metadata and drawing functions.
# Do not add graphics.off() to the routine export block as a required step.
# No density PNG was present in the correct folder when inspected. The
# boxplot PNG was present; the user still needs to run the corrected density
# export and provide the density plot for review.

# STEP 9E-9F RESULT (user's density PNG, 27 October 2026):
# The saved Glaucoma_GSE241782_logCPM_density.png is now present in the correct
# glaucoma output folder (about 236 KB). The author supplied it for interpretation.
# All eight curves share a large left peak near the zero-count floor (-4.14)
# and broadly overlapping distributions through the positive-expression range,
# including a broad smaller peak around log2-CPM 5. R-08's green left peak is
# somewhat higher, consistent with its greater fraction of zero rows; peak
# heights also depend on the density smoothing and are not percentages.
# R-11 has a modestly larger low-positive shoulder, consistent with its higher
# detected-gene total. These are visual observations, not gene-level tests.
# No library is grossly shifted relative to the rest in this marginal view.
# This is reassuring for the initial distribution screen. Retain all eight
# libraries provisionally for the remaining planned prechecks, including R-08.
# These plots do not prove full sample quality, equal gene-specific profiles,
# absence of disease responses, or disease resemblance to the flight reference.
# The small number of pooled libraries per stage remains a design limitation.
# Initial library-total, detected-gene and distribution summaries are reviewed.
# Next: gene identifiers/annotation, then later normalized sample-profile checks.
# No new executable analysis or checkpoint refresh was added during this review.

# STEP 10A: Inspect gene labels and count mouse-Ensembl-prefix rows.
# Added: 27 October 2026. Run manually in the active R session.
# R-08 remains included based on the completed initial checks; final
# gene-specific sample-profile QC is still part of the planned precheck.
#
# Why: later cross-study comparisons must match the same genes, not just
# samples. Observed labels include ENSMUSG00000025902|Sox17: an identifier
# joined to a gene symbol. Inventory formats before extracting/mapping IDs.
# Ensembl documents ENSMUSG as the mouse gene-ID prefix:
# https://asia.ensembl.org/info/genome/stable_ids/prefixes.html
# Source SOFT processing records specify GRCm38 and GENCODE M25 and report
# reference plus de novo transcripts. Official release record:
# https://www.gencodegenes.org/mouse/release_M25.html
# Do not assume each assembled row has a standard reference-gene identifier.

# rownames() retrieves the original labels in count-matrix gene order.
glaucoma_gene_labels <- rownames(glaucoma_counts_selected)
head(glaucoma_gene_labels, n = 10)

# grepl() returns one TRUE/FALSE value per label. ^ means starts with.
# This inventories a prefix; it is not full ID validation or annotation mapping.
glaucoma_ensembl_prefix <- grepl(
  pattern = "^ENSMUSG",
  x = glaucoma_gene_labels
)

# TRUE: rows starting with the standard mouse Ensembl gene-ID prefix.
# FALSE: rows with another naming format, to inspect before deciding handling.
# Counts here refer to gene rows, not sequenced libraries.
table(Ensembl_prefix = glaucoma_ensembl_prefix)

# ! reverses TRUE/FALSE. Show up to ten labels without the Ensembl prefix.
# character(0) would mean there are no such labels; it would not be an error.
head(glaucoma_gene_labels[!glaucoma_ensembl_prefix], n = 10)

# Await the author's label examples and prefix totals before extracting IDs and
# symbols, checking duplicates in parsed identifiers, or choosing annotation
# handling. Original gene labels and counts remain intact; no genes are removed,
# no reference annotation is downloaded and no cross-study matching is run here.

# STEP 10A RESULT (user's corrected R screenshot, 27 October 2026):
# The correctly spelled ^ENSMUSG pattern returned TRUE for all 55,487 labels.
# head(glaucoma_gene_labels[!glaucoma_ensembl_prefix], n = 10) returned
# character(0): an empty character vector because there are no prefix failures.
# The earlier all-FALSE output used the mistyped ^ESNMUSG pattern and is not
# a valid identifier inventory. The canonical variable includes ensembl in full.
# ! reverses the match flag; the display command searches for exceptions.
# It neither strips the ENSMUSG prefix nor removes rows from the count matrix.
# This confirms prefix format only, not full identifier validity or successful
# reference annotation. Next: inspect/extract identifier and symbol components,
# then check parsed IDs before choosing any gene-matching or annotation handling.
# No new executable analysis was added during this conceptual explanation.

# STEP 10B: Separate the original labels into Ensembl IDs and gene symbols.
# Added: 27 October 2026 at the author's request to proceed with extraction.
# This constructs a separate gene-annotation table using depositor-provided
# labels. It does not yet validate IDs against an external reference, replace
# count-matrix row names, merge genes, filter rows or alter any counts.

# strsplit() divides each label into components. fixed = TRUE treats | as a
# literal separator rather than a regular-expression operator.
# The result is a list with one item per gene row, in the same order.
glaucoma_gene_parts <- strsplit(
  x = glaucoma_gene_labels,
  split = "|",
  fixed = TRUE
)

# lengths() counts components within each list item. The observed ID|symbol
# format should give two components; inspect the table for any other formats.
table(Fields_per_label = lengths(glaucoma_gene_parts))

# For each gene's components, [1] selects the ID and [2] the symbol.
# sapply() applies that selection to each list item and returns a vector.
# If a label lacks a second component, [2] gives NA rather than inventing a
# symbol. The field-count table above will reveal labels with other lengths.
# Keep original_label so any missing or extra components can be reviewed.
glaucoma_gene_annotation <- data.frame(
  original_label = glaucoma_gene_labels,
  ensembl_id = sapply(glaucoma_gene_parts, function(x) x[1]),
  gene_symbol = sapply(glaucoma_gene_parts, function(x) x[2]),
  stringsAsFactors = FALSE
)

# Preview the original label and its extracted components for five gene rows.
head(glaucoma_gene_annotation, n = 5)

# Expected: 55487 rows and 3 columns, preserving count-matrix gene-row order.
dim(glaucoma_gene_annotation)

# Review field counts, the preview and dimensions before full ID-format,
# missing-field and duplicate-ID checks. Extra delimiters or missing symbols
# need interpretation; do not discard or combine gene rows automatically.

# STEP 10B OBSERVATION (recorded interactive R output, 27 October 2026):
# glaucoma_gene_annotation has 55,487 rows and three columns. The first five
# rows show the expected separation into original_label, ensembl_id and symbol.
# The screenshot does not include the corrected per-label field-count table.
# Extraction is complete; full ID structure, missingness, uniqueness and
# row correspondence have not yet been validated from the author's outputs.

# STEP 10C: Validate the extracted gene-annotation structure before sample QC.
# Added: 27 October 2026. Run manually in the active R session.
# These checks describe structural consistency, not external database matching.
# Full functional annotation is not required to inspect sample relationships;
# reference/version harmonization must be addressed before cross-study matching.

# Confirm every label's number of split components. Expected: 2 for 55487 rows.
# This uses lengths(), not length(): count components within each list item.
table(Fields_per_label = lengths(glaucoma_gene_parts))

# Count missing values in each annotation column. Expected: three zeros.
colSums(is.na(glaucoma_gene_annotation))

# Count empty/whitespace-only fields per column. Expected: three zeros.
# trimws() removes surrounding whitespace for this check only.
# na.rm = TRUE skips NAs here; their counts are reported separately above.
sapply(glaucoma_gene_annotation, function(x) {
  sum(trimws(x) == "", na.rm = TRUE)
})

# Check repeated extracted Ensembl IDs, rather than repeated original labels.
# Expected: 0. A positive result is the first duplicate's position, not a count.
anyDuplicated(glaucoma_gene_annotation$ensembl_id)

# Check the full ID format: mouse gene prefix plus eleven digits, allowing
# an optional .version suffix. ^ and $ anchor the beginning and end.
# This formal pattern does not establish that an ID exists in a given release.
# https://mart.ensembl.org/info/genome/stable_ids/index.html
glaucoma_gene_id_format_ok <- grepl(
  pattern = "^ENSMUSG[0-9]{11}(\\.[0-9]+)?$",
  x = glaucoma_gene_annotation$ensembl_id
)
table(Complete_ID_format = glaucoma_gene_id_format_ok)
# Expected: TRUE for all 55487 rows, if all extracted IDs use this format.

# Confirm every annotation row still corresponds to the same count-matrix row.
# Expected: TRUE.
identical(
  glaucoma_gene_annotation$original_label,
  rownames(glaucoma_counts_selected)
)

# Await and interpret these results before moving on. A missing symbol or
# other anomaly needs inspection, not automatic gene/sample removal. Symbols
# are not required to be unique; the Ensembl ID is the intended matching key.
# Next: prepare appropriately filtered/scaled expression values for the planned
# sample-profile comparisons. Current all-row log-CPM is an initial plotting
# view, not a statement that final filtering/normalization is complete.
# No external annotation lookup, gene filtering, final normalization, sample
# correlation, PCA/MDS or checkpoint export is executed in this block.

# STEP 10C RESULT (recorded interactive R output, 27 October 2026):
# NA counts: original_label 0, ensembl_id 0, gene_symbol 0.
# Empty/whitespace-only counts: 0 in every column.
# anyDuplicated(extracted ensembl_id): 0.
# Complete ID-format check: TRUE for all 55,487 rows.
# Original-label/count-row alignment: TRUE.
# These shown structural checks pass. This is not external database membership
# validation or completed OSD-758 gene-ID harmonization.

# STEP 11A: Prepare an expressed-gene QC copy for sample-profile comparisons.
# Added: 27 October 2026 at the author's request to inspect sample relationships.
# Start with a standard edgeR preparation before normalization and MDS.
# Rationale and method reference:
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# Very low-count genes can contribute little reliable profile information.
# filterByExpr() considers library sizes and the smallest biological group size.
# Its use of group sample size can retain genes supported only at one stage;
# it does not require expression in every group or select DE-significant genes.
# The adaptation here is to use four treatment/time groups of two libraries each,
# rather than treating all six glaucoma libraries as one collection-time group.
# Sex is retained in sample metadata for later visual interpretation.

# paste() joins treatment and collection time; factor() marks categorical groups.
# Take names directly from the metadata, preserving the accented naive label.
glaucoma_qc_group <- factor(
  paste(
    glaucoma_metadata_selected$treatment,
    glaucoma_metadata_selected$time,
    sep = "_"
  )
)

# Expected: four groups, each containing two pooled libraries.
table(glaucoma_qc_group)

# DGEList stores numeric counts, sample metadata and gene annotation in separate
# components of one object, with the previously validated ordering retained.
# This does not join descriptive text into the numeric count matrix.
glaucoma_dge <- edgeR::DGEList(
  counts = glaucoma_counts_selected,
  samples = glaucoma_metadata_selected,
  group = glaucoma_qc_group,
  genes = glaucoma_gene_annotation
)

# Use edgeR 4.8.2's standard filterByExpr defaults for this initial QC copy:
# min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7.
# With two libraries per group, the minimum group size is two. The expression
# threshold is library-size-aware; it is not a fixed 10 counts in every sample.
# More precisely: require CPM >= 10 / median(lib.size) * 1e6 in at least two
# of the eight libraries, plus at least 15 total counts across all libraries.
# The two qualifying libraries need not belong to the same named group.
# Group labels set the minimum support size; no treatment difference is tested.
glaucoma_keep_genes <- edgeR::filterByExpr(
  glaucoma_dge,
  group = glaucoma_qc_group
)

# TRUE: retain the gene for profile QC; FALSE: insufficient expression support.
# FALSE does not mean a bad gene identifier or a failed library.
# Counts across these two categories should sum to 55,487.
table(Retain_for_QC = glaucoma_keep_genes)

# Select gene rows while retaining every sample column. keep.lib.sizes = FALSE
# recalculates library totals from the retained genes for the QC object.
# DGEList subsetting also keeps the corresponding gene-annotation rows aligned.
glaucoma_dge_qc <- glaucoma_dge[
  glaucoma_keep_genes,
  ,
  keep.lib.sizes = FALSE
]

# Expected: number of retained genes by 8 libraries. Review the actual gene
# count from the author's output; no retained-gene count has been assumed.
dim(glaucoma_dge_qc)

# Await the author's group totals, keep table and dimensions before proceeding
# to composition normalization and the first MDS view (one point per pooled
# library). The original selected counts, full annotation table and all eight
# libraries remain available. No final DE model is chosen by this QC step.
# No normalization, MDS/PCA, sample correlations or export is executed here.

# PAUSE CHECKPOINT — 27 October 2026
# Run ONLY this save block in the current R console when taking a break.
# Do not source the whole script or reload the older workspace first.
# Latest confirmed results: selected retinal counts/metadata are aligned;
# library totals, detected genes, boxplot/density and shown structural-ID checks
# have been reviewed. All eight selected libraries, including R-08, are retained.
# Step 11A execution/results have not been confirmed. MDS/PCA interpretation
# still needs a simpler explanation before continuing sample-relationship QC.

glaucoma_pause_file <- file.path("results", "glaucoma", "Glaucoma_GSE241782_precheck_checkpoint_2026-10-27.RData")

# Save all objects currently in the global workspace, including counts,
# metadata, annotation, calculated summaries and any QC objects already made.
# The author decided to keep a backup too. Save the current workspace to this
# dated file while retaining the original named checkpoint as currently saved.
# Rerunning this block replaces the dated file. Plot PNGs are separate files.
save.image(file = glaucoma_pause_file)
file.exists(glaucoma_pause_file)

# AFTER REOPENING R: run the following load command manually (not now).
# load(file.path("results", "glaucoma", "Glaucoma_GSE241782_precheck_checkpoint_2026-10-27.RData"))
# Resume with an explanation of the purpose of sample-relationship QC, then
# check which preparation objects were saved before repeating any calculations.
# The assistant saved these instructions, not the author's live R workspace;
# updating the checkpoint requires the user to run the save block above.

# STEP 11 RESUME — 28 October 2026
# The author reports reloading the saved workspace and wants to start MDS.
# First check which preparation objects were restored, without repeating work.
# exists() checks for an object in R memory; it does not check a disk filename.
# Expected: the first three TRUE; filtered_QC may be TRUE or FALSE because
# execution of Step 11A was not confirmed before the pause.
c(
  selected_counts = exists("glaucoma_counts_selected"),
  selected_metadata = exists("glaucoma_metadata_selected"),
  gene_annotation = exists("glaucoma_gene_annotation"),
  filtered_QC = exists("glaucoma_dge_qc")
)

# Review the author's output before continuing: if filtered_QC is FALSE, prepare
# the expressed-gene copy using Step 11A. If TRUE, inspect its dimensions and
# alignment, then apply/review appropriate composition normalization for MDS.
# A TRUE here reports object presence, not validation or completed QC.
# Do not run the entire script or save/reload a workspace during this check.

# STEP 11A RESULT AND EXPLANATION — 28 October 2026
# User's screenshot confirms four treatment/time categories with two libraries
# each: Glaucoma_2W, Glaucoma_3D, Glaucoma_6W and naive_0D (original accent kept).
# filterByExpr result: FALSE 32,186; TRUE 23,301 (55,487 original gene rows).
# The initial subsetting call misspelled keep.lib.sizes as kep.lib.sizes and
# failed; the corrected call succeeded. dim(glaucoma_dge_qc): 23,301 by 8.
# These are gene rows excluded/retained in the QC copy; no library was removed.
# Original glaucoma_counts_selected, glaucoma_dge and full annotation remain.
# paste() joins treatment/time row by row; factor() marks categories. No counts
# are averaged, combined, reordered or tested by making glaucoma_qc_group.
# DGEList bundles separate counts, sample metadata/group and gene annotation.
# filterByExpr makes a logical row mask; [mask, , keep.lib.sizes = FALSE]
# applies it to counts and annotation while retaining all sample columns and
# recalculating totals from retained gene counts. It does not normalize counts.
# Filtering removes all-zero and insufficiently supported low-expression genes;
# zeros in some libraries do not by themselves exclude a gene.
# Group sizes set a two-library support requirement anywhere among the eight,
# not specifically both members of a condition or expression in all groups.
# This is preparation for sample-profile QC, not a glaucoma-specific gene list.
# Explain this block before adding normalization/MDS instructions. No further
# expression analysis or checkpoint save was performed by the assistant.

# STEP 11B: Estimate composition-normalization factors before MDS.
# Added: 28 October 2026. Run this small block manually, then review the table.
# User's latest screenshot confirms selected counts, metadata, annotation and
# filtered QC objects all exist. The preceding result was 23,301 genes by 8.
# Basic CPM adjusts for total count differences. If a subset of genes takes a
# larger share of a library's counts, relative counts of other genes can shift.
# TMM (trimmed mean of M-values) estimates a composition scaling factor per
# library using many genes while trimming extreme ratios/abundances.
# It assumes most shared genes are not differentially expressed; it does not
# force matching sample profiles or determine treatment effects.
# Method: https://doi.org/10.1186/gb-2010-11-3-r25
# Workflow: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# normLibSizes() is the current name for calcNormFactors(), and is available in
# the installed edgeR 4.8.2. This call stores factors in samples$norm.factors;
# the retained raw-count values and the eight library identities are preserved.
glaucoma_dge_qc <- edgeR::normLibSizes(
  glaucoma_dge_qc,
  method = "TMM"
)

# Display the library identifier, retained-gene total and estimated factor.
# lib.size is the total after filtering, so it can differ from Step 9A totals.
# A factor near 1 implies a small additional composition adjustment; factors
# above/below 1 are scaling information, not library pass/fail scores.
glaucoma_dge_qc$samples[
  , c("sample_id", "lib.size", "norm.factors")
]

# Await the author's eight factor values before preparing the first MDS plot.
# There is no automatic exclusion based on a normalization factor. The next
# MDS step will use this DGEList so its effective library sizes include TMM.
# No MDS/PCA, gene tests, plot export or checkpoint refresh is executed here.

# STEP 11B RESULT — 28 October 2026 (recorded interactive R output)
# sample_id  retained lib.size  norm.factors
# R-03       40917804           1.0329817
# R-04       32492167           0.9984731
# R-07       46912130           0.9757287
# R-08       19917066           1.0112024
# R-09       26581736           1.0054352
# R-10       34184040           1.0387825
# R-11       43724024           0.9665911
# R-12       37130025           0.9733799
# All factors are positive and near 1 (range about 0.967 to 1.039).
# Their neutral/reference value is 1, not 0. A factor below 1 is not negative
# and is not a failed-quality score; factors scale effective library sizes:
# effective library size = retained lib.size * norm.factors.
# Factor < 1 decreases that denominator relative to basic library-size CPM;
# factor > 1 increases it. Retained raw-count entries remain unchanged.
# Additional composition adjustments are modest. R-08's factor is about 1.011;
# this does not independently prove its full quality or explain its lower total.
# TMM estimates these factors using weighted, trimmed gene log-ratios; trimming
# for estimation does not remove further rows from this DGEList count matrix.
# User requested clarification of TMM and scaling vs quality. Explain before
# advancing; no MDS or checkpoint refresh is confirmed yet.

# STEP 11C: Draw the first MDS view of the eight retinal libraries.
# Added: 28 October 2026, after reviewing the author's TMM factors.
# Run manually in the current R session. Input: 23,301 retained gene rows,
# eight libraries, with the reviewed TMM factors stored in the DGEList.
# Method reference (section 5.4):
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

# Give each treatment/time category a colour. setNames() associates colours
# with the four existing group names. Indexing by each library's group name
# then gives eight colours in the same order as the DGEList's sample columns.
glaucoma_mds_palette <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  levels(glaucoma_dge_qc$samples$group)
)
glaucoma_mds_colours <- glaucoma_mds_palette[
  as.character(glaucoma_dge_qc$samples$group)
]

# The exported generic is limma::plotMDS(), which dispatches to edgeR's DGEList
# method. edgeR::plotMDS is not exported and must not be used.
# The default logFC method
# computes log2-CPM using the stored TMM-adjusted effective library sizes.
# It uses the default library-scaled prior.count = 2 for finite/stable values.
# top = 500 and pairwise are the established defaults: use the 500 most
# divergent retained genes for each pair of libraries when measuring distance.
# These are not genes selected by a statistical test; the set may differ by
# library pair. No group label or sex is used to position libraries.
# labels joins each sample ID and sex for reading the plot; colours identify
# treatment/time only. Distances in two dimensions approximate the profile
# differences. This call displays the plot and stores the returned MDS object.
glaucoma_mds <- limma::plotMDS(
  glaucoma_dge_qc,
  top = 500,
  gene.selection = "pairwise",
  labels = paste(
    glaucoma_dge_qc$samples$sample_id,
    glaucoma_dge_qc$samples$sex
  ),
  col = glaucoma_mds_colours,
  cex = 0.8,
  main = "Retinal libraries: MDS after TMM"
)

# Add a text-colour key for the four treatment/time categories.
legend(
  "topright",
  legend = names(glaucoma_mds_palette),
  text.col = glaucoma_mds_palette,
  bty = "n",
  cex = 0.8
)

# Expected display: eight library labels, with IDs/sex and four group colours.
# Review the author's plot before interpreting separation or possible unusual
# profiles. Sex/time effects are plausible biology; an isolated library is not
# automatically excluded, and overlap does not prove absence of DE genes.
# User execution/output review is pending. No MDS calculation was performed
# by the assistant; no figure export, gene testing or checkpoint save here.

# STEP 11D: Save the MDS figure as a PNG.
# Added: 28 October 2026. User preference: always include explicit figure-save
# code with future plotting instructions. Export is run manually by the user.
# Recovery revision: run this entire block after Step 11B. It calculates the
# MDS and saves it together; a pre-existing glaucoma_mds object is not needed.
# It recreates the display colours, without re-filtering or repeating TMM.
glaucoma_mds_palette <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  levels(glaucoma_dge_qc$samples$group)
)
glaucoma_mds_colours <- glaucoma_mds_palette[
  as.character(glaucoma_dge_qc$samples$group)
]

glaucoma_mds_file <- file.path(
  file.path("results", "glaucoma"),
  "Glaucoma_GSE241782_TMM_MDS.png"
)

# Open the PNG output device: 2400 x 1800 pixels at 300 dots per inch.
png(
  filename = glaucoma_mds_file,
  width = 2400,
  height = 1800,
  res = 300
)

# Calculate MDS directly from the filtered/TMM DGEList and draw into the PNG.
# This creates glaucoma_mds instead of requiring that it already exists.
glaucoma_mds <- limma::plotMDS(
  glaucoma_dge_qc,
  top = 500,
  gene.selection = "pairwise",
  labels = paste(
    glaucoma_dge_qc$samples$sample_id,
    glaucoma_dge_qc$samples$sex
  ),
  col = glaucoma_mds_colours,
  cex = 0.8,
  main = "Retinal libraries: MDS after TMM"
)
legend(
  "topright",
  legend = names(glaucoma_mds_palette),
  text.col = glaucoma_mds_palette,
  bty = "n",
  cex = 0.8
)

# Close the PNG device to finish writing the file. Parentheses are required.
dev.off()

# Expected: TRUE after a successful export. Rerunning replaces this PNG.
file.exists(glaucoma_mds_file)

# Only plotting/export instructions were added by the assistant. Actual MDS
# calculation, the PNG file and its interpretation await user execution/review.
# This PNG export does not update the saved R workspace checkpoint.

# USER EXPORT ERROR — 28 October 2026:
# The first export attempt used the original redraw-only block before
# glaucoma_mds had been created: object-not-found, then no plot for legend().
# dev.off() returned null device 1, so that failed PNG device is already closed.
# file.exists() returned FALSE; filesystem inspection confirmed no MDS PNG.
# Use the revised complete Step 11D above. MDS/PNG success is not yet confirmed.

# NAMESPACE CORRECTION — 28 October 2026:
# The assistant's earlier edgeR::plotMDS prefix was wrong: plotMDS is exported
# by limma; edgeR registers the DGEList method. Both calls above are corrected
# to limma::plotMDS. This uses the same intended edgeR DGEList method, not a
# change in normalization or distance settings. Verified namespace exports,
# registered method and calculation/legend/PNG export with synthetic dummy
# counts only; no user data were loaded or analysed by the assistant.
# User's failed device was closed with dev.off(); rerun the corrected Step 11D
# to reopen PNG, calculate/draw MDS, add legend, close, and check file existence.

# STEP 11C/11D RESULT — 28 October 2026:
# User supplied Glaucoma_GSE241782_TMM_MDS.png; file presence confirmed (~207 KB).
# MDS calculation/export now succeeded. The displayed axes are dimensions 1
# (19%) and 2 (18%), a partial two-dimensional view (~37% of MDS variation).
# Visual observations only; no coordinates or expression data recomputed here:
# - All four male libraries appear to the right of all four female libraries
#   along dimension 1. Each male is right of the female in its condition/time.
#   This supports sex-associated structure, not proof of biological causation
#   or a claim that sex alone explains the observed distances.
# - R-04 (2W female) and R-07 (3D male) appear close in these two dimensions;
#   condition/time groupmates do not form tight nearest-neighbour pairs.
#   This is not evidence that collection time has no gene-level effects.
# - R-08 (3D female) is upper left and R-09 (6W male) upper right. Both have high
#   dimension-2 coordinates; no single sex explanation accounts for all layout.
# - One pooled library per sex/condition-time cell limits separation of sources
#   of variability. These male/female groupmates are not same-animal repeated
#   measurements. Time labels describe tissue collection after induction.
# Percentages label MDS axes, not percentages attributed to sex or treatment.
# Preserve sex/time for model specification. Original study adjusted for sex:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10531004/ (Methods 4.9).
# Keep all eight libraries pending the remaining planned profile checks; R-08's
# distinct profile warrants follow-up, not automatic removal or failed status.
# No new analysis/model change or checkpoint refresh during this explanation.

# MDS PERCENTAGES CLARIFICATION — 28 October 2026:
# User asks whether ~37% was manually chosen, whether it is standard, and what
# this figure alone establishes. No new executable analysis was added.
# We used the default first two dimensions; their displayed ~19% and ~18%
# were calculated from the MDS distances/eigenvalues, not set by us.
# Installed limma implementation uses nonnegative eigenvalues divided by their
# sum for var.explained; plotted percentages are rounded. Thus ~37% describes
# the spread represented by these axes in this particular distance-based MDS,
# not a percentage of genes, raw counts, sample quality or causal sex effects.
# Showing dimensions 1 and 2 is a standard starting view. There is no universal
# requirement that the displayed fraction equal 37%, nor a pass/fail cutoff.
# The remaining directions are not automatically noise or irrelevant biology.
# Defensible figure-only observations: sex-associated horizontal structure;
# no tight condition/time pairs in this projection; R-08/R-09 distinct in dim2.
# Causes, statistical significance, absence of time effects and sample removal
# do not follow from the plot alone. Keep the agreed design and remaining QC.

# STEP 11E: PCA using one common set of 500 variable genes; draw and save.
# Added: 28 October 2026. The author explicitly requests PCA as a complementary
# view of internal expression patterns. Retain the existing MDS result.
# Input is the same filtered/TMM DGEList: 23,301 eligible genes, eight libraries.
# In limma's documented plotMDS() implementation, gene.selection = "common"
# constructs PCA using the top genes with largest standard deviations across
# all libraries. This differs from the pair-specific selection in our MDS.
# The DGEList method supplies TMM-adjusted log2-CPM with prior.count = 2.
# Genes are selected without treatment/sex labels, centred for PCA, and are
# not scaled individually to unit variance. 500 matches the standard top-gene
# setting already used for the MDS; it is not tuned to force group separation.
# Reference: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
# (plotMDS Details: gene.selection = "common" constructs a PCA.)

# Use the same four group colours for comparability with the MDS view.
glaucoma_pca_palette <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
  levels(glaucoma_dge_qc$samples$group)
)
glaucoma_pca_colours <- glaucoma_pca_palette[
  as.character(glaucoma_dge_qc$samples$group)
]

# Save this complementary figure under its own explicit filename.
glaucoma_pca_file <- file.path(
  file.path("results", "glaucoma"),
  "Glaucoma_GSE241782_TMM_PCA_top500.png"
)

png(
  filename = glaucoma_pca_file,
  width = 2400,
  height = 1800,
  res = 300
)

# Despite the function name, "common" selects its PCA mode. Each displayed
# label remains one pooled library; the same selected 500 genes inform all.
glaucoma_pca <- limma::plotMDS(
  glaucoma_dge_qc,
  top = 500,
  gene.selection = "common",
  labels = paste(
    glaucoma_dge_qc$samples$sample_id,
    glaucoma_dge_qc$samples$sex
  ),
  col = glaucoma_pca_colours,
  cex = 0.8,
  main = "Retinal libraries: PCA after TMM"
)

legend(
  "topright",
  legend = names(glaucoma_pca_palette),
  text.col = glaucoma_pca_palette,
  bty = "n",
  cex = 0.8
)

dev.off()
file.exists(glaucoma_pca_file)

# Expected: TRUE after successful export; open the PNG for interpretation.
# Rerunning replaces this PCA PNG. This does not update the RData checkpoint.
# Implementation note: glaucoma_pca has limma's MDS container class even in
# PCA mode. $x and $y are PC1/PC2 coordinate vectors, not a prcomp score matrix;
# $var.explained contains the fractions per component. Gene loadings are not
# returned by this plotting interface. Do not mistake it for a prcomp object.
# This PCA describes within-glaucoma relationships. Future OSD-758 QC may use
# its own PCA; separate PCA coordinates are not shared disease/flight scores.
# Cross-study resemblance remains based on within-study response comparisons.
# No user expression data were analysed by the assistant. This PCA/export code
# was verified with dummy data only; user execution and results are pending.

# STEP 11E RESULT — 28 October 2026:
# The author supplied a screenshot of the PCA; its PNG file exists (~211 KB).
# PC1 shows 23% and PC2 20% (rounded): about 43% of variation in the common
# 500 selected genes' log-expression values is represented in this view.
# Visual observations, not a new calculation from user data:
# - All four male libraries lie to the right of all four female libraries on
#   PC1, including within every treatment/time category, as in the MDS view.
# - Some neighbourhoods recur: R-10/R-11 are near each other, and R-04/R-07
#   remain relatively near. Do not assert exact full-profile nearest-neighbour
#   rankings from these text-label projections alone.
# - Same-condition/time male/female libraries do not form tight pairs. This
#   does not establish absence of time-dependent or glaucoma-specific changes.
# - R-08 remains distinctive on the second axis. Its bottom position vs the
#   top in MDS is not a biological reversal; component signs are arbitrary.
# Similarity of these views shows persistence under the gene-selection change;
# they use the same eight libraries and are not independent validation data.
# Preserve sex and time in later design considerations; the original paper
# adjusted for sex. Retain all eight pending planned sample-correlation QC.
# Neither PCA nor MDS identifies the cause of R-08's profile or diagnoses a
# failed library; no sample is removed from these plots alone.
# No real-data recomputation, new plots/model fitting or checkpoint save was
# performed by the assistant during this interpretation.

# PRECHECK CLOSED BY USER DECISION — 28 October 2026
# Retain all eight selected retinal libraries, including R-08. Its distinct
# profile remains documented; it has not been established as a failed sample.
# Omit the optional sample-correlation check and move to the next dataset.
# Do not reopen QC or add required follow-ups without a specific new reason.
# Precheck summary saved separately as:
# docs/glaucoma-precheck.md
# No differential-expression or disease/OSD-758 comparison was performed here.

# Run this in the active R session BEFORE beginning the next dataset, to save
# the current objects including filtered/TMM data, annotation, MDS and PCA.
# Keep both older checkpoints as backups; this dated file is the new target.
glaucoma_final_checkpoint <- file.path(
  file.path("results", "glaucoma"),
  "Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData"
)
save.image(file = glaucoma_final_checkpoint)
file.exists(glaucoma_final_checkpoint)

# Restore this newest checkpoint in a future session only AFTER the user has
# confirmed saving it successfully. The assistant has not saved live R objects.
# load(file.path("results", "glaucoma", "Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData"))

# HANDOFF STATUS — 28 October 2026:
# The newest checkpoint now exists (305,345,347 bytes); filesystem presence
# verified while preparing the handoff, without loading or inspecting objects.
# Full record for the new RP thread:
# docs/glaucoma-precheck.md (completed precheck record)
# Glaucoma remains closed with R-08 retained and no required correlation check.

# FIGURE STYLE REFRESH REQUESTED — 29 October 2026:
# User requests glaucoma MDS/PCA in the approved RP PCA style: actual sample
# points, nonoverlapping labels with connector lines, and an outside legend.
# Glaucoma's four treatment/time groups and existing colours are preserved.
# Run the standalone restyling code in the active RP session:
# source(file.path("scripts", "04_Restyle_Retinal_MDS_PCA.R"))
# It reuses the stored projections; if glaucoma objects are absent, it loads
# the dated 28 October checkpoint into an isolated environment. It does not
# load the checkpoint into the current RP workspace or rerun expression QC.
# The script replaces glaucoma MDS/PCA PNGs plus RP MDS. All eight glaucoma
# samples, including R-08, remain as agreed. Exports/visual QA await the user.

# FIGURE STYLE CONFIRMED — 29 October 2026:
# User reports completion; the two restyled PNGs exist with 13:46 timestamps.
# Assistant opened both for layout QA: eight points/sex labels and four
# treatment/time groups remain visible, with outside legends and readable
# labels. MDS axes remain 19%/18%; PCA remains 23%/20%. No reanalysis or
# revised sample decision occurred. The saved RData was not refreshed here.
