# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma"), recursive = TRUE, showWarnings = FALSE)

# Glaucoma GSE241782: within-study analysis
# Started: 3 November 2026.
# Run one small block at a time and review the output before continuing.

# Set the folder where glaucoma within-study outputs will be saved.
glaucoma_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma")
setwd(glaucoma_within_dir)

# Locate the completed glaucoma precheck checkpoint for the next block.
glaucoma_precheck_file <- Sys.getenv("GLAUCOMA_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "glaucoma", "Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData"))

# Confirm the working folder and whether R can find the checkpoint.
getwd()
file.exists(glaucoma_precheck_file)

# User confirmed the correct working folder and TRUE for checkpoint presence.

# Create a container for the restored precheck objects.
glaucoma_precheck <- new.env()

# Restore the checkpoint and record the names of the objects it contains.
glaucoma_loaded_objects <- load(
  glaucoma_precheck_file,
  envir = glaucoma_precheck
)

# Check that the selected metadata and both prepared edgeR objects were restored.
c(
  selected_metadata = "glaucoma_metadata_selected" %in% glaucoma_loaded_objects,
  unfiltered_DGEList = "glaucoma_dge" %in% glaucoma_loaded_objects,
  filtered_TMM_DGEList = "glaucoma_dge_qc" %in% glaucoma_loaded_objects
)

# User confirmed TRUE for all three restored-object checks.

# Make working copies of the filtered/TMM object and selected sample metadata.
glaucoma_within_dge <- glaucoma_precheck$glaucoma_dge_qc
glaucoma_within_metadata <- glaucoma_precheck$glaucoma_metadata_selected

# Display the number of retained gene rows and library columns.
dim(glaucoma_within_dge$counts)

# Review the sample identifiers, treatments, collection times and sexes.
glaucoma_within_metadata[, c("sample_id", "treatment", "time", "sex")]

# Confirm that metadata rows match the count columns in exactly the same order.
identical(
  glaucoma_within_metadata$sample_id,
  colnames(glaucoma_within_dge$counts)
)

# User output confirmed 23,301 gene rows, eight libraries and sample-order TRUE.
# The reviewed metadata has one male and one female library in each group.

# Combine treatment and time into four categories, with naive controls first.
glaucoma_within_metadata$group <- factor(
  paste(
    glaucoma_within_metadata$treatment,
    glaucoma_within_metadata$time,
    sep = "_"
  ),
  levels = c("Naïve_0D", "Glaucoma_3D", "Glaucoma_2W", "Glaucoma_6W")
)

# Treat sex as a category, with Female as its reference level.
glaucoma_within_metadata$sex <- factor(
  glaucoma_within_metadata$sex,
  levels = c("Female", "Male")
)

# Confirm the number of libraries in every group and sex combination.
with(glaucoma_within_metadata, table(group, sex, useNA = "ifany"))

# User table confirmed one female and one male pooled library in all four groups.
# Female remains the sex reference; Naive_0D remains the group reference.

# Build an additive model with a shared sex effect and stage-specific group effects.
# Method: edgeR User's Guide, additive models and blocking.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
glaucoma_design <- model.matrix(
  ~ sex + group,
  data = glaucoma_within_metadata
)

# Label the design rows with their library identifiers and display the matrix.
rownames(glaucoma_design) <- glaucoma_within_metadata$sample_id
glaucoma_design

# Check that every model coefficient can be estimated separately.
qr(glaucoma_design)$rank == ncol(glaucoma_design)

# Calculate the design's residual degrees of freedom for estimating variability.
nrow(glaucoma_design) - qr(glaucoma_design)$rank

# User output confirmed the expected eight-row, five-column design, full rank
# TRUE, treatment-reference coding and three residual degrees of freedom.

# Use the current edgeR quasi-likelihood pipeline (installed edgeR 4.8.2).
# glmQLFit estimates the required dispersions directly in this pipeline.
# Method: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Fit the quasi-likelihood count model using saved TMM factors and the design.
glaucoma_fit <- edgeR::glmQLFit(
  y = glaucoma_within_dge,
  design = glaucoma_design,
  robust = TRUE
)

# Check the number of fitted gene rows and model terms, then display term names.
dim(glaucoma_fit$coefficients)
colnames(glaucoma_fit$coefficients)

# User output confirmed 23,301 fitted gene rows and the five expected model terms.

# Display the gene-variability diagnostic from the fitted model.
edgeR::plotQLDisp(glaucoma_fit, main = "Glaucoma: QL dispersion")

# Build the diagnostic image path inside the within-study output folder.
glaucoma_ql_plot_file <- file.path(
  glaucoma_within_dir,
  "Glaucoma_GSE241782_QL_dispersion.png"
)

# Draw the same diagnostic into a high-resolution PNG file.
png(glaucoma_ql_plot_file, width = 2400, height = 1800, res = 300)
edgeR::plotQLDisp(glaucoma_fit, main = "Glaucoma: QL dispersion")
dev.off()

# Confirm that the diagnostic image was saved.
file.exists(glaucoma_ql_plot_file)

# The user supplied the saved diagnostic for joint interpretation.
# Red estimates form a narrower band around the trend than the initial black
# estimates; the trend generally decreases with abundance, and some genes retain
# high variability. This is consistent with the intended shrinkage behaviour.
# No gene-level significance tests or disease-flight comparisons have been run.
# The user requested saving here before continuing.

# Record the R and package versions used for this saved model.
glaucoma_within_session <- sessionInfo()

# Set a checkpoint filename inside the glaucoma within-study output folder.
glaucoma_within_checkpoint <- file.path(
  glaucoma_within_dir,
  "Glaucoma_GSE241782_withinstudy_checkpoint_2026-11-03.RData"
)

# Save the objects needed to resume without repeating prechecks or model fitting.
save(
  glaucoma_within_dir,
  glaucoma_precheck_file,
  glaucoma_within_dge,
  glaucoma_within_metadata,
  glaucoma_design,
  glaucoma_fit,
  glaucoma_ql_plot_file,
  glaucoma_within_session,
  glaucoma_within_checkpoint,
  file = glaucoma_within_checkpoint
)

# Confirm that the checkpoint exists and report its size in bytes.
file.exists(glaucoma_within_checkpoint)
file.info(glaucoma_within_checkpoint)$size

# The save block was not executed by the assistant. On 4 November 2026, after
# MaxineDrive became available, the dated checkpoint was found at 7,243,016 bytes.
# Current-session restoration remains pending review of the user's R output.
# Use 02_Glaucoma_GSE241782_withinstudy_reload.R to resume from this saved fit.

# Resumed 4 November 2026: the user requested proceeding after the reload block.
# Run this new block in the restored session and review its output together.

# Test the 3-day glaucoma-versus-naive effect while accounting for sex.
glaucoma_3D_test <- edgeR::glmQLFTest(
  glmfit = glaucoma_fit,
  coef = "groupGlaucoma_3D"
)

# Keep results for every tested gene, with BH-adjusted P values, ordered by P value.
glaucoma_3D_results <- edgeR::topTags(
  glaucoma_3D_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Check the total number of tested gene rows and preview the first ten results.
nrow(glaucoma_3D_results)
head(glaucoma_3D_results, n = 10)

# User supplied the first ten 3-day test rows on 4 November 2026; the full-table
# row-count output was not included in the screenshot. The console wrapped
# the annotation/logFC/logCPM columns and F/PValue/FDR columns into two panels.
# Among these displayed rows, FDR < 0.05 is met by Purb (lower), H3f3b, Dpm1,
# 2610001J05Rik, Edn2 and Sprr1a (the latter five higher in 3-day glaucoma).
# Displayed Hipk1, Adarb1, Gfap and Rpl19-ps11 do not meet FDR < 0.05.
# Large estimated log2 changes for Purb (-16.792208) and Sprr1a (11.636210)
# warrant library-level expression review before literal magnitude claims.
# Full-table counts, result exports and the other stage tests remain pending.

# Select the 3-day genes meeting the agreed BH-adjusted P-value criterion.
glaucoma_3D_degs <- subset(glaucoma_3D_results, FDR < 0.05)

# Count all tested genes, significant genes and their expression-change directions.
glaucoma_3D_summary <- c(
  Tested = nrow(glaucoma_3D_results),
  DEGs = nrow(glaucoma_3D_degs),
  Higher = sum(glaucoma_3D_degs$logFC > 0),
  Lower = sum(glaucoma_3D_degs$logFC < 0)
)
glaucoma_3D_summary

# Set the filenames for the complete result table and its significant-gene subset.
glaucoma_3D_results_file <- file.path(
  glaucoma_within_dir,
  "Glaucoma_GSE241782_3D_vs_Naive_all_genes.csv"
)
glaucoma_3D_degs_file <- file.path(
  glaucoma_within_dir,
  "Glaucoma_GSE241782_3D_vs_Naive_DEGs_FDR005.csv"
)

# Save both tables with gene identifiers retained in their existing columns.
write.csv(glaucoma_3D_results, file = glaucoma_3D_results_file, row.names = FALSE)
write.csv(glaucoma_3D_degs, file = glaucoma_3D_degs_file, row.names = FALSE)

# Confirm that both result files were created in the within-study output folder.
file.exists(c(glaucoma_3D_results_file, glaucoma_3D_degs_file))

# User summary confirmed Tested = 23,301, DEGs = 6, Higher = 5, Lower = 1.
# On 4 November, filesystem inspection also confirmed nonempty result exports:
# all_genes.csv = 3,443,146 bytes; DEGs_FDR005.csv = 970 bytes. Export contents
# were not independently imported or recalculated by the assistant.
# The six-gene list applies to 3-day glaucoma versus naive at BH FDR < 0.05.
# All tested-gene estimates remain available for the planned profile comparisons.

# Test the 2-week and 6-week effects against naive controls, accounting for sex.
glaucoma_2W_test <- edgeR::glmQLFTest(glaucoma_fit, coef = "groupGlaucoma_2W")
glaucoma_6W_test <- edgeR::glmQLFTest(glaucoma_fit, coef = "groupGlaucoma_6W")

# Retain every tested gene and apply BH adjustment within each comparison.
glaucoma_2W_results <- edgeR::topTags(
  glaucoma_2W_test, n = Inf, adjust.method = "BH", sort.by = "PValue"
)$table
glaucoma_6W_results <- edgeR::topTags(
  glaucoma_6W_test, n = Inf, adjust.method = "BH", sort.by = "PValue"
)$table

# Select the significant genes for each stage using the same FDR threshold.
glaucoma_2W_degs <- subset(glaucoma_2W_results, FDR < 0.05)
glaucoma_6W_degs <- subset(glaucoma_6W_results, FDR < 0.05)

# Make a count summary with one row per glaucoma stage.
glaucoma_stage_summary <- data.frame(
  Stage = c("3D", "2W", "6W"),
  Tested = c(nrow(glaucoma_3D_results), nrow(glaucoma_2W_results), nrow(glaucoma_6W_results)),
  DEGs = c(nrow(glaucoma_3D_degs), nrow(glaucoma_2W_degs), nrow(glaucoma_6W_degs)),
  Higher = c(sum(glaucoma_3D_degs$logFC > 0), sum(glaucoma_2W_degs$logFC > 0), sum(glaucoma_6W_degs$logFC > 0)),
  Lower = c(sum(glaucoma_3D_degs$logFC < 0), sum(glaucoma_2W_degs$logFC < 0), sum(glaucoma_6W_degs$logFC < 0))
)
glaucoma_stage_summary

# Set paths for the four new gene tables and the three-stage count summary.
glaucoma_stage_files <- file.path(
  glaucoma_within_dir,
  c(
    "Glaucoma_GSE241782_2W_vs_Naive_all_genes.csv",
    "Glaucoma_GSE241782_2W_vs_Naive_DEGs_FDR005.csv",
    "Glaucoma_GSE241782_6W_vs_Naive_all_genes.csv",
    "Glaucoma_GSE241782_6W_vs_Naive_DEGs_FDR005.csv",
    "Glaucoma_GSE241782_DEG_summary_by_stage.csv"
  )
)

# Save the tables to those filenames, in the same order as the list above.
write.csv(glaucoma_2W_results, file = glaucoma_stage_files[1], row.names = FALSE)
write.csv(glaucoma_2W_degs, file = glaucoma_stage_files[2], row.names = FALSE)
write.csv(glaucoma_6W_results, file = glaucoma_stage_files[3], row.names = FALSE)
write.csv(glaucoma_6W_degs, file = glaucoma_stage_files[4], row.names = FALSE)
write.csv(glaucoma_stage_summary, file = glaucoma_stage_files[5], row.names = FALSE)

# Confirm that all five output files were created.
file.exists(glaucoma_stage_files)

# User summary confirmed 23,301 tested genes at every stage: 3D has six DEGs
# (five higher, one lower), 2W has two (both higher), and 6W has three
# (one higher, two lower). Five nonempty new CSV files were found on disk.
# The assistant has not imported these CSVs or inspected the new gene lists.

# Display the 2-week DEG names, estimated changes and adjusted P values.
glaucoma_2W_degs[, c("gene_symbol", "logFC", "FDR")]

# Display the 6-week DEG names, estimated changes and adjusted P values.
glaucoma_6W_degs[, c("gene_symbol", "logFC", "FDR")]

# User output confirmed 2W DEGs: Gm28439 (logFC 16.256101, FDR 0.001602864)
# and Edn2 (4.330388, 0.016650254), both higher than naive controls.
# The 6W DEGs are Gm28439 (15.875659, 0.001950255; higher), Grk1
# (-2.181896, 0.002938366; lower), and Ak4 (-3.476917, 0.003156838; lower).
# Edn2 meets the threshold at 3D and 2W; Gm28439 at 2W and 6W. Both shared
# genes have positive estimated changes at those stages. Nine distinct genes
# meet FDR < 0.05 at one or more stages; none meets it at all three stages.
# Significance in one stage and not another does not itself establish a
# between-stage biological difference. Gene-function interpretation is pending.

# Extract original counts for the three genes with especially large estimated changes.
glaucoma_large_change_counts <- glaucoma_within_dge$counts[
  glaucoma_within_dge$genes$gene_symbol %in% c("Purb", "Sprr1a", "Gm28439"),
  ,
  drop = FALSE
]
glaucoma_large_change_counts

# Show which stage and sex correspond to each library column, in the same order.
glaucoma_within_metadata[, c("sample_id", "group", "sex")]

# User count output established the following original counts in library order:
# R-03, R-04, R-07, R-08, R-09, R-10, R-11, R-12.
# Gm28439: 25,668, 0, 0, 1, 12,456, 0, 0, 0.
# Purb:     8,968, 0, 0, 0,  5,704, 0, 0, 21,032.
# Sprr1a:      12, 125, 90, 381, 0, 0, 0, 0.
# These counts show near-zero female observations for Gm28439 and all-zero
# female observations for Purb; that does not establish sex-specific biology.
# Sprr1a has positive counts in both 3D libraries and zeros in both naive
# libraries. Library-size normalization is needed for magnitude interpretation.
# Male-only support and zero baselines constrain interpretation of the very
# large estimated changes; no automatic exclusion or model change is made here.
# The user next requested inspecting the other DEGs' counts as well.

# Collect distinct Ensembl IDs from the DEG lists for all three stages.
glaucoma_deg_ids <- unique(c(
  glaucoma_3D_degs$ensembl_id,
  glaucoma_2W_degs$ensembl_id,
  glaucoma_6W_degs$ensembl_id
))

# Extract original counts for every DEG, keeping all eight library columns.
glaucoma_deg_counts <- glaucoma_within_dge$counts[
  glaucoma_within_dge$genes$ensembl_id %in% glaucoma_deg_ids,
  ,
  drop = FALSE
]
glaucoma_deg_counts

# The user supplied the complete nine-gene count table. H3f3b, Dpm1 and
# 2610001J05Rik have higher raw counts in both 3D libraries than their same-sex
# naive controls; Edn2 has higher counts in both 2W libraries. Grk1 and Ak4
# have lower raw counts in both 6W libraries. These are raw-count observations;
# library size and composition must be allowed for when comparing magnitudes.
# The previously noted Gm28439/Purb male-concentrated patterns are confirmed.

# Calculate observed CPM using the full DGEList's saved library sizes and TMM factors.
glaucoma_cpm_within <- edgeR::cpm(
  glaucoma_within_dge,
  normalized.lib.sizes = TRUE,
  log = FALSE
)

# Select the same nine DEG rows and display their normalized measurements.
glaucoma_deg_cpm <- glaucoma_cpm_within[
  rownames(glaucoma_deg_counts),
  ,
  drop = FALSE
]
round(glaucoma_deg_cpm, 2)

# Set paths for saving both reviewed expression matrices in the analysis folder.
glaucoma_deg_expression_files <- file.path(
  glaucoma_within_dir,
  c(
    "Glaucoma_GSE241782_DEGs_raw_counts.csv",
    "Glaucoma_GSE241782_DEGs_TMM_CPM.csv"
  )
)

# Keep row names in these files because they contain the native gene labels.
write.csv(glaucoma_deg_counts, file = glaucoma_deg_expression_files[1], row.names = TRUE)
write.csv(glaucoma_deg_cpm, file = glaucoma_deg_expression_files[2], row.names = TRUE)

# Confirm that both expression-review tables were saved.
file.exists(glaucoma_deg_expression_files)

# User supplied the nine-gene CPM output and reported successful file checks.
# Figure code is in three teaching scripts with descriptive comments. Run one
# source block at a time, review its figure together, then continue to the next plot.

# Draw and save the nine-gene expression dot plot in the current R session.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "glaucoma", "03_Glaucoma_GSE241782_DEG_dotplots.R"))

# Draw and save the relative-expression heatmap after reviewing the dot plot.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "glaucoma", "04_Glaucoma_GSE241782_DEG_heatmap.R"))

# Draw and save the three MD panels after reviewing the heatmap.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "glaucoma", "05_Glaucoma_GSE241782_MD_plots.R"))

# Save the current model, stage results and plot objects at this stopping point.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "glaucoma", "06_Glaucoma_GSE241782_save_results.R"))

# Prepare and save full signed gene rankings after reviewing the DEG figures.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "glaucoma", "07_Glaucoma_GSE241782_ranked_gene_lists.R"))

# Save a final gene-level checkpoint including all three reviewed ranked tables.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "glaucoma", "09_Glaucoma_GSE241782_save_final.R"))
