# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# OSD-758: manual within-study analysis, started 6 November 2026.
# Run one small block at a time and review the user's output before advancing.
# After a restart, use 02_OSD758_withinstudy_reload.R for the input setup,
# or 03_OSD758_reload_model.R after the model checkpoint has been saved.
# Primary comparison: flight 0G versus habitat ground control; retain all 17
# reviewed libraries, including HGC05 with its recorded projection-review flag.
# Planned model: ~ group, with HGC as reference, using robust edgeR direct QL.
# Kremsky et al. (2024), Methods 4.1-4.2, reports male C57BL/6J mice from the
# same cohort and individually preserved right-retinal tissue. HGC procedures
# ran three days behind flight, including euthanasia after the return of flight
# mice. Interpret the contrast as an overall flight-associated retinal response.
# Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC11593819/
# Existing precheck filtering/TMM are carried forward rather than repeated.

# Set a separate folder for the OSD-758 within-study scripts and outputs.
osd758_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758")
dir.create(osd758_within_dir, recursive = TRUE, showWarnings = FALSE)

# Locate the saved precheck checkpoint and confirm that it is available.
osd758_precheck_file <- Sys.getenv("OSD758_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "osdr758", "OSD758_precheck_checkpoint_2026-11-06.RData"))
stopifnot(file.exists(osd758_precheck_file))

# Restore the precheck into a separate environment, as in the RP workflow.
osd758_precheck <- new.env()
osd758_loaded_objects <- load(osd758_precheck_file, envir = osd758_precheck)

# Make within-study working copies of the filtered/TMM counts and sample metadata.
osd758_within_dge <- osd758_precheck$osd758_dge_qc
osd758_within_metadata <- osd758_precheck$osd758_metadata_selected

# Display the input dimensions and numbers of libraries in each condition.
dim(osd758_within_dge$counts)
table(osd758_within_metadata$condition)

# Confirm that metadata rows match the count columns in exactly the same order.
identical(
  colnames(osd758_within_dge$counts), osd758_within_metadata$sample
)

# User confirmed the restored input: 21,335 x 17, FLT_0G 5 / HGC 12,
# and TRUE sample alignment. The edgeR namespace-loading message is normal.

# Set habitat ground control as the reference for the flight comparison.
osd758_within_metadata$group <- factor(
  osd758_within_metadata$condition, levels = c("HGC", "FLT_0G")
)

# Build a group-only design; each row represents one selected sample.
osd758_design <- model.matrix(~ group, data = osd758_within_metadata)
rownames(osd758_design) <- osd758_within_metadata$sample
osd758_design

# Check that both model terms are estimable and sample order matches the counts.
c(
  full_rank = qr(osd758_design)$rank == ncol(osd758_design),
  sample_order = identical(rownames(osd758_design), colnames(osd758_within_dge$counts))
)

# Calculate nominal residual degrees of freedom: samples minus model rank.
nrow(osd758_design) - qr(osd758_design)$rank

# User confirmed treatment coding, full_rank TRUE, sample_order TRUE and
# nominal residual df 15. The design has (Intercept) and groupFLT_0G terms.

# Fit the robust quasi-likelihood model using the saved counts and TMM factors.
osd758_fit <- edgeR::glmQLFit(
  y = osd758_within_dge,
  design = osd758_design,
  robust = TRUE
)

# Check how many gene rows and model terms were fitted.
dim(osd758_fit$coefficients)
colnames(osd758_fit$coefficients)

# Save the fitted model and its inputs so this stage can be restored directly.
osd758_model_checkpoint_file <- file.path(
  osd758_within_dir, "OSD758_withinstudy_model_checkpoint_2026-11-06.RData"
)
save(
  osd758_within_dge, osd758_within_metadata,
  osd758_design, osd758_fit,
  osd758_within_dir, osd758_precheck_file,
  file = osd758_model_checkpoint_file
)
file.exists(osd758_model_checkpoint_file)

# User confirmed 21,335 x 2 fitted coefficients, the expected (Intercept) /
# groupFLT_0G terms and TRUE model-checkpoint file presence on 6 November 2026.

# Display the gene-wise variability diagnostic from the fitted QL model.
edgeR::plotQLDisp(osd758_fit, main = "OSD-758: QL dispersion")

# Set the diagnostic filename in the within-study output folder.
osd758_ql_plot_file <- file.path(osd758_within_dir, "OSD758_QL_dispersion.png")

# Save the same diagnostic as a high-resolution PNG and close the device.
png(osd758_ql_plot_file, width = 2400, height = 1800, res = 300)
edgeR::plotQLDisp(osd758_fit, main = "OSD-758: QL dispersion")
dev.off()
file.exists(osd758_ql_plot_file)

# User supplied the saved QL diagnostic and requested interpretation on 6 November.
# Black: unmoderated gene-wise variability; red: moderated estimates; blue: trend.
# The red cloud is generally less dispersed, while some high-variability genes
# remain well above the trend. This is partial moderation, not convergence of
# every estimate to the blue line. The trend rises slightly at low abundance,
# declines around log2 CPM 5-7, then rises modestly at the highest abundances.
# The y-axis is a quarter-root mean-deviance scale, not an effect-size scale.
# These observations do not establish DE genes, a failed sample, or complete
# model validity. Retain the current input and proceed with the planned QL test
# after this conceptual interpretation; no additional filtering is introduced.

# Test flight 0G versus HGC for every fitted gene using the confirmed coefficient.
osd758_test <- edgeR::glmQLFTest(
  glmfit = osd758_fit,
  coef = "groupFLT_0G"
)

# Keep every tested gene and calculate BH-adjusted P values in the FDR column.
osd758_results <- edgeR::topTags(
  osd758_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Preserve the original source gene identifiers as an explicit result column.
osd758_results$gene_id <- rownames(osd758_results)

# Display the number of tested genes and the ten lowest-P results.
nrow(osd758_results)
head(
  osd758_results[, c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")],
  n = 10
)

# User confirmed 21,335 tested rows. All ten lowest-P preview rows meet
# FDR < 0.05: nine higher in flight and ENSMUSG00000020288 lower (logFC -0.941382).
# This preview does not establish the total number of DEGs.

# Select genes meeting the agreed FDR criterion without an extra fold-change cutoff.
osd758_degs <- subset(osd758_results, FDR < 0.05)

# Summarize all tested genes, significant genes and directions relative to HGC.
osd758_summary <- data.frame(
  Comparison = "FLT_0G_vs_HGC",
  Tested = nrow(osd758_results),
  DEGs = nrow(osd758_degs),
  Higher = sum(osd758_degs$logFC > 0),
  Lower = sum(osd758_degs$logFC < 0)
)
osd758_summary

# Set filenames for the full results, DEG list and count summary.
osd758_result_files <- file.path(
  osd758_within_dir,
  c(
    "OSD758_FLT_0G_vs_HGC_all_genes.csv",
    "OSD758_FLT_0G_vs_HGC_DEGs_FDR005.csv",
    "OSD758_DEG_summary.csv"
  )
)

# Save all three tables with the original gene IDs preserved in gene_id.
write.csv(osd758_results, osd758_result_files[1], row.names = FALSE)
write.csv(osd758_degs, osd758_result_files[2], row.names = FALSE)
write.csv(osd758_summary, osd758_result_files[3], row.names = FALSE)
file.exists(osd758_result_files)

# User confirmed 21,335 tested genes, 373 DEGs (286 higher and 87 lower in flight),
# and three TRUE checks for the all-gene, DEG and summary CSV files.

# Draw and save the genome-wide MD plot using the reviewed results.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "osdr758", "04_OSD758_MD_plot.R"))

# User supplied the saved MD plot: most genes are near zero, with more significant
# increases than decreases. Some large changes are grey because they fail FDR < 0.05.
# The plot does not establish pathways or an overall increase in total cellular RNA.

# Draw and save relative expression for 15 higher and 15 lower DEGs.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "osdr758", "05_OSD758_DEG_heatmap.R"))

# Draw and save individual-library expression for five DEGs per direction.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "osdr758", "06_OSD758_DEG_dotplots.R"))

# User supplied both expression figures: the selected genes show clear directional
# group differences with within-group variation. These are selected DEG displays,
# not independent validation or evidence that all tested genes separate the groups.
# Supporting CSV existence checks were not shown in this message.

# Rank and save every tested gene using the same signed score as RP/glaucoma.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "osdr758", "07_OSD758_ranked_gene_list.R"))

# User confirmed 21,335 ranked genes and a TRUE export check. The original F values
# have minimum 1.142974e-08, no negatives and no nonfinite values, so none were clipped.
# The shown scores span 10.716504 at the positive end to -7.611517 at the negative end.

# Save the completed within-study analysis and its R/package version record.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "osdr758", "08_OSD758_save_final.R"))
