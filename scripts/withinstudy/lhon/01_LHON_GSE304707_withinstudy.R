# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "lhon"), recursive = TRUE, showWarnings = FALSE)

# LHON-linked mitochondrial optic neuropathy model GSE304707: within-study analysis.
# Started 5 November 2026. Run one small block at a time and review its output.

# Set the folder where LHON within-study outputs will be saved.
lhon_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "lhon")
stopifnot(dir.exists(lhon_within_dir))
setwd(lhon_within_dir)

# Locate the completed LHON precheck checkpoint from 30 October.
lhon_precheck_file <- Sys.getenv("LHON_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "lhon", "LHON_GSE304707_precheck_checkpoint_2026-10-30.RData"))

# Confirm the output folder and check that R can find the saved precheck.
print(getwd())
print(file.exists(lhon_precheck_file))

# User confirmed the correct LHON output folder and TRUE checkpoint existence.

# Restore the precheck into its own container, following the RP/glaucoma pattern.
lhon_precheck <- new.env()
lhon_loaded_objects <- load(lhon_precheck_file, envir = lhon_precheck)

# Make working copies of the filtered/TMM count object and matched metadata.
lhon_within_dge <- lhon_precheck$lhon_dge_qc
lhon_within_metadata <- lhon_precheck$lhon_metadata_matched

# Display the number of retained gene rows and library columns.
print(dim(lhon_within_dge$counts))

# Review sample identifiers, genotypes, sexes and ages before modelling.
print(lhon_within_metadata[, c("sample_id", "genotype", "sex", "age")])

# Confirm that metadata rows match the count columns in exactly the same order.
print(identical(
  lhon_within_metadata$sample_id,
  colnames(lhon_within_dge$counts)
))

# User output confirmed 17,860 genes, 12 libraries and TRUE sample alignment.
# Each genotype has three libraries; all samples are male and 12 months old.

# Define the four genotype groups, with WT first as the reference group.
# Later comparisons will be YARS2 vs WT, COI vs WT and YARS2/COI vs WT.
# Method: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#treatments-versus-control
lhon_within_metadata$group <- factor(
  lhon_within_metadata$genotype,
  levels = c("WT", "YARS2", "COI", "YARS2/COI")
)

# Check how many libraries belong to each group, including any missing labels.
print(table(lhon_within_metadata$group, useNA = "ifany"))

# User output confirmed three libraries per group and no missing group labels.

# Build the genotype design, using WT as the reference for the other groups.
lhon_design <- model.matrix(
  ~ group,
  data = lhon_within_metadata
)

# Label the design rows with their library identifiers and display the matrix.
rownames(lhon_design) <- lhon_within_metadata$sample_id
print(lhon_design)

# Check that every model coefficient can be estimated separately.
print(qr(lhon_design)$rank == ncol(lhon_design))

# Calculate the residual degrees of freedom available for estimating variability.
print(nrow(lhon_design) - qr(lhon_design)$rank)

# User output confirmed the 12-row, four-column treatment-coded design,
# full rank TRUE and eight residual degrees of freedom.

# Fit the robust quasi-likelihood count model using saved TMM factors and the design.
# Method: https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
lhon_fit <- edgeR::glmQLFit(
  y = lhon_within_dge,
  design = lhon_design,
  robust = TRUE
)

# Check the number of fitted gene rows and model terms, then display term names.
print(dim(lhon_fit$coefficients))
print(colnames(lhon_fit$coefficients))

# User output confirmed 17,860 fitted gene rows and the four expected model terms.

# Display the gene-variability diagnostic from the fitted model.
edgeR::plotQLDisp(lhon_fit, main = "LHON: QL dispersion")

# Build the diagnostic image path inside the LHON within-study output folder.
lhon_ql_plot_file <- file.path(
  lhon_within_dir,
  "LHON_GSE304707_QL_dispersion.png"
)

# Draw the same diagnostic into a high-resolution PNG file.
png(lhon_ql_plot_file, width = 2400, height = 1800, res = 300)
edgeR::plotQLDisp(lhon_fit, main = "LHON: QL dispersion")
dev.off()

# Confirm that the diagnostic image was saved.
print(file.exists(lhon_ql_plot_file))

# The user supplied the QL diagnostic: red estimates remain broadly spread,
# with less visible contraction toward the blue trend than in RP/glaucoma.
# The extent of moderation is data-dependent; this alone is not a QC failure.

# Test YARS2 versus WT expression using the existing fitted model.
# Method: https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
lhon_YARS2_test <- edgeR::glmQLFTest(
  glmfit = lhon_fit,
  coef = "groupYARS2"
)

# Keep every tested gene and calculate BH-adjusted P values in the FDR column.
lhon_YARS2_results <- edgeR::topTags(
  lhon_YARS2_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Check the number of tested genes and preview the first ten result rows.
print(nrow(lhon_YARS2_results))
print(head(
  lhon_YARS2_results[, c("gene_name", "logFC", "logCPM", "F", "PValue", "FDR")],
  n = 10
))

# User preview shows mt-Tg, Lysmd3 and Fbn2 at FDR < 0.05, all with negative logFC.

# Select every gene meeting the agreed BH-adjusted P-value criterion.
lhon_YARS2_degs <- subset(lhon_YARS2_results, FDR < 0.05)

# Count tested genes, DEGs and their expression-change directions versus WT.
lhon_YARS2_summary <- data.frame(
  Comparison = "YARS2_vs_WT",
  Tested = nrow(lhon_YARS2_results),
  DEGs = nrow(lhon_YARS2_degs),
  Higher = sum(lhon_YARS2_degs$logFC > 0),
  Lower = sum(lhon_YARS2_degs$logFC < 0)
)
print(lhon_YARS2_summary)

# Set filenames for the complete results, significant genes and count summary.
lhon_YARS2_result_files <- file.path(
  lhon_within_dir,
  c(
    "LHON_GSE304707_YARS2_vs_WT_all_genes.csv",
    "LHON_GSE304707_YARS2_vs_WT_DEGs_FDR005.csv",
    "LHON_GSE304707_YARS2_vs_WT_DEG_summary.csv"
  )
)

# Save all three tables, retaining the supplied gene_name column.
write.csv(lhon_YARS2_results, file = lhon_YARS2_result_files[1], row.names = FALSE)
write.csv(lhon_YARS2_degs, file = lhon_YARS2_result_files[2], row.names = FALSE)
write.csv(lhon_YARS2_summary, file = lhon_YARS2_result_files[3], row.names = FALSE)

# Confirm that all three files were created in the LHON output folder.
print(file.exists(lhon_YARS2_result_files))

# User summary confirmed 17,860 tested genes and three DEGs, all lower in YARS2.

# Test COI versus WT expression using the same fitted model.
lhon_COI_test <- edgeR::glmQLFTest(
  glmfit = lhon_fit,
  coef = "groupCOI"
)

# Keep every tested gene and calculate BH-adjusted P values in the FDR column.
lhon_COI_results <- edgeR::topTags(
  lhon_COI_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Check the number of tested genes and preview the first ten result rows.
print(nrow(lhon_COI_results))
print(head(
  lhon_COI_results[, c("gene_name", "logFC", "logCPM", "F", "PValue", "FDR")],
  n = 10
))

# User preview shows mt-Tg at FDR < 0.05 with negative logFC for COI versus WT.

# Select every gene meeting the agreed BH-adjusted P-value criterion.
lhon_COI_degs <- subset(lhon_COI_results, FDR < 0.05)

# Count tested genes, DEGs and their expression-change directions versus WT.
lhon_COI_summary <- data.frame(
  Comparison = "COI_vs_WT",
  Tested = nrow(lhon_COI_results),
  DEGs = nrow(lhon_COI_degs),
  Higher = sum(lhon_COI_degs$logFC > 0),
  Lower = sum(lhon_COI_degs$logFC < 0)
)
print(lhon_COI_summary)

# Set filenames for the complete results, significant genes and count summary.
lhon_COI_result_files <- file.path(
  lhon_within_dir,
  c(
    "LHON_GSE304707_COI_vs_WT_all_genes.csv",
    "LHON_GSE304707_COI_vs_WT_DEGs_FDR005.csv",
    "LHON_GSE304707_COI_vs_WT_DEG_summary.csv"
  )
)

# Save all three tables, retaining the supplied gene_name column.
write.csv(lhon_COI_results, file = lhon_COI_result_files[1], row.names = FALSE)
write.csv(lhon_COI_degs, file = lhon_COI_result_files[2], row.names = FALSE)
write.csv(lhon_COI_summary, file = lhon_COI_result_files[3], row.names = FALSE)

# Confirm that all three files were created in the LHON output folder.
print(file.exists(lhon_COI_result_files))

# User summary confirmed 17,860 tested genes and one DEG, lower in COI versus WT.

# Test the combined YARS2/COI genotype versus WT using the same fitted model.
lhon_YARS2_COI_test <- edgeR::glmQLFTest(
  glmfit = lhon_fit,
  coef = "groupYARS2/COI"
)

# Keep every tested gene and calculate BH-adjusted P values in the FDR column.
lhon_YARS2_COI_results <- edgeR::topTags(
  lhon_YARS2_COI_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Check the number of tested genes and preview the first ten result rows.
print(nrow(lhon_YARS2_COI_results))
print(head(
  lhon_YARS2_COI_results[, c("gene_name", "logFC", "logCPM", "F", "PValue", "FDR")],
  n = 10
))

# User preview shows all ten displayed genes at FDR < 0.05 with negative logFC.

# Select every gene meeting the agreed BH-adjusted P-value criterion.
lhon_YARS2_COI_degs <- subset(lhon_YARS2_COI_results, FDR < 0.05)

# Count tested genes, DEGs and their expression-change directions versus WT.
lhon_YARS2_COI_summary <- data.frame(
  Comparison = "YARS2_COI_vs_WT",
  Tested = nrow(lhon_YARS2_COI_results),
  DEGs = nrow(lhon_YARS2_COI_degs),
  Higher = sum(lhon_YARS2_COI_degs$logFC > 0),
  Lower = sum(lhon_YARS2_COI_degs$logFC < 0)
)
print(lhon_YARS2_COI_summary)

# Set filenames for the complete results, significant genes and count summary.
lhon_YARS2_COI_result_files <- file.path(
  lhon_within_dir,
  c(
    "LHON_GSE304707_YARS2_COI_vs_WT_all_genes.csv",
    "LHON_GSE304707_YARS2_COI_vs_WT_DEGs_FDR005.csv",
    "LHON_GSE304707_YARS2_COI_vs_WT_DEG_summary.csv"
  )
)

# Save all three tables, retaining the supplied gene_name column.
write.csv(lhon_YARS2_COI_results, file = lhon_YARS2_COI_result_files[1], row.names = FALSE)
write.csv(lhon_YARS2_COI_degs, file = lhon_YARS2_COI_result_files[2], row.names = FALSE)
write.csv(lhon_YARS2_COI_summary, file = lhon_YARS2_COI_result_files[3], row.names = FALSE)

# Confirm that all three files were created in the LHON output folder.
print(file.exists(lhon_YARS2_COI_result_files))

# User summary confirmed 17,860 tested genes and 1,209 combined-genotype DEGs:
# 324 higher and 885 lower than WT at BH FDR < 0.05.

# Combine the three reviewed genotype-versus-WT summaries into one table.
lhon_genotype_summary <- rbind(
  lhon_YARS2_summary,
  lhon_COI_summary,
  lhon_YARS2_COI_summary
)
print(lhon_genotype_summary)

# Save the comparison overview in the LHON within-study output folder.
lhon_summary_file <- file.path(
  lhon_within_dir,
  "LHON_GSE304707_DEG_summary_all_genotypes.csv"
)
write.csv(lhon_genotype_summary, file = lhon_summary_file, row.names = FALSE)

# Confirm that the comparison overview was saved.
print(file.exists(lhon_summary_file))

# User output confirmed the three-comparison summary: 3, 1 and 1,209 DEGs.

# Draw and save the genome-wide MD overview for all three genotype comparisons.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "lhon", "02_LHON_GSE304707_MD_plots.R"))

# Draw and save the selected combined-genotype DEG heatmap across all four groups.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "lhon", "03_LHON_GSE304707_DEG_heatmap.R"))

# Draw and save individual-library expression plots covering the selected DEGs.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "lhon", "04_LHON_GSE304707_DEG_dotplots.R"))

# Create and save the three full gene rankings after reviewing the expression figures.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "lhon", "05_LHON_GSE304707_ranked_gene_lists.R"))

# Save the completed gene-level analysis after reviewing all three rankings.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "lhon", "06_LHON_GSE304707_save_final.R"))
