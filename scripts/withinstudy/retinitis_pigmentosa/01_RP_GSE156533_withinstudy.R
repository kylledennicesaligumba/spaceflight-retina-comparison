# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "retinitis_pigmentosa"), recursive = TRUE, showWarnings = FALSE)

# Retinitis pigmentosa GSE156533: within-study analysis.
# Started 4 November 2026. Run one small block at a time and review its output.
# Planned comparison: P23H mutant versus WT retina at P90, accounting for sex.

# Set the folder where RP within-study outputs will be saved.
rp_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "retinitis_pigmentosa")
dir.create(rp_within_dir, recursive = TRUE, showWarnings = FALSE)
setwd(rp_within_dir)

# Locate the final RP precheck checkpoint and confirm it is available.
rp_precheck_file <- Sys.getenv("RP_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "retinitis_pigmentosa", "RP_GSE156533_precheck_checkpoint_2026-10-29.RData"))
stopifnot(file.exists(rp_precheck_file))

# Restore the precheck into its own container, following the glaucoma pattern.
rp_precheck <- new.env()
rp_loaded_objects <- load(rp_precheck_file, envir = rp_precheck)

# Make working copies of the filtered count object and aligned sample metadata.
rp_within_dge <- rp_precheck$rp_dge_qc
rp_within_metadata <- rp_precheck$rp_metadata_matched

# Display the output folder and the number of gene rows and library columns.
print(getwd())
print(dim(rp_within_dge$counts))

# Review the sample identifiers, genotypes, sexes and ages before modelling.
print(rp_within_metadata[, c("sample_id", "genotype", "sex", "age")])

# Confirm that metadata rows match the count columns in exactly the same order.
print(identical(
  rp_within_metadata$sample_id,
  colnames(rp_within_dge$counts)
))

# User output confirmed 17,351 gene rows, seven libraries and TRUE alignment.
# All seven samples are P90; the reviewed genotype labels identify four mutants
# and three WT libraries, with both sexes represented in each genotype.

# Give the two genotype groups short labels, with WT as the reference group.
rp_within_metadata$group <- factor(
  rp_within_metadata$genotype,
  levels = c("RhoWT (wild-type)", "RhoP23H/WT (mutant)"),
  labels = c("WT", "P23H")
)

# Treat sex as a category, with Female as its reference, as in glaucoma.
rp_within_metadata$sex <- factor(
  rp_within_metadata$sex,
  levels = c("Female", "Male")
)

# Check the number of libraries in each genotype and sex combination.
print(with(rp_within_metadata, table(group, sex, useNA = "ifany")))

# User confirmed WT has two female and one male library; P23H has two of each.

# Build an additive design for the genotype comparison while accounting for sex.
# Method: edgeR User's Guide, additive models and blocking.
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
rp_design <- model.matrix(
  ~ sex + group,
  data = rp_within_metadata
)

# Label the design rows with their library identifiers and display the matrix.
rownames(rp_design) <- rp_within_metadata$sample_id
print(rp_design)

# Check that every model coefficient can be estimated separately.
print(qr(rp_design)$rank == ncol(rp_design))

# Calculate the residual degrees of freedom available for estimating variability.
print(nrow(rp_design) - qr(rp_design)$rank)

# User output confirmed the seven-row, three-column treatment-coded design,
# full rank TRUE and four residual degrees of freedom.

# Fit the robust quasi-likelihood count model using saved TMM factors and the design.
# This follows the same current edgeR pipeline used for glaucoma (edgeR 4.8.2).
# Method: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
rp_fit <- edgeR::glmQLFit(
  y = rp_within_dge,
  design = rp_design,
  robust = TRUE
)

# Check the number of fitted gene rows and model terms, then display term names.
print(dim(rp_fit$coefficients))
print(colnames(rp_fit$coefficients))

# User output confirmed 17,351 fitted gene rows and the three expected model terms.

# Display the gene-variability diagnostic from the fitted model.
edgeR::plotQLDisp(rp_fit, main = "RP: QL dispersion")

# Build the diagnostic image path inside the RP within-study output folder.
rp_ql_plot_file <- file.path(
  rp_within_dir,
  "RP_GSE156533_QL_dispersion.png"
)

# Draw the same diagnostic into a high-resolution PNG file.
png(rp_ql_plot_file, width = 2400, height = 1800, res = 300)
edgeR::plotQLDisp(rp_fit, main = "RP: QL dispersion")
dev.off()

# Confirm that the diagnostic image was saved.
print(file.exists(rp_ql_plot_file))

# The user supplied the saved QL diagnostic. Red estimates form a tighter band
# around the trend than the initial black estimates; some high estimates remain.
# The blue trend decreases and then rises slightly at the highest abundances.
# This shows the moderation pattern, not gene-level significance or model validation.

# Test P23H mutant versus WT expression while accounting for sex.
rp_test <- edgeR::glmQLFTest(
  glmfit = rp_fit,
  coef = "groupP23H"
)

# Keep every tested gene and calculate BH-adjusted P values in the FDR column.
rp_results <- edgeR::topTags(
  rp_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Check the number of tested genes and preview their names and statistical results.
# RP's supplied gene-symbol field is gene_name; all annotations remain in rp_results.
print(nrow(rp_results))
print(head(
  rp_results[, c("gene_name", "logFC", "logCPM", "F", "PValue", "FDR")],
  n = 10
))

# The user supplied the first ten results: all meet FDR < 0.05.
# Nine have positive logFC; Glb1l3 has negative logFC. The total DEG count
# cannot be inferred from this preview alone.

# Select all genes meeting the agreed BH-adjusted P-value criterion.
rp_degs <- subset(rp_results, FDR < 0.05)

# Count all tested genes, DEGs and their expression-change directions versus WT.
rp_summary <- data.frame(
  Comparison = "P23H_vs_WT",
  Tested = nrow(rp_results),
  DEGs = nrow(rp_degs),
  Higher = sum(rp_degs$logFC > 0),
  Lower = sum(rp_degs$logFC < 0)
)
print(rp_summary)

# Set filenames for the complete results, significant genes and count summary.
rp_result_files <- file.path(
  rp_within_dir,
  c(
    "RP_GSE156533_P23H_vs_WT_all_genes.csv",
    "RP_GSE156533_P23H_vs_WT_DEGs_FDR005.csv",
    "RP_GSE156533_DEG_summary.csv"
  )
)

# Save all three tables; gene identifiers are retained in the gene_id column.
write.csv(rp_results, file = rp_result_files[1], row.names = FALSE)
write.csv(rp_degs, file = rp_result_files[2], row.names = FALSE)
write.csv(rp_summary, file = rp_result_files[3], row.names = FALSE)

# Confirm that the three result files were created in the RP output folder.
print(file.exists(rp_result_files))

# User summary confirmed 17,351 tested genes and 2,576 DEGs: 1,361 higher
# and 1,215 lower in P23H mutants versus WT at BH FDR < 0.05.

# Draw and save the genome-wide MD overview after reviewing the DEG summary.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "retinitis_pigmentosa", "02_RP_GSE156533_MD_plot.R"))

# Draw and save a readable heatmap of 15 higher and 15 lower DEGs.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "retinitis_pigmentosa", "03_RP_GSE156533_DEG_heatmap.R"))

# Draw and save individual-library expression plots for five higher and five lower DEGs.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "retinitis_pigmentosa", "04_RP_GSE156533_DEG_dotplots.R"))

# Prepare and save the full signed gene ranking after reviewing the result figures.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "retinitis_pigmentosa", "05_RP_GSE156533_ranked_gene_list.R"))

# Save the completed gene-level analysis after reviewing the full ranking.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "retinitis_pigmentosa", "06_RP_GSE156533_save_final.R"))
