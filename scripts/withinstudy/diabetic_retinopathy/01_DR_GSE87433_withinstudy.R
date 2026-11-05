# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "diabetic_retinopathy"), recursive = TRUE, showWarnings = FALSE)

# Diabetic retinopathy GSE87433: within-study analysis.
# Started 5 November 2026. Run one small block at a time and review its output.
# Input: deposited processed microarray expression values from the completed precheck.

# Set the folder where DR within-study outputs will be saved.
dr_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "diabetic_retinopathy")
stopifnot(dir.exists(dr_within_dir))
setwd(dr_within_dir)

# Locate the completed DR precheck checkpoint from 31 October.
dr_precheck_file <- Sys.getenv("DR_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "diabetic_retinopathy", "DR_GSE87433_precheck_checkpoint_2026-10-31.RData"))

# Confirm the output folder and check that R can find the saved precheck.
print(getwd())
print(file.exists(dr_precheck_file))

# User confirmed the correct DR output folder and TRUE checkpoint existence.

# Restore the precheck into its own container, following the previous study pattern.
dr_precheck <- new.env()
dr_loaded_objects <- load(dr_precheck_file, envir = dr_precheck)

# Make working copies of the processed expression matrix, metadata and gene annotation.
dr_within_expr <- dr_precheck$dr_expr_all
dr_within_metadata <- dr_precheck$dr_metadata_matched
dr_within_annotation <- dr_precheck$dr_gene_annotation

# Display the expression-matrix and gene-annotation dimensions.
print(dim(dr_within_expr))
print(dim(dr_within_annotation))

# Review sample identifiers, disease status, transplantation status and saved groups.
print(dr_within_metadata[, c("geo_accession", "disease", "transplantation", "group")])

# Confirm that metadata rows match the expression columns in exactly the same order.
print(identical(
  dr_within_metadata$geo_accession,
  colnames(dr_within_expr)
))

# Confirm that gene-annotation rows match the expression rows in the same order.
print(identical(
  dr_within_annotation$ID,
  rownames(dr_within_expr)
))

# User output confirmed expression 24,973 x 17, annotation 24,973 x 15,
# both order checks TRUE, and normal/diabetic/transplanted group sizes 6/6/5.
# The primary disease response compares non-transplanted diabetic and normal controls.
# Treatment context: https://doi.org/10.1007/s00125-017-4254-y

# Select the normal and diabetic controls with no transplantation.
dr_primary_metadata <- subset(
  dr_within_metadata,
  transplantation == "no" & disease %in% c("normoglycemic", "diabetic")
)

# Create the primary expression matrix in the selected metadata's sample order.
dr_primary_expr <- dr_within_expr[
  , dr_primary_metadata$geo_accession, drop = FALSE
]

# Give the two groups short labels, with Normal as the reference group.
dr_primary_metadata$group <- factor(
  dr_primary_metadata$disease,
  levels = c("normoglycemic", "diabetic"),
  labels = c("Normal", "Diabetic")
)

# Check group sizes and the number of retained gene rows and selected arrays.
print(table(dr_primary_metadata$group, useNA = "ifany"))
print(dim(dr_primary_expr))

# Confirm that selected metadata rows match the primary expression columns.
print(identical(
  dr_primary_metadata$geo_accession,
  colnames(dr_primary_expr)
))

# User output confirmed Normal 6, Diabetic 6, expression 24,973 x 12 and TRUE order.

# Build the disease-group design, using Normal as the reference group.
dr_design <- model.matrix(
  ~ group,
  data = dr_primary_metadata
)

# Label the design rows with their array identifiers and display the matrix.
rownames(dr_design) <- dr_primary_metadata$geo_accession
print(dr_design)

# Check that both model coefficients can be estimated separately.
print(qr(dr_design)$rank == ncol(dr_design))

# Calculate the residual degrees of freedom available for estimating variability.
print(nrow(dr_design) - qr(dr_design)$rank)

# User output confirmed the 12-row, two-column treatment-coded design,
# full rank TRUE and 10 residual degrees of freedom.

# Fit a linear model for every gene using the processed microarray expression values.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
dr_fit <- limma::lmFit(
  object = dr_primary_expr,
  design = dr_design
)

# Attach the already aligned gene annotation so identifiers accompany later results.
dr_fit$genes <- dr_within_annotation

# Check the number of fitted gene rows and model terms, then display term names.
print(dim(dr_fit$coefficients))
print(colnames(dr_fit$coefficients))

# User output confirmed 24,973 fitted gene rows and the two expected model terms.

# Moderate gene variances and calculate comparison statistics using limma-trend.
# Allow an expression-dependent variability trend and robust estimation across genes.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
dr_fit <- limma::eBayes(dr_fit, trend = TRUE, robust = TRUE)

# Display and save the limma diagnostic with the requested light blue variance-outlier points.
source(file.path(portfolio_repo_dir, "scripts", "withinstudy", "diabetic_retinopathy", "02_DR_GSE87433_variability_plot.R"))

# Extract all Diabetic-versus-Normal results with BH-adjusted P values.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
dr_results <- limma::topTable(
  dr_fit,
  coef = "groupDiabetic",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Select DEGs using FDR < 0.05; limma names this column adj.P.Val.
dr_degs <- subset(dr_results, adj.P.Val < 0.05)

# Count tested genes and DEGs with higher or lower expression in Diabetic.
dr_summary <- data.frame(
  Comparison = "Diabetic_vs_Normal",
  Tested = nrow(dr_results),
  DEGs = nrow(dr_degs),
  Higher = sum(dr_degs$logFC > 0),
  Lower = sum(dr_degs$logFC < 0)
)
print(dr_summary)

# Preview the ten genes with the lowest unadjusted P values.
print(head(
  dr_results[, c("Symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
  n = 10
))

# User output confirmed 24,973 tested genes and 3 DEGs: 2 higher and 1 lower.
# Samd7 and Oxct1 are higher; Rprl2 is lower. Each has adj.P.Val 0.03625559.

# Set output paths for the complete results, DEG subset and summary.
dr_result_files <- file.path(dr_within_dir, c(
  "DR_GSE87433_Diabetic_vs_Normal_all_genes.csv",
  "DR_GSE87433_Diabetic_vs_Normal_DEGs_FDR005.csv",
  "DR_GSE87433_DEG_summary.csv"
))

# Save each table; gene identifiers are already included as annotation columns.
write.csv(dr_results, dr_result_files[1], row.names = FALSE)
write.csv(dr_degs, dr_result_files[2], row.names = FALSE)
write.csv(dr_summary, dr_result_files[3], row.names = FALSE)

# Confirm that all three result files exist.
print(file.exists(dr_result_files))

# Draw and save the MD plot using the existing results and FDR < 0.05 criterion.
source(
  file.path(portfolio_repo_dir, "scripts", "withinstudy", "diabetic_retinopathy", "03_DR_GSE87433_MD_plot.R")
)

# User supplied the MD figure: most differences near zero, two higher and one lower DEG.
# Grey points beyond the DEGs illustrate that magnitude alone does not determine FDR.

# Draw and save the heatmap of all three DEGs across the 12 selected arrays.
source(
  file.path(portfolio_repo_dir, "scripts", "withinstudy", "diabetic_retinopathy", "04_DR_GSE87433_DEG_heatmap.R")
)

# Draw and save the dot plots showing each sample's processed expression for these genes.
source(
  file.path(portfolio_repo_dir, "scripts", "withinstudy", "diabetic_retinopathy", "05_DR_GSE87433_DEG_dotplots.R")
)

# User supplied both expression figures; the higher/lower patterns agree across the displays.
# These are descriptive views of the same selected DEGs and samples.

# Rank all tested DR genes using their signed moderated t statistics and save the CSV.
source(
  file.path(portfolio_repo_dir, "scripts", "withinstudy", "diabetic_retinopathy", "06_DR_GSE87433_ranked_gene_list.R")
)

# Save the completed DR analysis and a record of the R and package versions.
source(
  file.path(portfolio_repo_dir, "scripts", "withinstudy", "diabetic_retinopathy", "07_DR_GSE87433_save_final.R")
)
