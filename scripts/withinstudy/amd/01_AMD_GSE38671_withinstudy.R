# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "amd"), recursive = TRUE, showWarnings = FALSE)

# AMD-related CFH-deficiency model GSE38671: within-study analysis.
# Started 5 November 2026. Run one small commented block at a time and review its output.
# Preserve the selected neural-retina scope and separate young/aged genotype comparisons.

# Set the folder where AMD within-study scripts and outputs will be saved.
amd_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "amd")
dir.create(amd_within_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(amd_within_dir))
setwd(amd_within_dir)

# Point to the latest saved precheck identified by the existing checkpoint pointer.
amd_precheck_file <- Sys.getenv("AMD_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "amd", "AMD_GSE38671_precheck_checkpoint_2026-11-03_130914.RData"))

# Confirm the working folder and check that R can find the saved precheck.
print(getwd())
print(file.exists(amd_precheck_file))

# User confirmed the AMD within-study folder and TRUE checkpoint existence.

# Restore the saved precheck into its own container, following the previous study pattern.
amd_precheck <- new.env()
amd_loaded_objects <- load(amd_precheck_file, envir = amd_precheck)

# Make working copies of the saved log2 neural-retina expression, metadata and annotation.
amd_within_expr <- amd_precheck$amd_expr_log2
amd_within_metadata <- amd_precheck$amd_metadata_selected
amd_within_annotation <- amd_precheck$amd_feature_annotation

# Display the expression-matrix and feature-annotation dimensions.
print(dim(amd_within_expr))
print(dim(amd_within_annotation))

# Review the selected sample identifiers, ages, genotypes and tissue before modelling.
print(amd_within_metadata[, c("geo_accession", "age_group", "genotype", "tissue")])

# Confirm that metadata rows match the expression columns in exactly the same order.
print(identical(
  amd_within_metadata$geo_accession,
  colnames(amd_within_expr)
))

# Confirm that annotation rows match the expression feature IDs in the same order.
print(identical(
  amd_within_annotation$ID,
  rownames(amd_within_expr)
))

# User output confirmed expression 35,556 x 16, annotation 35,556 x 12,
# all selected samples neuroretina, four per age/genotype group and both order checks TRUE.

# Give the four age/genotype groups short labels and set their order for modelling.
amd_within_metadata$group <- factor(
  paste(amd_within_metadata$age_group, amd_within_metadata$genotype, sep = "_"),
  levels = c("young_WT", "young_Cfh-/-", "aged_WT", "aged_Cfh-/-"),
  labels = c("Young_WT", "Young_CfhKO", "Aged_WT", "Aged_CfhKO")
)

# Count the samples in each group and show any unmatched labels.
print(table(amd_within_metadata$group, useNA = "ifany"))
stopifnot(!anyNA(amd_within_metadata$group))

# User output confirmed four samples in each of the four age/genotype groups.

# Build the four-group design, using Young_WT as the reference group.
# Method: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
amd_design <- model.matrix(
  ~ group,
  data = amd_within_metadata
)

# Label the design rows with their array identifiers and display the matrix.
rownames(amd_design) <- amd_within_metadata$geo_accession
print(amd_design)

# Check that every model coefficient can be estimated separately.
print(qr(amd_design)$rank == ncol(amd_design))

# Calculate the residual degrees of freedom available for estimating variability.
print(nrow(amd_design) - qr(amd_design)$rank)

# The later aged genotype contrast must compare Aged_CfhKO with Aged_WT.
# In this reference-coded design, subtract groupAged_WT from groupAged_CfhKO.

# User output confirmed 16 named rows, four expected columns, treatment coding,
# full rank TRUE and 12 residual degrees of freedom.

# Fit a linear model for each array feature using the saved log2 expression values.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
amd_fit <- limma::lmFit(
  object = amd_within_expr,
  design = amd_design
)

# Attach the already aligned feature annotation so identifiers accompany later results.
amd_fit$genes <- amd_within_annotation

# Check the number of fitted feature rows and model terms, then display term names.
print(dim(amd_fit$coefficients))
print(colnames(amd_fit$coefficients))

# User output confirmed 35,556 fitted feature rows and the four expected model terms.

# Define the young and aged CfhKO-versus-WT comparisons in the confirmed coefficient order.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
amd_contrasts <- cbind(
  Young_CfhKO_vs_WT = c(0, 1, 0, 0),
  Aged_CfhKO_vs_WT = c(0, 0, -1, 1)
)

# Label each row with its fitted coefficient and display the comparison instructions.
rownames(amd_contrasts) <- colnames(amd_fit$coefficients)
print(amd_contrasts)

# Calculate both age-matched expression differences from the existing fitted model.
amd_contrast_fit <- limma::contrasts.fit(
  fit = amd_fit,
  contrasts = amd_contrasts
)

# Confirm that the comparison fit contains every feature and the two intended comparisons.
print(dim(amd_contrast_fit$coefficients))
print(colnames(amd_contrast_fit$coefficients))

# User output confirmed the contrast instructions and 35,556 x 2 comparison coefficients.

# Moderate feature variances and calculate statistics for both comparisons, following DR.
# Allow expression-dependent variability and robust estimation across features.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
amd_contrast_fit <- limma::eBayes(
  amd_contrast_fit,
  trend = TRUE,
  robust = TRUE
)

# Display and save the limma mean-variance diagnostic using the chosen blue outlier colour.
source(
  file.path(portfolio_repo_dir, "scripts", "withinstudy", "amd", "02_AMD_GSE38671_variability_plot.R")
)

# User supplied the SA diagnostic and asks whether any variance outliers are highlighted.
# Reproduce installed limma 3.66.0 plotSA's colour rule for this unweighted, positive-df fit.
amd_sa_variance_ratio <- amd_contrast_fit$sigma^2 / amd_contrast_fit$s2.prior
amd_sa_tail_p <- 2 * pmin(
  stats::pf(
    amd_sa_variance_ratio, df1 = amd_contrast_fit$df.residual,
    df2 = max(amd_contrast_fit$df.prior)
  ),
  stats::pf(
    amd_sa_variance_ratio, df1 = amd_contrast_fit$df.residual,
    df2 = max(amd_contrast_fit$df.prior), lower.tail = FALSE
  )
)

# Match plotSA's internal variance-display cutoff of 0.5; the DEG cutoff remains 0.05.
amd_sa_variance_outlier_fdr <- stats::p.adjust(amd_sa_tail_p, method = "BH")
amd_sa_variance_outliers <- amd_sa_variance_outlier_fdr <= 0.5

# Count the features assigned the outlier colour under that existing diagnostic rule.
print(c(
  Highlighted_variance_outliers = sum(amd_sa_variance_outliers, na.rm = TRUE)
))

# User output confirmed zero highlighted feature-variance outliers under plotSA's rule.

# Extract every feature's young CfhKO-versus-WT result with BH-adjusted P values.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
amd_young_results <- limma::topTable(
  amd_contrast_fit,
  coef = "Young_CfhKO_vs_WT",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Extract every feature's aged CfhKO-versus-WT result with BH-adjusted P values.
amd_aged_results <- limma::topTable(
  amd_contrast_fit,
  coef = "Aged_CfhKO_vs_WT",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Preview the ten lowest-P feature results from the young comparison.
print(head(
  amd_young_results[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
  n = 10
))

# Preview the ten lowest-P feature results from the aged comparison.
print(head(
  amd_aged_results[, c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
  n = 10
))

# User reviewed both previews; passing rows are significant array features pending gene mapping.

# Select significant features from each complete result table using FDR < 0.05.
amd_young_de_features <- subset(amd_young_results, adj.P.Val < 0.05)
amd_aged_de_features <- subset(amd_aged_results, adj.P.Val < 0.05)

# Count tested features and significant features with higher or lower expression in CfhKO.
amd_feature_summary <- data.frame(
  Comparison = c("Young_CfhKO_vs_WT", "Aged_CfhKO_vs_WT"),
  Tested_features = c(nrow(amd_young_results), nrow(amd_aged_results)),
  Significant_features = c(nrow(amd_young_de_features), nrow(amd_aged_de_features)),
  Higher = c(
    sum(amd_young_de_features$logFC > 0),
    sum(amd_aged_de_features$logFC > 0)
  ),
  Lower = c(
    sum(amd_young_de_features$logFC < 0),
    sum(amd_aged_de_features$logFC < 0)
  )
)
print(amd_feature_summary)

# Set paths for both complete tables, both significant-feature lists and the summary.
amd_feature_result_files <- file.path(amd_within_dir, c(
  "AMD_GSE38671_Young_CfhKO_vs_WT_all_features.csv",
  "AMD_GSE38671_Young_CfhKO_vs_WT_significant_features_FDR005.csv",
  "AMD_GSE38671_Aged_CfhKO_vs_WT_all_features.csv",
  "AMD_GSE38671_Aged_CfhKO_vs_WT_significant_features_FDR005.csv",
  "AMD_GSE38671_significant_feature_summary.csv"
))

# Save the results with their supplied feature IDs and annotation columns preserved.
write.csv(amd_young_results, amd_feature_result_files[1], row.names = FALSE)
write.csv(amd_young_de_features, amd_feature_result_files[2], row.names = FALSE)
write.csv(amd_aged_results, amd_feature_result_files[3], row.names = FALSE)
write.csv(amd_aged_de_features, amd_feature_result_files[4], row.names = FALSE)
write.csv(amd_feature_summary, amd_feature_result_files[5], row.names = FALSE)

# Confirm that all five result files exist in the AMD within-study folder.
print(file.exists(amd_feature_result_files))

# User confirmed five significant young features, 51 aged features and all five saved tables.
# Begin the requested gene-mapping approach using every measured feature.

# Reuse the readable annotation assignments saved during the AMD precheck.
amd_within_mapping <- amd_precheck$amd_annotation_clean
stopifnot(is.data.frame(amd_within_mapping))

# Keep annotation records for features present in the current expression matrix.
amd_within_mapping <- amd_within_mapping[
  amd_within_mapping$ID %in% rownames(amd_within_expr), , drop = FALSE
]

# Keep one copy of each recorded feature-gene pair, merging repeated transcript assignments.
amd_feature_gene_pairs <- unique(amd_within_mapping[
  !is.na(amd_within_mapping$entrez_id) &
    nzchar(trimws(amd_within_mapping$entrez_id)),
  c("ID", "entrez_id"), drop = FALSE
])

# Check the recorded gene-ID format before counting; database validity is a later check.
stopifnot(all(grepl("^[1-9][0-9]*$", amd_feature_gene_pairs$entrez_id)))

# Count distinct recorded gene IDs per feature, including features with no recorded gene ID.
amd_genes_per_feature <- table(factor(
  amd_feature_gene_pairs$ID,
  levels = rownames(amd_within_expr)
))
print(table(Recorded_gene_IDs_per_feature = as.integer(amd_genes_per_feature)))

# Report coverage and truncation flags; one recorded ID alone does not prove a complete mapping.
print(c(
  Measured_features = nrow(amd_within_expr),
  Features_in_mapping = length(unique(amd_within_mapping$ID)),
  Features_with_truncated_annotation = length(unique(
    amd_within_mapping$ID[amd_within_mapping$annotation_truncated]
  ))
))

# User confirmed 10,158 features with no recorded gene, 22,610 with one and 2,788 with several.
# Identify features with exactly one distinct recorded gene ID.
amd_single_gene_ids <- names(amd_genes_per_feature)[amd_genes_per_feature == 1L]
amd_single_gene_mapping <- amd_feature_gene_pairs[
  amd_feature_gene_pairs$ID %in% amd_single_gene_ids, , drop = FALSE
]

# Identify features whose supplied annotation text was cut short.
amd_truncated_feature_ids <- unique(amd_within_mapping$ID[
  amd_within_mapping$annotation_truncated %in% TRUE
])

# Identify features with at least one annotation assignment lacking a gene ID.
amd_missing_gene_id_features <- unique(amd_within_mapping$ID[
  is.na(amd_within_mapping$entrez_id) |
    !nzchar(trimws(amd_within_mapping$entrez_id))
])

# Mark these two possible problems in the single-gene mapping table.
amd_single_gene_mapping$truncated_annotation <-
  amd_single_gene_mapping$ID %in% amd_truncated_feature_ids
amd_single_gene_mapping$missing_gene_ID_assignment <-
  amd_single_gene_mapping$ID %in% amd_missing_gene_id_features

# Count single-gene features with each flag and those with neither flag; flags can overlap.
print(c(
  Single_gene_features = nrow(amd_single_gene_mapping),
  With_truncated_annotation = sum(amd_single_gene_mapping$truncated_annotation),
  With_missing_gene_ID_assignments = sum(amd_single_gene_mapping$missing_gene_ID_assignment),
  Without_either_flag = sum(
    !amd_single_gene_mapping$truncated_annotation &
      !amd_single_gene_mapping$missing_gene_ID_assignment
  )
))

# User confirmed 22,610 single-ID features: zero truncated, 1,594 with missing-ID records, 21,016 with neither flag.
# Choose three incomplete examples so we can review their recorded gene assignments.
amd_incomplete_example_ids <- head(
  amd_single_gene_mapping$ID[
    amd_single_gene_mapping$missing_gene_ID_assignment
  ],
  n = 3
)

# Show distinct gene-label/ID combinations, including records with missing IDs.
amd_incomplete_examples <- unique(amd_within_mapping[
  amd_within_mapping$ID %in% amd_incomplete_example_ids,
  c("ID", "gene_symbol", "entrez_id"), drop = FALSE
])

# Display up to 15 combinations; these examples do not establish all 1,594 mappings.
print(head(amd_incomplete_examples, n = 15), row.names = FALSE)

# User reviewed one feature with several gene names and two with repeated names but missing IDs.
# For the main gene analysis, require one recorded gene ID with no incomplete/truncated assignments.
# This conservative project rule can omit usable features until their annotation is independently resolved.
# Keep the original feature-level expression, annotation and statistical results intact.
amd_eligible_mapping <- amd_single_gene_mapping[
  !amd_single_gene_mapping$truncated_annotation &
    !amd_single_gene_mapping$missing_gene_ID_assignment,
  c("ID", "entrez_id"), drop = FALSE
]

# Count how many eligible array features represent each distinct recorded gene ID.
amd_features_per_gene <- table(amd_eligible_mapping$entrez_id)

# Distinguish feature totals from gene totals before selecting one representative per gene.
print(c(
  Eligible_array_features = nrow(amd_eligible_mapping),
  Distinct_recorded_genes = length(amd_features_per_gene),
  Genes_with_one_feature = sum(amd_features_per_gene == 1L),
  Genes_with_multiple_features = sum(amd_features_per_gene > 1L)
))

# User confirmed 21,016 eligible features, 19,825 gene IDs, and 598 genes with multiple features.
# Copy the eligible mapping and calculate each feature's mean log2 expression across all 16 arrays.
# Representative-selection method: https://pmc.ncbi.nlm.nih.gov/articles/PMC3166942/
amd_feature_selection <- amd_eligible_mapping
amd_feature_selection$mean_log2_expression <- rowMeans(
  amd_within_expr[amd_feature_selection$ID, , drop = FALSE]
)

# Put each gene's highest-mean feature first; break exact ties by ascending feature ID.
amd_feature_selection <- amd_feature_selection[order(
  amd_feature_selection$entrez_id,
  -amd_feature_selection$mean_log2_expression,
  amd_feature_selection$ID,
  method = "radix"
), , drop = FALSE]

# Keep one representative per gene, using the same selection for both age comparisons.
amd_gene_representatives <- amd_feature_selection[
  !duplicated(amd_feature_selection$entrez_id), , drop = FALSE
]

# Copy the selected expression rows and label them with their recorded Entrez gene IDs.
amd_gene_expr <- amd_within_expr[
  amd_gene_representatives$ID, , drop = FALSE
]
rownames(amd_gene_expr) <- amd_gene_representatives$entrez_id

# Check the gene/sample dimensions, unique gene IDs and unchanged sample-column order.
print(dim(amd_gene_expr))
print(anyDuplicated(rownames(amd_gene_expr)) == 0L)
print(identical(
  amd_within_metadata$geo_accession,
  colnames(amd_gene_expr)
))
