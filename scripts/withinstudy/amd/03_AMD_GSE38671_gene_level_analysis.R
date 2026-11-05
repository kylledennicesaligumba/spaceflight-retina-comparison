# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "amd"), recursive = TRUE, showWarnings = FALSE)

# AMD-related CFH-deficiency model: reproducible gene-level within-study analysis.
# Authorized for autonomous execution by the user on 5 November 2026 (AMD only).
# Preserve the manual feature-level analysis and its original result files.
# Methods: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
# Representative selection: https://pmc.ncbi.nlm.nih.gov/articles/PMC3166942/

# Restore the completed precheck without changing its saved contents.
amd_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "amd")
amd_precheck_file <- Sys.getenv("AMD_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "amd", "AMD_GSE38671_precheck_checkpoint_2026-11-03_130914.RData"))
stopifnot(dir.exists(amd_within_dir), file.exists(amd_precheck_file),
          requireNamespace("limma", quietly = TRUE), requireNamespace("statmod", quietly = TRUE))
amd_precheck <- new.env()
load(amd_precheck_file, envir = amd_precheck)
amd_within_expr <- amd_precheck$amd_expr_log2
amd_within_metadata <- amd_precheck$amd_metadata_selected
amd_within_mapping <- amd_precheck$amd_annotation_clean
stopifnot(identical(dim(amd_within_expr), c(35556L, 16L)),
          all(is.finite(amd_within_expr)), anyDuplicated(rownames(amd_within_expr)) == 0L,
          identical(colnames(amd_within_expr), amd_within_metadata$geo_accession),
          all(amd_within_metadata$tissue == "neuroretina"))

# Retain all 16 arrays and reuse their existing log2 values without another transform.
amd_within_metadata$group <- factor(
  paste(amd_within_metadata$age_group, amd_within_metadata$genotype, sep = "_"),
  levels = c("young_WT", "young_Cfh-/-", "aged_WT", "aged_Cfh-/-"),
  labels = c("Young_WT", "Young_CfhKO", "Aged_WT", "Aged_CfhKO")
)
stopifnot(!anyNA(amd_within_metadata$group), all(table(amd_within_metadata$group) == 4L))

# Reproduce the reviewed annotation audit, counting each feature/gene pair once.
amd_within_mapping <- amd_within_mapping[
  amd_within_mapping$ID %in% rownames(amd_within_expr), , drop = FALSE
]
amd_missing_id <- is.na(amd_within_mapping$entrez_id) | !nzchar(trimws(amd_within_mapping$entrez_id))
amd_feature_gene_pairs <- unique(amd_within_mapping[!amd_missing_id, c("ID", "entrez_id")])
stopifnot(all(grepl("^[1-9][0-9]*$", amd_feature_gene_pairs$entrez_id)))
amd_genes_per_feature <- table(factor(amd_feature_gene_pairs$ID, levels = rownames(amd_within_expr)))
amd_mapping_audit <- data.frame(
  feature_ID = rownames(amd_within_expr),
  recorded_gene_IDs = as.integer(amd_genes_per_feature),
  truncated_annotation = rownames(amd_within_expr) %in% unique(
    amd_within_mapping$ID[amd_within_mapping$annotation_truncated %in% TRUE]),
  missing_gene_ID_assignment = rownames(amd_within_expr) %in% unique(amd_within_mapping$ID[amd_missing_id]),
  stringsAsFactors = FALSE
)

# Use the agreed conservative completeness rule; preserve excluded features in the audit.
amd_mapping_audit$eligible <- with(amd_mapping_audit,
  recorded_gene_IDs == 1L & !truncated_annotation & !missing_gene_ID_assignment)
amd_eligible_mapping <- amd_feature_gene_pairs[
  amd_feature_gene_pairs$ID %in% amd_mapping_audit$feature_ID[amd_mapping_audit$eligible], , drop = FALSE
]
amd_features_per_gene <- table(amd_eligible_mapping$entrez_id)
stopifnot(nrow(amd_eligible_mapping) == 21016L, length(amd_features_per_gene) == 19825L,
          sum(amd_features_per_gene > 1L) == 598L)

# Choose the maximum-mean representative once across all arrays, independently of significance.
# Ties use ascending feature ID with radix ordering, fixed for both age comparisons.
amd_feature_selection <- amd_eligible_mapping
amd_feature_selection$mean_log2_expression <- rowMeans(amd_within_expr[amd_feature_selection$ID, , drop = FALSE])
amd_feature_selection <- amd_feature_selection[order(
  amd_feature_selection$entrez_id, -amd_feature_selection$mean_log2_expression,
  amd_feature_selection$ID, method = "radix"), , drop = FALSE]
amd_feature_selection$selected <- !duplicated(amd_feature_selection$entrez_id)
amd_gene_representatives <- amd_feature_selection[amd_feature_selection$selected, , drop = FALSE]
amd_mapping_audit$selected_representative <- amd_mapping_audit$feature_ID %in% amd_gene_representatives$ID
amd_gene_expr <- amd_within_expr[amd_gene_representatives$ID, , drop = FALSE]
rownames(amd_gene_expr) <- amd_gene_representatives$entrez_id

# Carry supplied symbols forward as labels; preserve all distinct names rather than inventing a canonical name.
amd_symbols_by_feature <- split(amd_within_mapping$gene_symbol, amd_within_mapping$ID)
amd_gene_symbols <- vapply(amd_symbols_by_feature[amd_gene_representatives$ID], function(x) {
  paste(sort(unique(x[!is.na(x) & nzchar(trimws(x))])), collapse = " / ")
}, character(1))
amd_gene_annotation <- data.frame(
  entrez_id = amd_gene_representatives$entrez_id,
  feature_ID = amd_gene_representatives$ID,
  gene_symbol = unname(amd_gene_symbols),
  representative_mean_log2_expression = amd_gene_representatives$mean_log2_expression,
  row.names = amd_gene_representatives$entrez_id, stringsAsFactors = FALSE
)
amd_gene_annotation$gene_label <- ifelse(nzchar(amd_gene_annotation$gene_symbol),
                                      amd_gene_annotation$gene_symbol, amd_gene_annotation$entrez_id)
amd_duplicate_labels <- duplicated(amd_gene_annotation$gene_label) |
  duplicated(amd_gene_annotation$gene_label, fromLast = TRUE)
amd_gene_annotation$gene_label[amd_duplicate_labels] <- paste0(
  amd_gene_annotation$gene_label[amd_duplicate_labels], " [", amd_gene_annotation$entrez_id[amd_duplicate_labels], "]")
stopifnot(identical(dim(amd_gene_expr), c(19825L, 16L)),
          anyDuplicated(rownames(amd_gene_expr)) == 0L,
          identical(rownames(amd_gene_expr), amd_gene_annotation$entrez_id),
          identical(colnames(amd_gene_expr), amd_within_metadata$geo_accession))

# Reuse the four-group design and its two age-matched CfhKO-versus-WT contrasts.
amd_design <- model.matrix(~ group, data = amd_within_metadata)
rownames(amd_design) <- amd_within_metadata$geo_accession
stopifnot(qr(amd_design)$rank == 4L,
          identical(colnames(amd_design), c("(Intercept)", "groupYoung_CfhKO", "groupAged_WT", "groupAged_CfhKO")))
amd_contrasts <- cbind(Young_CfhKO_vs_WT = c(0, 1, 0, 0), Aged_CfhKO_vs_WT = c(0, 0, -1, 1))
rownames(amd_contrasts) <- colnames(amd_design)
amd_gene_fit <- limma::lmFit(amd_gene_expr, design = amd_design)
amd_gene_fit$genes <- amd_gene_annotation
amd_gene_contrast_fit <- limma::contrasts.fit(amd_gene_fit, contrasts = amd_contrasts)
amd_gene_contrast_fit <- limma::eBayes(amd_gene_contrast_fit, trend = TRUE, robust = TRUE)

# Extract all genes, applying BH adjustment separately within each predefined comparison.
amd_gene_results <- lapply(colnames(amd_contrasts), function(comparison) {
  limma::topTable(amd_gene_contrast_fit, coef = comparison, number = Inf, adjust.method = "BH", sort.by = "P")
})
names(amd_gene_results) <- colnames(amd_contrasts)
amd_gene_degs <- lapply(amd_gene_results, function(x) x[x$adj.P.Val < 0.05, , drop = FALSE])
amd_gene_summary <- do.call(rbind, lapply(names(amd_gene_results), function(comparison) {
  hits <- amd_gene_degs[[comparison]]
  data.frame(Comparison = comparison, Tested = nrow(amd_gene_results[[comparison]]),
             DEGs = nrow(hits), Higher = sum(hits$logFC > 0), Lower = sum(hits$logFC < 0))
}))
rownames(amd_gene_summary) <- NULL

# Save traceable mapping, sample/design keys and complete/significant gene-level tables.
write.csv(amd_mapping_audit, file.path(amd_within_dir, "AMD_GSE38671_gene_level_feature_mapping_audit.csv"), row.names = FALSE)
write.csv(amd_feature_selection, file.path(amd_within_dir, "AMD_GSE38671_gene_level_representative_selection.csv"), row.names = FALSE)
write.csv(amd_gene_annotation, file.path(amd_within_dir, "AMD_GSE38671_gene_level_annotation.csv"), row.names = FALSE)
write.csv(amd_within_metadata, file.path(amd_within_dir, "AMD_GSE38671_gene_level_sample_key.csv"), row.names = FALSE)
write.csv(amd_design, file.path(amd_within_dir, "AMD_GSE38671_gene_level_design.csv"))
write.csv(amd_contrasts, file.path(amd_within_dir, "AMD_GSE38671_gene_level_contrasts.csv"))
write.csv(amd_gene_summary, file.path(amd_within_dir, "AMD_GSE38671_gene_level_DEG_summary.csv"), row.names = FALSE)
for (comparison in names(amd_gene_results)) {
  write.csv(amd_gene_results[[comparison]], file.path(amd_within_dir, paste0("AMD_GSE38671_gene_level_", comparison, "_all_genes.csv")), row.names = FALSE)
  write.csv(amd_gene_degs[[comparison]], file.path(amd_within_dir, paste0("AMD_GSE38671_gene_level_", comparison, "_DEGs_FDR005.csv")), row.names = FALSE)
}
amd_input_checksum <- tools::md5sum(amd_precheck_file)
write.csv(data.frame(file = names(amd_input_checksum), md5 = unname(amd_input_checksum)),
          file.path(amd_within_dir, "AMD_GSE38671_gene_level_input_checksum.csv"), row.names = FALSE)
print(amd_gene_summary)
