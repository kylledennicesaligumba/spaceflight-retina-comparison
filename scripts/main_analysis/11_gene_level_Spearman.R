# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Compare every terrestrial contrast with the fixed OSD-758 reference on common genes.
# Use the prospectively proposed signed-score Spearman measure; retain every contrast.
# https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html

# Calculate rank correlations using every aligned gene, with R's default average tie ranks.
main_reference_spearman <- function(score_matrix, reference_key) {
  stopifnot(
    is.matrix(score_matrix), is.numeric(score_matrix), nrow(score_matrix) >= 3L,
    !is.null(rownames(score_matrix)), !is.null(colnames(score_matrix)),
    anyDuplicated(rownames(score_matrix)) == 0L,
    anyDuplicated(colnames(score_matrix)) == 0L,
    reference_key %in% colnames(score_matrix), all(is.finite(score_matrix)),
    all(apply(score_matrix, 2, function(x) length(unique(x)) > 1L))
  )
  comparison_keys <- setdiff(colnames(score_matrix), reference_key)
  stopifnot(length(comparison_keys) > 0L)
  rho <- vapply(comparison_keys, function(key) {
    unname(stats::cor(
      score_matrix[, reference_key], score_matrix[, key],
      method = "spearman", use = "everything"
    ))
  }, numeric(1))
  stopifnot(all(is.finite(rho)), all(rho >= -1 & rho <= 1))
  rho
}

# Restore the exact aligned response checkpoint saved in the reviewed preceding block.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_gene_similarity_input_file <- file.path(
  main_analysis_dir, "Main_analysis_common_responses_checkpoint_2026-11-06.rds"
)
stopifnot(file.exists(main_gene_similarity_input_file))
main_gene_similarity_inputs <- readRDS(main_gene_similarity_input_file)
main_gene_reference_key <- "OSD758__FLT_0G_vs_HGC"
main_gene_score_matrix <- main_gene_similarity_inputs$rank_scores
main_gene_contrast_info <- main_gene_similarity_inputs$inputs
main_gene_contrast_keys <- paste(
  main_gene_contrast_info$Study, main_gene_contrast_info$Comparison, sep = "__"
)
stopifnot(
  length(main_gene_similarity_inputs$checks) > 0L,
  all(main_gene_similarity_inputs$checks),
  identical(rownames(main_gene_score_matrix), main_gene_similarity_inputs$common_ids),
  identical(colnames(main_gene_score_matrix), main_gene_contrast_keys),
  sum(main_gene_contrast_info$Role == "Reference") == 1L,
  main_gene_contrast_info$Role[match(main_gene_reference_key, main_gene_contrast_keys)] == "Reference",
  all(main_gene_contrast_info$Role[main_gene_contrast_keys != main_gene_reference_key] == "Terrestrial")
)

# Record the comparison rule before calculating the first resemblance outcomes.
main_gene_similarity_method <- list(
  measure = "Spearman rank correlation of existing signed within-study gene scores",
  reference = main_gene_reference_key,
  common_genes = nrow(main_gene_score_matrix),
  gene_basis = "Same all-13 common genes in identical order; all included regardless of FDR",
  direction = "Positive input scores denote higher expression in flight or disease versus its own controls",
  ties = "Average ranks for tied gene scores, as used by stats::cor(method = 'spearman')",
  contrast_handling = "Report every terrestrial stage/genotype/age separately; no pooling or best-stage selection",
  inference = "Descriptive point estimates; no correlation P values or confidence intervals calculated",
  interpretation = "Positive rho: agreement in signed-score ordering; negative: opposing ordering; near zero: weak overall ordering agreement",
  observed_order = "Sort rho decreasing; exactly equal rho values share an observed rank",
  planned_followup = "Review magnitudes, available-pairwise coverage sensitivity and biological programs before final conclusions",
  R_version = as.character(getRversion()),
  stats_version = as.character(utils::packageVersion("stats")),
  documentation = "https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html"
)

# Compare the complete common-gene score vector with the same fixed reference in each case.
main_common_gene_rho <- main_reference_spearman(main_gene_score_matrix, main_gene_reference_key)
main_gene_terrestrial_index <- match(names(main_common_gene_rho), main_gene_contrast_keys)
stopifnot(!anyNA(main_gene_terrestrial_index))
main_gene_similarity_results <- data.frame(
  Study = main_gene_contrast_info$Study[main_gene_terrestrial_index],
  Comparison = main_gene_contrast_info$Comparison[main_gene_terrestrial_index],
  Common_genes = nrow(main_gene_score_matrix),
  Spearman_rho = unname(main_common_gene_rho),
  row.names = NULL, stringsAsFactors = FALSE
)
main_gene_similarity_results$Observed_rank <- rank(
  -main_gene_similarity_results$Spearman_rho, ties.method = "min"
)
main_gene_similarity_results <- main_gene_similarity_results[
  order(-main_gene_similarity_results$Spearman_rho,
    main_gene_similarity_results$Study, main_gene_similarity_results$Comparison),
  c("Observed_rank", "Study", "Comparison", "Common_genes", "Spearman_rho"), drop = FALSE
]
rownames(main_gene_similarity_results) <- NULL
stopifnot(
  nrow(main_gene_similarity_results) == sum(main_gene_contrast_info$Role == "Terrestrial"),
  all(main_gene_similarity_results$Common_genes == nrow(main_gene_score_matrix))
)
print(main_gene_similarity_results, digits = 5, row.names = FALSE)

# Save full numerical precision, the prospective method record and a plain-data checkpoint.
main_gene_similarity_files <- file.path(
  main_analysis_dir,
  c(
    "Main_analysis_OSD758_gene_Spearman_common.csv",
    "Main_analysis_OSD758_gene_Spearman_method_2026-11-06.txt",
    "Main_analysis_gene_Spearman_checkpoint_2026-11-06.rds"
  )
)
utils::write.csv(main_gene_similarity_results, main_gene_similarity_files[1], row.names = FALSE)
utils::capture.output(print(main_gene_similarity_method), file = main_gene_similarity_files[2])
main_gene_similarity_checkpoint <- list(
  results = main_gene_similarity_results, rho_in_manifest_order = main_common_gene_rho,
  method = main_gene_similarity_method, common_ids = main_gene_similarity_inputs$common_ids,
  contrast_inputs = main_gene_contrast_info,
  response_checkpoint = main_gene_similarity_input_file,
  response_checkpoint_md5 = tools::md5sum(main_gene_similarity_input_file),
  original_input_md5 = main_gene_similarity_inputs$input_md5,
  annotation_versions = main_gene_similarity_inputs$annotation_versions,
  mapping_policy = main_gene_similarity_inputs$mapping_policy
)
saveRDS(main_gene_similarity_checkpoint, main_gene_similarity_files[3])
print(file.exists(main_gene_similarity_files))
