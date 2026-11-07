# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Compare the same pathway NES profiles against the fixed OSD-758 reference.
# Retain every pathway for correlation; use FDR flags separately for interpretation.
# https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html

# Align one complete result table by pathway ID and reject incomplete estimates.
main_align_pathway_result <- function(result, pathway_ids) {
  stopifnot(
    all(c("Study", "Comparison", "Role", "Pathway_ID", "Pathway_name", "size",
      "NES", "FDR", "Estimation_status") %in% names(result)),
    is.character(pathway_ids), length(pathway_ids) > 1L,
    !anyNA(pathway_ids), anyDuplicated(pathway_ids) == 0L,
    !anyNA(result$Pathway_ID), anyDuplicated(result$Pathway_ID) == 0L,
    nrow(result) == length(pathway_ids), setequal(result$Pathway_ID, pathway_ids),
    is.numeric(result$NES), all(is.finite(result$NES)),
    is.numeric(result$FDR), all(is.finite(result$FDR) & result$FDR >= 0 & result$FDR <= 1),
    all(result$Estimation_status == "Estimated")
  )
  result <- result[match(pathway_ids, result$Pathway_ID), , drop = FALSE]
  rownames(result) <- NULL
  result
}

# Calculate descriptive reference correlations without dropping or FDR-selecting pathways.
main_pathway_reference_spearman <- function(nes_matrix, reference_key) {
  stopifnot(
    is.matrix(nes_matrix), is.numeric(nes_matrix), nrow(nes_matrix) >= 3L,
    !is.null(rownames(nes_matrix)), !is.null(colnames(nes_matrix)),
    anyDuplicated(rownames(nes_matrix)) == 0L,
    anyDuplicated(colnames(nes_matrix)) == 0L,
    reference_key %in% colnames(nes_matrix), all(is.finite(nes_matrix)),
    all(apply(nes_matrix, 2, function(x) length(unique(x)) > 1L))
  )
  other_keys <- setdiff(colnames(nes_matrix), reference_key)
  stopifnot(length(other_keys) > 0L)
  rho <- vapply(other_keys, function(key) unname(stats::cor(
    nes_matrix[, reference_key], nes_matrix[, key],
    method = "spearman", use = "everything"
  )), numeric(1))
  stopifnot(all(is.finite(rho)), all(rho >= -1 & rho <= 1))
  rho
}

# Restore the complete enrichment results reviewed by the user, without rerunning enrichment.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_pathway_comparison_input <- file.path(
  main_analysis_dir, "Main_analysis_Reactome_fgsea_checkpoint_2026-11-06.rds"
)
stopifnot(file.exists(main_pathway_comparison_input))
main_pathway_comparison_context <- readRDS(main_pathway_comparison_input)
main_pathway_profile_ids <- main_pathway_comparison_context$pathway_ids
main_pathway_profile_info <- main_pathway_comparison_context$contrasts
main_pathway_profile_keys <- paste(
  main_pathway_profile_info$Study, main_pathway_profile_info$Comparison, sep = "__"
)
main_pathway_reference_key <- "OSD758__FLT_0G_vs_HGC"
stopifnot(
  setequal(names(main_pathway_comparison_context$results), main_pathway_profile_keys),
  anyDuplicated(main_pathway_profile_keys) == 0L,
  main_pathway_reference_key %in% main_pathway_profile_keys,
  sum(main_pathway_profile_info$Role == "Reference") == 1L,
  main_pathway_profile_info$Role[match(main_pathway_reference_key, main_pathway_profile_keys)] == "Reference",
  all(main_pathway_profile_info$Role[main_pathway_profile_keys != main_pathway_reference_key] == "Terrestrial"),
  is.character(main_pathway_comparison_context$method$planned_profile_comparison),
  nzchar(main_pathway_comparison_context$method$planned_profile_comparison),
  main_pathway_comparison_context$method$FDR_cutoff == 0.05
)
main_pathway_aligned_results <- stats::setNames(lapply(main_pathway_profile_keys, function(key) {
  main_align_pathway_result(main_pathway_comparison_context$results[[key]], main_pathway_profile_ids)
}), main_pathway_profile_keys)

# Build matching NES/FDR matrices and verify every column's pathway definitions and labels.
main_pathway_reference_table <- main_pathway_aligned_results[[main_pathway_reference_key]]
main_pathway_profile_key <- main_pathway_reference_table[, c("Pathway_ID", "Pathway_name", "size"), drop = FALSE]
main_pathway_NES <- vapply(main_pathway_aligned_results, function(x) x$NES, numeric(length(main_pathway_profile_ids)))
main_pathway_FDR <- vapply(main_pathway_aligned_results, function(x) x$FDR, numeric(length(main_pathway_profile_ids)))
dimnames(main_pathway_NES) <- list(main_pathway_profile_ids, main_pathway_profile_keys)
dimnames(main_pathway_FDR) <- dimnames(main_pathway_NES)
main_pathway_definition_checks <- vapply(seq_along(main_pathway_profile_keys), function(i) {
  result <- main_pathway_aligned_results[[main_pathway_profile_keys[i]]]
  identical(result$Pathway_name, main_pathway_profile_key$Pathway_name) &&
    identical(result$size, main_pathway_profile_key$size) &&
    all(result$Study == main_pathway_profile_info$Study[i]) &&
    all(result$Comparison == main_pathway_profile_info$Comparison[i]) &&
    all(result$Role == main_pathway_profile_info$Role[i])
}, logical(1))
main_pathway_profile_checks <- c(
  Pathway_order = identical(rownames(main_pathway_NES), main_pathway_profile_ids) &&
    identical(rownames(main_pathway_FDR), main_pathway_profile_ids),
  Contrast_order = identical(colnames(main_pathway_NES), main_pathway_profile_keys) &&
    identical(colnames(main_pathway_FDR), main_pathway_profile_keys),
  Shared_definitions = all(main_pathway_definition_checks),
  Finite_NES = all(is.finite(main_pathway_NES)),
  Valid_FDR = all(is.finite(main_pathway_FDR) & main_pathway_FDR >= 0 & main_pathway_FDR <= 1),
  Nonconstant_profiles = all(apply(main_pathway_NES, 2, function(x) length(unique(x)) > 1L))
)
main_pathway_profile_dimensions <- rbind(NES = dim(main_pathway_NES), FDR = dim(main_pathway_FDR))
colnames(main_pathway_profile_dimensions) <- c("Pathways", "Contrasts")
print(main_pathway_profile_dimensions)
print(main_pathway_profile_checks)
stopifnot(all(main_pathway_profile_checks))

# Record the predeclared profile comparison and descriptive evidence-label rules.
main_pathway_profile_method <- list(
  measure = "Descriptive Spearman correlation of the complete retained-pathway NES profiles",
  reference = main_pathway_reference_key,
  pathways = length(main_pathway_profile_ids),
  pathway_selection = "Same complete retained pathway set in every contrast; no pathway-FDR filtering",
  input_proposal = main_pathway_comparison_context$method$planned_profile_comparison,
  ties = "Average ranks as used by stats::cor(method = 'spearman')",
  inference = "Descriptive point estimates only; no correlation P values, confidence intervals or independent-pathway assumption",
  dependency = "Reactome pathways can overlap or be nested; enrichment counts are not independent mechanisms",
  evidence_flags = "FDR < 0.05 separately within each contrast, unchanged from enrichment",
  one_contrast_flags = "Passing FDR only in one contrast does not establish biological exclusivity or a significant between-study difference",
  direction_flags = "Same/opposite nonzero NES signs; a zero NES is labelled separately",
  overlap_interpretation = "Both-enriched flags intersect separately adjusted results; no joint FDR or test of differential enrichment between studies",
  observed_order = "Sort rho decreasing; exact ties share a minimum observed rank; keep every terrestrial contrast separate",
  R_version = as.character(getRversion()),
  documentation = "https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html"
)
main_pathway_profile_method_file <- file.path(
  main_analysis_dir, "Main_analysis_OSD758_pathway_Spearman_method_2026-11-06.txt"
)
utils::capture.output(print(main_pathway_profile_method), file = main_pathway_profile_method_file)

# Compare every complete terrestrial NES profile with the same fixed flight profile.
main_pathway_profile_rho <- main_pathway_reference_spearman(main_pathway_NES, main_pathway_reference_key)
main_pathway_terrestrial_index <- match(names(main_pathway_profile_rho), main_pathway_profile_keys)
main_pathway_similarity_results <- data.frame(
  Study = main_pathway_profile_info$Study[main_pathway_terrestrial_index],
  Comparison = main_pathway_profile_info$Comparison[main_pathway_terrestrial_index],
  Common_pathways = nrow(main_pathway_NES), Spearman_rho = unname(main_pathway_profile_rho),
  stringsAsFactors = FALSE
)
main_pathway_similarity_results$Observed_rank <- rank(-main_pathway_similarity_results$Spearman_rho, ties.method = "min")
main_pathway_similarity_results <- main_pathway_similarity_results[
  order(-main_pathway_similarity_results$Spearman_rho, main_pathway_similarity_results$Study,
    main_pathway_similarity_results$Comparison),
  c("Observed_rank", "Study", "Comparison", "Common_pathways", "Spearman_rho"), drop = FALSE
]
rownames(main_pathway_similarity_results) <- NULL
print(main_pathway_similarity_results, digits = 5, row.names = FALSE)

# Match pathway identities, NES directions and the original FDR flags for biological review.
main_pathway_pair_tables <- lapply(main_pathway_terrestrial_index, function(i) {
  osd_nes <- main_pathway_NES[, main_pathway_reference_key]
  disease_nes <- main_pathway_NES[, i]
  same <- sign(osd_nes) == sign(disease_nes) & osd_nes != 0 & disease_nes != 0
  opposite <- sign(osd_nes) * sign(disease_nes) < 0
  osd_enriched <- main_pathway_FDR[, main_pathway_reference_key] < 0.05
  disease_enriched <- main_pathway_FDR[, i] < 0.05
  data.frame(
    Study = main_pathway_profile_info$Study[i], Comparison = main_pathway_profile_info$Comparison[i],
    main_pathway_profile_key,
    OSD_NES = unname(osd_nes), OSD_FDR = unname(main_pathway_FDR[, main_pathway_reference_key]),
    Terrestrial_NES = unname(disease_nes), Terrestrial_FDR = unname(main_pathway_FDR[, i]),
    NES_direction = ifelse(same, "Same", ifelse(opposite, "Opposite", "Zero_NES")),
    OSD_enriched = unname(osd_enriched), Terrestrial_enriched = unname(disease_enriched),
    Both_enriched = unname(osd_enriched & disease_enriched),
    stringsAsFactors = FALSE, row.names = NULL
  )
})
names(main_pathway_pair_tables) <- names(main_pathway_profile_rho)
main_pathway_pair_results <- do.call(rbind, main_pathway_pair_tables)
rownames(main_pathway_pair_results) <- NULL
main_pathway_overlap_summary <- do.call(rbind, lapply(main_pathway_pair_tables, function(x) {
  data.frame(
    Study = x$Study[1], Comparison = x$Comparison[1],
    Both_FDR005 = sum(x$Both_enriched),
    Same_direction_both = sum(x$Both_enriched & x$NES_direction == "Same"),
    Opposite_direction_both = sum(x$Both_enriched & x$NES_direction == "Opposite"),
    Zero_NES_both = sum(x$Both_enriched & x$NES_direction == "Zero_NES"),
    OSD_only_FDR005 = sum(x$OSD_enriched & !x$Terrestrial_enriched),
    Terrestrial_only_FDR005 = sum(!x$OSD_enriched & x$Terrestrial_enriched),
    stringsAsFactors = FALSE
  )
}))
rownames(main_pathway_overlap_summary) <- NULL
stopifnot(
  nrow(main_pathway_pair_results) == length(main_pathway_profile_ids) * length(main_pathway_profile_rho),
  all(main_pathway_overlap_summary$Both_FDR005 == main_pathway_overlap_summary$Same_direction_both +
    main_pathway_overlap_summary$Opposite_direction_both + main_pathway_overlap_summary$Zero_NES_both)
)
print(main_pathway_overlap_summary, row.names = FALSE)

# Save aligned matrices, all matched pathway records, score summaries and a plain checkpoint.
main_pathway_comparison_files <- file.path(main_analysis_dir, c(
  "Main_analysis_Reactome_common_NES.csv", "Main_analysis_Reactome_common_FDR.csv",
  "Main_analysis_OSD758_pathway_Spearman_common.csv", "Main_analysis_OSD758_pathway_pairs.csv",
  "Main_analysis_OSD758_pathway_overlap_summary.csv",
  "Main_analysis_pathway_comparison_checkpoint_2026-11-06.rds"
))
utils::write.csv(data.frame(main_pathway_profile_key, main_pathway_NES, check.names = FALSE),
  main_pathway_comparison_files[1], row.names = FALSE)
utils::write.csv(data.frame(main_pathway_profile_key, main_pathway_FDR, check.names = FALSE),
  main_pathway_comparison_files[2], row.names = FALSE)
utils::write.csv(main_pathway_similarity_results, main_pathway_comparison_files[3], row.names = FALSE)
utils::write.csv(main_pathway_pair_results, main_pathway_comparison_files[4], row.names = FALSE)
utils::write.csv(main_pathway_overlap_summary, main_pathway_comparison_files[5], row.names = FALSE)
main_pathway_comparison_checkpoint <- list(
  NES = main_pathway_NES, FDR = main_pathway_FDR, pathway_key = main_pathway_profile_key,
  dimensions = main_pathway_profile_dimensions, checks = main_pathway_profile_checks,
  similarity = main_pathway_similarity_results, rho_in_manifest_order = main_pathway_profile_rho,
  pathway_pairs = main_pathway_pair_results, overlap_summary = main_pathway_overlap_summary,
  contrasts = main_pathway_profile_info, method = main_pathway_profile_method,
  enrichment_checkpoint = main_pathway_comparison_input,
  enrichment_checkpoint_md5 = tools::md5sum(main_pathway_comparison_input),
  original_input_md5 = main_pathway_comparison_context$input_md5,
  enrichment_method = main_pathway_comparison_context$method,
  common_gene_ids = main_pathway_comparison_context$common_ids,
  resource_versions = main_pathway_comparison_context$resource_versions,
  resource_metadata = main_pathway_comparison_context$resource_metadata
)
saveRDS(main_pathway_comparison_checkpoint, main_pathway_comparison_files[6])
print(file.exists(c(main_pathway_comparison_files, main_pathway_profile_method_file)))
