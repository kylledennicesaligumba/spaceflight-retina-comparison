# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Test the same retained Reactome pathways in every separately labelled contrast.
# Use the original signed scores for all common genes, with no gene-FDR selection.
# https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf

# Sort unchanged scores reproducibly, retaining all genes and their Entrez identifiers.
main_prepare_fgsea_ranks <- function(scores, gene_ids) {
  stopifnot(
    is.numeric(scores), is.character(gene_ids), length(scores) == length(gene_ids),
    length(scores) > 1L, all(is.finite(scores)), length(unique(scores)) > 1L,
    !anyNA(gene_ids), all(nzchar(gene_ids)), anyDuplicated(gene_ids) == 0L
  )
  ranks <- stats::setNames(as.numeric(scores), gene_ids)
  ranks[order(-ranks, names(ranks), method = "radix")]
}

# Restore the fixed pathway memberships and their exact common-response inputs.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
stopifnot(requireNamespace("fgsea", quietly = TRUE))
main_fgsea_match_file <- file.path(
  main_analysis_dir, "Main_analysis_Reactome_common_pathways_checkpoint_2026-11-06.rds"
)
main_fgsea_response_file <- file.path(
  main_analysis_dir, "Main_analysis_common_responses_checkpoint_2026-11-06.rds"
)
stopifnot(file.exists(main_fgsea_match_file), file.exists(main_fgsea_response_file))
main_fgsea_match_context <- readRDS(main_fgsea_match_file)
main_fgsea_response_context <- readRDS(main_fgsea_response_file)
main_fgsea_scores <- main_fgsea_response_context$rank_scores
main_fgsea_pathways <- main_fgsea_match_context$pathways
main_fgsea_common_ids <- main_fgsea_match_context$common_ids
main_fgsea_contrasts <- main_fgsea_response_context$inputs
main_fgsea_keys <- paste(main_fgsea_contrasts$Study, main_fgsea_contrasts$Comparison, sep = "__")
main_fgsea_reference <- "OSD758__FLT_0G_vs_HGC"
main_fgsea_catalog <- main_fgsea_match_context$coverage[
  match(names(main_fgsea_pathways), main_fgsea_match_context$coverage$Pathway_ID), , drop = FALSE
]
stopifnot(
  length(main_fgsea_match_context$checks) > 0L, all(main_fgsea_match_context$checks),
  length(main_fgsea_response_context$checks) > 0L, all(main_fgsea_response_context$checks),
  identical(unname(tools::md5sum(main_fgsea_response_file)),
    unname(main_fgsea_match_context$response_checkpoint_md5)),
  identical(main_fgsea_common_ids, main_fgsea_response_context$common_ids),
  identical(rownames(main_fgsea_scores), main_fgsea_common_ids),
  identical(colnames(main_fgsea_scores), main_fgsea_keys),
  anyDuplicated(main_fgsea_keys) == 0L, main_fgsea_reference %in% main_fgsea_keys,
  sum(main_fgsea_contrasts$Role == "Reference") == 1L,
  main_fgsea_contrasts$Role[match(main_fgsea_reference, main_fgsea_keys)] == "Reference",
  all(main_fgsea_contrasts$Role[main_fgsea_keys != main_fgsea_reference] == "Terrestrial"),
  identical(main_fgsea_catalog$Pathway_ID, names(main_fgsea_pathways)),
  all(main_fgsea_catalog$Eligible_for_enrichment),
  all(is.finite(main_fgsea_scores)),
  identical(as.character(utils::packageVersion("fgsea")),
    main_fgsea_match_context$resource_versions$Version[
      main_fgsea_match_context$resource_versions$Component == "fgsea"])
)

# Declare the complete enrichment rule and later comparison proposal before results.
main_fgsea_seeds <- stats::setNames(20261106L + seq_along(main_fgsea_keys), main_fgsea_keys)
main_fgsea_method <- list(
  algorithm = "fgsea::fgseaMultilevel",
  gene_universe = "All accepted all-13 common Entrez IDs; no FDR filtering or removal of unannotated genes",
  gene_count = length(main_fgsea_common_ids), pathway_count = length(main_fgsea_pathways),
  inputs = "Original signed within-study scores: sign(logFC)*sqrt(F) for edgeR and moderated t for limma",
  ties = "Presort equal scores by Entrez ID; retain exact values without jitter; record fgsea warnings",
  minSize = main_fgsea_match_context$size_rule$minSize,
  maxSize = main_fgsea_match_context$size_rule$maxSize,
  scoreType = "std", gseaParam = 1, sampleSize = 101L, nPermSimple = 10000L,
  eps = 0, nproc = 1L, seeds = main_fgsea_seeds,
  RNG = c("Mersenne-Twister", "Inversion", "Rejection"),
  multiple_testing = "BH across the full retained pathway family within each contrast, including unestimated tests in the family size",
  FDR_cutoff = 0.05,
  interpretation = "Positive/negative NES denotes enrichment toward the higher/lower-expression end; not direct proof of pathway activation/inhibition",
  uncertainty = "Gene-set randomization conditional on stored gene scores; no biological-sample resampling or disease-ranking confidence intervals",
  contrast_handling = "Keep every reference/terrestrial stage, genotype and age separate",
  planned_profile_comparison = "Descriptive Spearman correlation of all retained-pathway NES values against the fixed OSD-758 profile, without pathway-FDR filtering; resolve missing estimates before that comparison",
  versions = c(R = as.character(getRversion()),
    fgsea = as.character(utils::packageVersion("fgsea")),
    BiocParallel = as.character(utils::packageVersion("BiocParallel"))),
  script_md5 = unname(tools::md5sum(file.path(main_analysis_dir, "17_run_Reactome_fgsea.R")))
)
main_fgsea_input_md5 <- tools::md5sum(c(main_fgsea_response_file, main_fgsea_match_file))
main_fgsea_progress_file <- file.path(
  main_analysis_dir, "Main_analysis_Reactome_fgsea_progress_checkpoint_2026-11-06.rds"
)
main_fgsea_method_file <- file.path(main_analysis_dir, "Main_analysis_Reactome_fgsea_method_2026-11-06.txt")
utils::capture.output(print(main_fgsea_method), file = main_fgsea_method_file)

# Reuse completed contrasts only when the saved inputs, code and settings are identical.
main_fgsea_results <- list()
main_fgsea_raw_results <- list()
main_fgsea_warning_lists <- list()
if (file.exists(main_fgsea_progress_file)) {
  main_fgsea_previous <- readRDS(main_fgsea_progress_file)
  if (identical(main_fgsea_previous$input_md5, main_fgsea_input_md5) &&
      identical(main_fgsea_previous$method, main_fgsea_method)) {
    main_fgsea_results <- main_fgsea_previous$results
    main_fgsea_raw_results <- main_fgsea_previous$raw_results
    main_fgsea_warning_lists <- main_fgsea_previous$warnings
  }
}

# Run enrichment separately and save progress after each completed contrast.
for (i in seq_along(main_fgsea_keys)) {
  contrast_key <- main_fgsea_keys[i]
  if (contrast_key %in% names(main_fgsea_results)) {
    message("Restored completed enrichment: ", contrast_key)
    next
  }
  message("Enrichment ", i, "/", length(main_fgsea_keys), ": ", contrast_key)
  contrast_ranks <- main_prepare_fgsea_ranks(main_fgsea_scores[, i], main_fgsea_common_ids)
  contrast_warnings <- character()
  set.seed(main_fgsea_seeds[[contrast_key]], kind = "Mersenne-Twister",
    normal.kind = "Inversion", sample.kind = "Rejection")
  contrast_raw <- withCallingHandlers(
    as.data.frame(fgsea::fgseaMultilevel(
      pathways = main_fgsea_pathways, stats = contrast_ranks,
      minSize = main_fgsea_method$minSize, maxSize = main_fgsea_method$maxSize,
      scoreType = main_fgsea_method$scoreType, gseaParam = main_fgsea_method$gseaParam,
      sampleSize = main_fgsea_method$sampleSize, nPermSimple = main_fgsea_method$nPermSimple,
      eps = main_fgsea_method$eps, nproc = main_fgsea_method$nproc
    )),
    warning = function(w) { contrast_warnings <<- c(contrast_warnings, conditionMessage(w)) }
  )
  stopifnot(
    all(c("pathway", "size", "ES", "NES", "pval", "padj", "log2err", "leadingEdge") %in% names(contrast_raw)),
    nrow(contrast_raw) == length(main_fgsea_pathways),
    anyDuplicated(contrast_raw$pathway) == 0L,
    setequal(contrast_raw$pathway, names(main_fgsea_pathways))
  )
  contrast_raw <- contrast_raw[match(names(main_fgsea_pathways), contrast_raw$pathway), , drop = FALSE]
  stopifnot(
    identical(as.integer(contrast_raw$size), unname(lengths(main_fgsea_pathways))),
    all(is.na(contrast_raw$pval) | (is.finite(contrast_raw$pval) & contrast_raw$pval >= 0 & contrast_raw$pval <= 1))
  )
  contrast_result <- data.frame(
    Study = main_fgsea_contrasts$Study[i], Comparison = main_fgsea_contrasts$Comparison[i],
    Role = main_fgsea_contrasts$Role[i], Pathway_ID = contrast_raw$pathway,
    Pathway_name = main_fgsea_catalog$Pathway_name, size = contrast_raw$size,
    ES = contrast_raw$ES, NES = contrast_raw$NES, pval = contrast_raw$pval,
    FDR = stats::p.adjust(contrast_raw$pval, method = "BH", n = length(main_fgsea_pathways)),
    log2err = contrast_raw$log2err, stringsAsFactors = FALSE
  )
  contrast_result$Estimation_status <- ifelse(
    is.finite(contrast_result$NES) & is.finite(contrast_result$pval) & is.finite(contrast_result$FDR),
    "Estimated", "Unestimated"
  )
  contrast_result$leadingEdge <- contrast_raw$leadingEdge
  main_fgsea_results[[contrast_key]] <- contrast_result
  main_fgsea_raw_results[[contrast_key]] <- contrast_raw
  main_fgsea_warning_lists[[contrast_key]] <- unique(contrast_warnings)
  saveRDS(list(results = main_fgsea_results, raw_results = main_fgsea_raw_results,
    warnings = main_fgsea_warning_lists, method = main_fgsea_method,
    input_md5 = main_fgsea_input_md5), main_fgsea_progress_file)
}

# Summarize estimated and significant pathways without treating missing estimates as zero.
main_fgsea_summary <- do.call(rbind, lapply(seq_along(main_fgsea_keys), function(i) {
  result <- main_fgsea_results[[main_fgsea_keys[i]]]
  significant <- result$Estimation_status == "Estimated" & !is.na(result$FDR) & result$FDR < 0.05
  data.frame(
    Study = main_fgsea_contrasts$Study[i], Comparison = main_fgsea_contrasts$Comparison[i],
    Pathways = nrow(result), Unestimated = sum(result$Estimation_status == "Unestimated"),
    Enriched_FDR005 = sum(significant),
    Positive_NES_FDR005 = sum(significant & result$NES > 0, na.rm = TRUE),
    Negative_NES_FDR005 = sum(significant & result$NES < 0, na.rm = TRUE),
    Tied_score_genes = sum(duplicated(main_fgsea_scores[, i]) |
      duplicated(main_fgsea_scores[, i], fromLast = TRUE)),
    Warnings = length(main_fgsea_warning_lists[[main_fgsea_keys[i]]])
  )
}))
main_fgsea_warning_table <- do.call(rbind, lapply(seq_along(main_fgsea_keys), function(i) {
  warnings <- main_fgsea_warning_lists[[main_fgsea_keys[i]]]
  data.frame(Study = rep(main_fgsea_contrasts$Study[i], length(warnings)),
    Comparison = rep(main_fgsea_contrasts$Comparison[i], length(warnings)), Warning = warnings)
}))
main_fgsea_all_results <- do.call(rbind, main_fgsea_results[main_fgsea_keys])
rownames(main_fgsea_all_results) <- NULL
print(main_fgsea_summary, row.names = FALSE)

# Show the lowest-P reference pathways for the next manual interpretation.
main_fgsea_reference_results <- main_fgsea_results[[main_fgsea_reference]]
main_fgsea_reference_preview <- main_fgsea_reference_results[
  order(main_fgsea_reference_results$pval, na.last = TRUE),
  c("Pathway_name", "size", "NES", "pval", "FDR", "Estimation_status"), drop = FALSE
]
print(utils::head(main_fgsea_reference_preview, 10L), digits = 4, row.names = FALSE)

# Save all pathway rows, leading-edge IDs, warnings, methods and restorable results.
main_fgsea_output_files <- file.path(main_analysis_dir, c(
  "Main_analysis_Reactome_fgsea_all_contrasts.csv",
  "Main_analysis_Reactome_fgsea_summary.csv",
  "Main_analysis_Reactome_fgsea_warnings.csv",
  "Main_analysis_Reactome_fgsea_checkpoint_2026-11-06.rds",
  "Main_analysis_Reactome_fgsea_sessionInfo_2026-11-06.txt"
))
main_fgsea_flat_results <- main_fgsea_all_results
main_fgsea_flat_results$leadingEdge <- vapply(main_fgsea_flat_results$leadingEdge,
  function(ids) paste(ids, collapse = ";"), character(1))
utils::write.csv(main_fgsea_flat_results, main_fgsea_output_files[1], row.names = FALSE)
utils::write.csv(main_fgsea_summary, main_fgsea_output_files[2], row.names = FALSE)
utils::write.csv(main_fgsea_warning_table, main_fgsea_output_files[3], row.names = FALSE)
main_fgsea_checkpoint <- list(
  results = main_fgsea_results, raw_results = main_fgsea_raw_results,
  summary = main_fgsea_summary, warnings = main_fgsea_warning_table,
  method = main_fgsea_method, input_md5 = main_fgsea_input_md5,
  common_ids = main_fgsea_common_ids, pathway_ids = names(main_fgsea_pathways),
  contrasts = main_fgsea_contrasts, response_checkpoint = main_fgsea_response_file,
  pathway_checkpoint = main_fgsea_match_file,
  resource_versions = main_fgsea_match_context$resource_versions,
  resource_metadata = main_fgsea_match_context$resource_metadata
)
saveRDS(main_fgsea_checkpoint, main_fgsea_output_files[4])
utils::capture.output(utils::sessionInfo(), file = main_fgsea_output_files[5])
print(file.exists(c(main_fgsea_output_files, main_fgsea_method_file, main_fgsea_progress_file)))
