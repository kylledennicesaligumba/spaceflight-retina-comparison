# Align the existing study responses on the same shared genes before comparing them.
# Preserve original logFC, signed rank scores and within-study adjusted P values.

# Retrieve common-gene values by verified Entrez-to-source-row links.
main_extract_common_response <- function(result_table, source_map, common_ids, source_column) {
  q_columns <- intersect(c("FDR", "adj.P.Val"), names(result_table))
  stopifnot(
    all(c(source_column, "logFC", "rank_score") %in% names(result_table)),
    length(q_columns) == 1L,
    anyDuplicated(source_map$Entrez_ID) == 0L,
    anyDuplicated(common_ids) == 0L
  )
  map_index <- match(common_ids, as.character(source_map$Entrez_ID))
  stopifnot(!anyNA(map_index))
  aligned_map <- source_map[map_index, , drop = FALSE]
  source_rows <- aligned_map$Source_row
  stopifnot(
    !anyNA(source_rows), all(source_rows == floor(source_rows)),
    all(source_rows >= 1L & source_rows <= nrow(result_table)),
    identical(as.character(aligned_map$Entrez_ID), common_ids)
  )
  source_id_alignment <- identical(
    as.character(result_table[[source_column]][source_rows]), as.character(aligned_map$Source_ID)
  )
  stopifnot(source_id_alignment)
  values <- result_table[source_rows, c("rank_score", "logFC", q_columns), drop = FALSE]
  stopifnot(
    all(vapply(values, is.numeric, logical(1))),
    all(vapply(values, function(x) all(is.finite(x)), logical(1)))
  )
  list(
    rank_score = values$rank_score, logFC = values$logFC, FDR = values[[q_columns]],
    source_map = aligned_map, source_id_alignment = source_id_alignment,
    q_source_column = q_columns
  )
}

# Confirm that the reviewed common-gene set, source tables and saved mapping context exist.
stopifnot(
  exists("main_candidate_common_ids"), exists("main_eligible_row_maps"),
  exists("main_ranked_tables"), exists("main_coverage_checkpoint"),
  exists("main_gene_coverage_policy"), dir.exists(main_analysis_dir)
)
main_common_ids <- as.character(main_candidate_common_ids)
main_response_keys <- paste(main_input_files$Study, main_input_files$Comparison, sep = "__")
stopifnot(
  length(main_common_ids) > 1L, anyDuplicated(main_common_ids) == 0L,
  identical(main_common_ids, main_coverage_checkpoint$candidate_common_ids),
  all(main_response_keys %in% names(main_eligible_row_maps)),
  all(main_response_keys %in% names(main_ranked_tables))
)

# Check that the saved source files still match the mapping checkpoint's fingerprints.
main_common_input_md5 <- tools::md5sum(main_input_files$File)
stopifnot(identical(main_common_input_md5, main_coverage_checkpoint$input_md5))
main_response_source_columns <- c(
  OSD758 = "gene_id", Glaucoma = "ensembl_id", RP = "gene_id", LHON = "gene_name",
  DR = "ENTREZ_GENE_ID", AMD = "entrez_id", Autoimmune_uveitis = "AccID", Light_damage = "Geneid"
)

# Give every matrix the same Entrez row order and labelled contrast columns.
main_common_rank_scores <- matrix(
  NA_real_, nrow = length(main_common_ids), ncol = length(main_response_keys),
  dimnames = list(main_common_ids, main_response_keys)
)
main_common_logFC <- main_common_rank_scores
main_common_FDR <- main_common_rank_scores
main_common_source_alignment <- stats::setNames(rep(FALSE, length(main_response_keys)), main_response_keys)
main_common_q_columns <- stats::setNames(character(length(main_response_keys)), main_response_keys)
main_common_source_keys <- vector("list", length(main_response_keys))

# Fill each column from the original result rows without refitting or rescaling responses.
for (i in seq_along(main_response_keys)) {
  table_key <- main_response_keys[i]
  source_column <- main_response_source_columns[[main_input_files$Study[i]]]
  extracted <- main_extract_common_response(
    main_ranked_tables[[table_key]], main_eligible_row_maps[[table_key]],
    main_common_ids, source_column
  )
  main_common_rank_scores[, i] <- extracted$rank_score
  main_common_logFC[, i] <- extracted$logFC
  main_common_FDR[, i] <- extracted$FDR
  main_common_source_alignment[i] <- extracted$source_id_alignment
  main_common_q_columns[i] <- extracted$q_source_column
  main_common_source_keys[[i]] <- data.frame(
    Study = main_input_files$Study[i], Comparison = main_input_files$Comparison[i],
    Entrez_ID = main_common_ids, Source_column = source_column,
    Source_ID = extracted$source_map$Source_ID,
    Original_record_ID = extracted$source_map$Original_record_ID,
    Source_row = extracted$source_map$Source_row,
    Lookup_route = extracted$source_map$Lookup_route,
    row.names = NULL, stringsAsFactors = FALSE
  )
}
main_common_source_row_key <- do.call(rbind, main_common_source_keys)
rownames(main_common_source_row_key) <- NULL

# Verify shared ordering, usable values and preserved response direction.
# Zero scores from the existing numerical guard are allowed to accompany nonzero logFC.
main_common_response_matrices <- list(
  Rank_scores = main_common_rank_scores, logFC = main_common_logFC, FDR = main_common_FDR
)
main_common_response_checks <- c(
  Gene_order = all(vapply(main_common_response_matrices, function(x) {
    identical(rownames(x), main_common_ids)
  }, logical(1))),
  Contrast_order = all(vapply(main_common_response_matrices, function(x) {
    identical(colnames(x), main_response_keys)
  }, logical(1))),
  Source_ID_alignment = all(main_common_source_alignment),
  Finite_scores = all(is.finite(main_common_rank_scores)),
  Finite_logFC = all(is.finite(main_common_logFC)),
  Valid_FDR = all(is.finite(main_common_FDR) & main_common_FDR >= 0 & main_common_FDR <= 1),
  Score_direction = all(main_common_rank_scores == 0 | sign(main_common_rank_scores) == sign(main_common_logFC)),
  Nonconstant_scores = all(apply(main_common_rank_scores, 2, function(x) length(unique(x)) > 1L))
)
main_common_response_dimensions <- do.call(rbind, lapply(main_common_response_matrices, dim))
colnames(main_common_response_dimensions) <- c("Genes", "Contrasts")
print(main_common_response_dimensions)
print(main_common_response_checks)
stopifnot(all(main_common_response_checks))

# Record the proposed first gene-level measure before any similarity result is calculated.
main_gene_comparison_proposal <- list(
  measure = "Descriptive Spearman correlation of the existing signed gene scores",
  genes = "The same accepted all-13 common set, with no FDR filtering",
  reference = "OSD758__FLT_0G_vs_HGC",
  contrasts = "Keep all 12 terrestrial stages/genotypes/ages separate; no post-result best-stage selection",
  coverage_sensitivity = "Revisit available reference-pairwise coverage using the same measure",
  stage = "Aligned inputs only; similarity, pathway and uncertainty calculations remain later"
)

# Save explicit gene-ID columns in the three matrices and preserve every source-row link.
main_common_response_files <- file.path(
  main_analysis_dir,
  c(
    "Main_analysis_common_rank_scores.csv", "Main_analysis_common_logFC.csv",
    "Main_analysis_common_FDR.csv", "Main_analysis_common_response_row_key.csv",
    "Main_analysis_common_responses_checkpoint_2026-11-06.rds"
  )
)
for (i in seq_along(main_common_response_matrices)) {
  output_table <- data.frame(
    Entrez_ID = main_common_ids, main_common_response_matrices[[i]],
    check.names = FALSE, row.names = NULL
  )
  utils::write.csv(output_table, main_common_response_files[i], row.names = FALSE)
}
utils::write.csv(main_common_source_row_key, main_common_response_files[4], row.names = FALSE)

# Preserve plain response matrices and provenance for restoration without database connections.
main_common_responses_checkpoint <- list(
  common_ids = main_common_ids, inputs = main_input_files, input_md5 = main_common_input_md5,
  dimensions = main_common_response_dimensions, checks = main_common_response_checks,
  rank_scores = main_common_rank_scores, logFC = main_common_logFC, FDR = main_common_FDR,
  source_row_key = main_common_source_row_key, adjusted_P_source_columns = main_common_q_columns,
  original_tested_counts = vapply(main_ranked_tables[main_response_keys], nrow, integer(1)),
  mapping_policy = main_gene_coverage_policy,
  annotation_versions = main_coverage_checkpoint$annotation_versions,
  annotation_metadata = main_coverage_checkpoint$annotation_metadata,
  comparison_proposal = main_gene_comparison_proposal
)
saveRDS(main_common_responses_checkpoint, main_common_response_files[5])
print(file.exists(main_common_response_files))
