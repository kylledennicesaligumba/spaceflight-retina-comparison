# Assess provisional common-gene coverage using the completed identifier audits.
# Apply explicit mapping rules to full tested tables, independently of response statistics.

# Classify every audit row while retaining unresolved cases and their original identifiers.
main_make_eligibility_decisions <- function(audit, known_entrez_ids) {
  stopifnot(
    all(c("Source_row", "Source_ID", "Entrez_ID", "Target_repeated") %in% names(audit)),
    !anyNA(audit$Target_repeated)
  )
  target_ids <- as.character(audit$Entrez_ID)
  if ("Mapping_status" %in% names(audit)) {
    status <- as.character(audit$Mapping_status)
    stopifnot(
      all(status %in% c("One_match", "Multiple_matches", "No_match", "Missing_source_ID")),
      all((status == "One_match") == (audit$N_matches == 1L))
    )
    reason <- status
    reason[status == "One_match"] <- "Eligible"
    default_route <- "Ensembl"
  } else {
    stopifnot("Validation_status" %in% names(audit))
    status <- as.character(audit$Validation_status)
    stopifnot(all(status %in% c(
      "In_database", "Not_in_database", "Invalid_format", "Missing_source_ID"
    )))
    reason <- status
    reason[status == "In_database"] <- "Eligible"
    default_route <- "Supplied_Entrez"
  }

  # Require a represented target and exclude all repeated-target rows from this baseline.
  target_supported <- !is.na(target_ids) & nzchar(target_ids) & target_ids %in% known_entrez_ids
  reason[reason == "Eligible" & !target_supported] <- "Target_not_in_database"
  reason[reason == "Eligible" & audit$Target_repeated] <- "Repeated_target"
  eligible <- reason == "Eligible"
  stopifnot(anyDuplicated(target_ids[eligible]) == 0L)
  data.frame(
    Source_row = as.integer(audit$Source_row), Source_ID = as.character(audit$Source_ID),
    Original_record_ID = if ("Original_record_ID" %in% names(audit)) {
      as.character(audit$Original_record_ID)
    } else as.character(audit$Source_ID),
    Candidate_Entrez = if ("Candidate_Entrez" %in% names(audit)) {
      as.character(audit$Candidate_Entrez)
    } else target_ids,
    Entrez_ID = target_ids,
    Lookup_route = if ("Lookup_route" %in% names(audit)) {
      as.character(audit$Lookup_route)
    } else rep(default_route, nrow(audit)),
    Eligibility_status = reason, Eligible = eligible,
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Confirm that all reviewed audit groups and their source tables are available.
stopifnot(
  exists("main_osd_mapping_audit"), exists("main_ensembl_mapping_audits"),
  exists("main_native_entrez_audits"), exists("main_symbol_mapping_audits"),
  exists("main_ranked_tables"), exists("main_db_entrez_keys"),
  exists("main_mapping_versions"), exists("main_mapping_metadata"),
  dir.exists(main_analysis_dir)
)
main_contrast_keys <- paste(main_input_files$Study, main_input_files$Comparison, sep = "__")
main_mapping_registry <- c(
  list(OSD758__FLT_0G_vs_HGC = main_osd_mapping_audit),
  main_ensembl_mapping_audits, main_native_entrez_audits, main_symbol_mapping_audits
)
stopifnot(
  anyDuplicated(names(main_mapping_registry)) == 0L,
  setequal(names(main_mapping_registry), main_contrast_keys)
)
main_mapping_registry <- main_mapping_registry[main_contrast_keys]
main_mapping_registry[["OSD758__FLT_0G_vs_HGC"]]$Source_column <- "gene_id"

# Verify that every audit still refers to the same original source-table rows.
for (table_key in main_contrast_keys) {
  audit <- main_mapping_registry[[table_key]]
  result_table <- main_ranked_tables[[table_key]]
  source_column <- unique(as.character(audit$Source_column))
  stopifnot(
    length(source_column) == 1L, source_column %in% names(result_table),
    nrow(audit) == nrow(result_table),
    identical(as.integer(audit$Source_row), seq_len(nrow(result_table))),
    identical(as.character(audit$Source_ID), as.character(result_table[[source_column]]))
  )
}

# Record eligibility decisions without changing the source results or audit objects.
main_mapping_row_decisions <- lapply(
  main_mapping_registry, main_make_eligibility_decisions,
  known_entrez_ids = main_db_entrez_keys
)
main_eligible_row_maps <- lapply(main_mapping_row_decisions, function(x) {
  x[x$Eligible, c("Source_row", "Source_ID", "Original_record_ID", "Entrez_ID", "Lookup_route"),
    drop = FALSE]
})
main_eligible_gene_ids <- lapply(main_eligible_row_maps, function(x) x$Entrez_ID)

# Find genes represented on the same eligible basis in every labelled contrast.
# Sort identifiers independently of expression, statistical significance or rank score.
main_candidate_common_ids <- sort(Reduce(intersect, main_eligible_gene_ids), method = "radix")
main_candidate_common_genes <- data.frame(
  Entrez_ID = main_candidate_common_ids, stringsAsFactors = FALSE
)
main_reference_eligible_ids <- main_eligible_gene_ids[["OSD758__FLT_0G_vs_HGC"]]

# Report per-contrast coverage and available pairwise overlap with the reference.
# Pairwise counts describe coverage; no resemblance measure is calculated here.
main_gene_coverage <- data.frame(
  Study = main_input_files$Study, Comparison = main_input_files$Comparison,
  Tested_rows = vapply(main_ranked_tables[main_contrast_keys], nrow, integer(1)),
  Eligible_unique = vapply(main_eligible_row_maps, nrow, integer(1)),
  Flagged_rows = vapply(main_mapping_row_decisions, function(x) sum(!x$Eligible), numeric(1)),
  Shared_with_OSD = vapply(main_eligible_gene_ids, function(x) {
    length(intersect(x, main_reference_eligible_ids))
  }, integer(1)),
  row.names = NULL
)
stopifnot(with(main_gene_coverage, all(Eligible_unique + Flagged_rows == Tested_rows)))
print(main_gene_coverage, row.names = FALSE)
print(c(Candidate_common_genes_all_13 = length(main_candidate_common_ids)))

# Preserve the count of each eligibility reason for every study and contrast.
main_eligibility_reason_counts <- do.call(rbind, lapply(seq_along(main_contrast_keys), function(i) {
  counts <- table(main_mapping_row_decisions[[main_contrast_keys[i]]]$Eligibility_status)
  data.frame(
    Study = main_input_files$Study[i], Comparison = main_input_files$Comparison[i],
    Eligibility_status = names(counts), Rows = as.integer(counts), row.names = NULL
  )
}))
main_gene_coverage_policy <- list(
  status = "Provisional baseline for coverage review",
  target = paste0(
    "Mouse Entrez Gene IDs represented in org.Mm.eg.db ",
    main_mapping_versions$Version[main_mapping_versions$Component == "org.Mm.eg.db"]
  ),
  converted_ids = "One candidate per row; exclude every repeated-target row",
  native_ids = "Supplied ID present in the same database; no repeated target",
  source_rows = "All tested genes considered independently of FDR, effect size and rank score",
  common_set = "Intersection across the reference and all 12 labelled terrestrial contrasts",
  unresolved_cases = "Preserved in original results, mapping audits and row decisions"
)

# Record source-file fingerprints so the saved row maps remain traceable after a restart.
main_coverage_input_md5 <- tools::md5sum(main_input_files$File)
stopifnot(!anyNA(main_coverage_input_md5))
main_coverage_checkpoint <- list(
  policy = main_gene_coverage_policy, inputs = main_input_files,
  input_md5 = main_coverage_input_md5,
  annotation_versions = main_mapping_versions, annotation_metadata = main_mapping_metadata,
  row_decisions = main_mapping_row_decisions, eligible_row_maps = main_eligible_row_maps,
  candidate_common_ids = main_candidate_common_ids,
  coverage = main_gene_coverage, reason_counts = main_eligibility_reason_counts
)

# Save coverage, the candidate common IDs, decision counts and a reloadable mapping checkpoint.
# The checkpoint stores plain tables and identifiers, without the live database connection.
main_gene_coverage_files <- file.path(
  main_analysis_dir,
  c(
    "Main_analysis_gene_coverage.csv",
    "Main_analysis_candidate_common_genes.csv",
    "Main_analysis_eligibility_reason_counts.csv",
    "Main_analysis_gene_coverage_checkpoint_2026-11-06.rds"
  )
)
utils::write.csv(main_gene_coverage, main_gene_coverage_files[1], row.names = FALSE)
utils::write.csv(main_candidate_common_genes, main_gene_coverage_files[2], row.names = FALSE)
utils::write.csv(main_eligibility_reason_counts, main_gene_coverage_files[3], row.names = FALSE)
saveRDS(main_coverage_checkpoint, main_gene_coverage_files[4])
print(file.exists(main_gene_coverage_files))
