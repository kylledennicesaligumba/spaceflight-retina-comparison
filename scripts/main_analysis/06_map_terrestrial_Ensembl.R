# Audit Ensembl-to-Entrez mappings for glaucoma, RP and light-damage responses.
# Apply the reference database and candidate-retention rule to every full table.

# Expand the existing unique-ID audit back to all original source rows.
# Retain missing IDs and repeated source labels as explicit audit records.
main_make_row_id_audit <- function(source_ids, candidates) {
  source_ids <- as.character(source_ids)
  has_id <- !is.na(source_ids) & nzchar(source_ids)
  audit <- data.frame(
    Source_row = seq_along(source_ids), Source_ID = source_ids,
    Candidate_Entrez = rep("", length(source_ids)),
    N_matches = rep(0L, length(source_ids)),
    Mapping_status = ifelse(has_id, "No_match", "Missing_source_ID"),
    Entrez_ID = rep(NA_character_, length(source_ids)),
    Target_repeated = rep(FALSE, length(source_ids)),
    stringsAsFactors = FALSE
  )
  if (any(has_id)) {
    unique_audit <- main_make_id_audit(unique(source_ids[has_id]), candidates)
    transfer_columns <- c("Candidate_Entrez", "N_matches", "Mapping_status", "Entrez_ID")
    audit[has_id, transfer_columns] <- unique_audit[
      match(source_ids[has_id], unique_audit$Source_ID), transfer_columns, drop = FALSE
    ]
  }

  # Recheck target repetition across all original rows after expansion.
  audit$Target_repeated <- !is.na(audit$Entrez_ID) & (
    duplicated(audit$Entrez_ID) | duplicated(audit$Entrez_ID, fromLast = TRUE)
  )
  audit
}

# Confirm that the reviewed source tables, database and reference helper are available.
stopifnot(
  exists("main_ranked_tables"), exists("main_mouse_db"),
  exists("main_make_id_audit"), dir.exists(main_analysis_dir)
)

# Use the exact source-identifier fields established in the user's previews.
main_ensembl_mapping_plan <- data.frame(
  Study = c(rep("Glaucoma", 3), "RP", "Light_damage"),
  Comparison = c("3D_vs_Naive", "2W_vs_Naive", "6W_vs_Naive", "P23H_vs_WT", "LD_vs_CTRL"),
  Source_column = c(rep("ensembl_id", 3), "gene_id", "Geneid"),
  Record_column = c(rep("original_label", 3), "gene_id", "Geneid"),
  stringsAsFactors = FALSE
)
main_ensembl_mapping_plan$Table <- paste(
  main_ensembl_mapping_plan$Study, main_ensembl_mapping_plan$Comparison, sep = "__"
)
stopifnot(all(main_ensembl_mapping_plan$Table %in% names(main_ranked_tables)))

# Read the database's supported Ensembl keys to handle unmatched inputs explicitly.
main_db_ensembl_keys <- AnnotationDbi::keys(main_mouse_db, keytype = "ENSEMBL")
main_ensembl_mapping_audits <- vector("list", nrow(main_ensembl_mapping_plan))
names(main_ensembl_mapping_audits) <- main_ensembl_mapping_plan$Table
main_ensembl_summary_rows <- vector("list", nrow(main_ensembl_mapping_plan))
main_ensembl_audit_files <- character(nrow(main_ensembl_mapping_plan))

# Retrieve every candidate while preserving each contrast and its original row order.
for (i in seq_len(nrow(main_ensembl_mapping_plan))) {
  table_key <- main_ensembl_mapping_plan$Table[i]
  source_column <- main_ensembl_mapping_plan$Source_column[i]
  record_column <- main_ensembl_mapping_plan$Record_column[i]
  result_table <- main_ranked_tables[[table_key]]
  stopifnot(all(c(source_column, record_column) %in% names(result_table)))
  source_ids <- as.character(result_table[[source_column]])
  usable_ids <- unique(source_ids[!is.na(source_ids) & nzchar(source_ids)])
  candidate_map <- stats::setNames(vector("list", length(usable_ids)), usable_ids)
  query_ids <- usable_ids[usable_ids %in% main_db_ensembl_keys]
  if (length(query_ids) > 0L) {
    query_matches <- AnnotationDbi::mapIds(
      main_mouse_db, keys = query_ids,
      column = "ENTREZID", keytype = "ENSEMBL", multiVals = "list"
    )
    candidate_map[query_ids] <- query_matches[query_ids]
  }
  audit <- main_make_row_id_audit(source_ids, candidate_map)
  audit$Original_record_ID <- as.character(result_table[[record_column]])
  audit$Source_column <- source_column
  stopifnot(nrow(audit) == nrow(result_table), identical(audit$Source_ID, source_ids))
  main_ensembl_mapping_audits[[table_key]] <- audit

  # Count the same mapping categories and separately report any missing source IDs.
  main_ensembl_summary_rows[[i]] <- data.frame(
    Study = main_ensembl_mapping_plan$Study[i],
    Comparison = main_ensembl_mapping_plan$Comparison[i],
    Source_rows = nrow(audit),
    One_match = sum(audit$Mapping_status == "One_match"),
    Multiple_matches = sum(audit$Mapping_status == "Multiple_matches"),
    No_match = sum(audit$Mapping_status == "No_match"),
    Missing_ID = sum(audit$Mapping_status == "Missing_source_ID"),
    Rows_sharing_target = sum(audit$Target_repeated)
  )

  # Save a full audit per contrast without modifying its source results.
  main_ensembl_audit_files[i] <- file.path(
    main_analysis_dir, paste0(table_key, "_Entrez_mapping_audit.csv")
  )
  write.csv(audit, main_ensembl_audit_files[i], row.names = FALSE)
}

# Display and save the five-contrast summary for review before eligibility decisions.
main_ensembl_mapping_summary <- do.call(rbind, main_ensembl_summary_rows)
stopifnot(with(main_ensembl_mapping_summary,
  all(One_match + Multiple_matches + No_match + Missing_ID == Source_rows)
))
print(main_ensembl_mapping_summary, row.names = FALSE)
main_ensembl_summary_file <- file.path(
  main_analysis_dir, "Terrestrial_Ensembl_Entrez_mapping_summary.csv"
)
write.csv(main_ensembl_mapping_summary, main_ensembl_summary_file, row.names = FALSE)
print(file.exists(c(main_ensembl_audit_files, main_ensembl_summary_file)))
