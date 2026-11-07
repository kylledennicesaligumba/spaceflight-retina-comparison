# Audit the Entrez IDs already supplied by DR and the final gene-level AMD tables.
# Preserve source rows, existing statistics and AMD's fixed representative features.

# Check the supplied identifiers without substituting IDs or deleting rows.
main_audit_native_entrez <- function(source_ids, known_ids) {
  source_text <- as.character(source_ids)
  source_text[is.na(source_ids)] <- NA_character_
  native_ids <- source_text

  # Render numeric IDs as ordinary decimal text; keep character labels unchanged.
  if (is.numeric(source_ids)) {
    native_ids <- vapply(source_ids, function(x) {
      if (is.na(x)) NA_character_ else format(x, scientific = FALSE, trim = TRUE, digits = 22)
    }, character(1), USE.NAMES = FALSE)
  }
  has_id <- !is.na(source_text) & nzchar(source_text)
  valid_format <- has_id & !is.na(native_ids) & grepl("^[1-9][0-9]*$", native_ids)
  in_database <- valid_format & native_ids %in% known_ids
  audit <- data.frame(
    Source_row = seq_along(source_text), Source_ID = source_text,
    Entrez_ID = ifelse(valid_format, native_ids, NA_character_),
    In_database = in_database,
    Validation_status = ifelse(
      !has_id, "Missing_source_ID",
      ifelse(!valid_format, "Invalid_format",
        ifelse(in_database, "In_database", "Not_in_database"))
    ),
    row.names = NULL, stringsAsFactors = FALSE
  )

  # Flag repeated supplied targets, including format-valid IDs absent from this snapshot.
  audit$Target_repeated <- !is.na(audit$Entrez_ID) & (
    duplicated(audit$Entrez_ID) | duplicated(audit$Entrez_ID, fromLast = TRUE)
  )
  audit
}

# Confirm that the imported full results and reviewed annotation database are available.
stopifnot(
  exists("main_ranked_tables"), exists("main_mouse_db"), dir.exists(main_analysis_dir)
)

# Use the reviewed Entrez fields and retain the original source feature/record identifiers.
main_native_entrez_plan <- data.frame(
  Study = c("DR", "AMD", "AMD"),
  Comparison = c("Diabetic_vs_Normal", "Young_CfhKO_vs_WT", "Aged_CfhKO_vs_WT"),
  Source_column = c("ENTREZ_GENE_ID", "entrez_id", "entrez_id"),
  Record_column = c("ID", "feature_ID", "feature_ID"),
  stringsAsFactors = FALSE
)
main_native_entrez_plan$Table <- paste(
  main_native_entrez_plan$Study, main_native_entrez_plan$Comparison, sep = "__"
)
stopifnot(all(main_native_entrez_plan$Table %in% names(main_ranked_tables)))
main_db_entrez_keys <- AnnotationDbi::keys(main_mouse_db, keytype = "ENTREZID")
main_native_entrez_audits <- vector("list", nrow(main_native_entrez_plan))
names(main_native_entrez_audits) <- main_native_entrez_plan$Table
main_native_summary_rows <- vector("list", nrow(main_native_entrez_plan))
main_native_audit_files <- character(nrow(main_native_entrez_plan))

# Audit every original row using the same database snapshot as the other studies.
for (i in seq_len(nrow(main_native_entrez_plan))) {
  table_key <- main_native_entrez_plan$Table[i]
  source_column <- main_native_entrez_plan$Source_column[i]
  record_column <- main_native_entrez_plan$Record_column[i]
  result_table <- main_ranked_tables[[table_key]]
  stopifnot(all(c(source_column, record_column) %in% names(result_table)))
  audit <- main_audit_native_entrez(result_table[[source_column]], main_db_entrez_keys)
  audit$Original_record_ID <- as.character(result_table[[record_column]])
  audit$Source_column <- source_column
  stopifnot(nrow(audit) == nrow(result_table))
  main_native_entrez_audits[[table_key]] <- audit

  # Report database coverage and format issues separately from repeated-target flags.
  main_native_summary_rows[[i]] <- data.frame(
    Study = main_native_entrez_plan$Study[i],
    Comparison = main_native_entrez_plan$Comparison[i],
    Source_rows = nrow(audit),
    In_database = sum(audit$Validation_status == "In_database"),
    Not_in_database = sum(audit$Validation_status == "Not_in_database"),
    Invalid_format = sum(audit$Validation_status == "Invalid_format"),
    Missing_ID = sum(audit$Validation_status == "Missing_source_ID"),
    Rows_sharing_target = sum(audit$Target_repeated)
  )
  main_native_audit_files[i] <- file.path(
    main_analysis_dir, paste0(table_key, "_Entrez_ID_audit.csv")
  )
  write.csv(audit, main_native_audit_files[i], row.names = FALSE)
}

# Display and save the three-contrast audit without making eligibility decisions.
main_native_entrez_summary <- do.call(rbind, main_native_summary_rows)
stopifnot(with(main_native_entrez_summary,
  all(In_database + Not_in_database + Invalid_format + Missing_ID == Source_rows)
))
print(main_native_entrez_summary, row.names = FALSE)
main_native_summary_file <- file.path(main_analysis_dir, "Native_Entrez_ID_validation_summary.csv")
write.csv(main_native_entrez_summary, main_native_summary_file, row.names = FALSE)
print(file.exists(c(main_native_audit_files, main_native_summary_file)))
