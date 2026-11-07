# Audit candidate Ensembl-to-Entrez mappings for the fixed OSD-758 reference.
# Keep every tested row and all candidate matches before deciding eligibility.
# https://bioconductor.org/packages/release/bioc/manuals/AnnotationDbi/man/AnnotationDbi.pdf

# Build an audit table without selecting a first match or combining source rows.
main_make_id_audit <- function(source_ids, candidates) {
  stopifnot(
    is.list(candidates), !is.null(names(candidates)),
    !anyNA(source_ids), all(nzchar(source_ids)),
    anyDuplicated(source_ids) == 0L,
    anyDuplicated(names(candidates)) == 0L
  )

  # Preserve source order and count distinct, nonmissing candidate IDs for each row.
  candidate_lists <- lapply(candidates[source_ids], function(x) {
    unique(as.character(x[!is.na(x) & nzchar(x)]))
  })
  n_matches <- lengths(candidate_lists)
  candidate_text <- vapply(candidate_lists, paste, character(1), collapse = ";")
  audit <- data.frame(
    Source_row = seq_along(source_ids),
    Source_ID = source_ids,
    Candidate_Entrez = candidate_text,
    N_matches = n_matches,
    Mapping_status = ifelse(
      n_matches == 0L, "No_match",
      ifelse(n_matches == 1L, "One_match", "Multiple_matches")
    ),
    Entrez_ID = ifelse(n_matches == 1L, candidate_text, NA_character_),
    row.names = NULL, stringsAsFactors = FALSE
  )

  # Flag repeated targets among one-match rows while retaining every source row.
  audit$Target_repeated <- !is.na(audit$Entrez_ID) & (
    duplicated(audit$Entrez_ID) | duplicated(audit$Entrez_ID, fromLast = TRUE)
  )
  audit
}

# Confirm that the imported reference and the reviewed mouse database are available.
stopifnot(
  exists("main_ranked_tables"), exists("main_mouse_db"),
  dir.exists(main_analysis_dir),
  "OSD758__FLT_0G_vs_HGC" %in% names(main_ranked_tables)
)
main_osd_source_ids <- as.character(
  main_ranked_tables[["OSD758__FLT_0G_vs_HGC"]]$gene_id
)
stopifnot(
  length(main_osd_source_ids) == nrow(main_ranked_tables[["OSD758__FLT_0G_vs_HGC"]]),
  !anyNA(main_osd_source_ids), all(nzchar(main_osd_source_ids)),
  anyDuplicated(main_osd_source_ids) == 0L
)

# Retrieve every candidate Entrez match for each original Ensembl identifier.
# Keeping list results preserves any one-to-many mappings for review.
main_osd_entrez_candidates <- AnnotationDbi::mapIds(
  main_mouse_db,
  keys = main_osd_source_ids,
  column = "ENTREZID", keytype = "ENSEMBL",
  multiVals = "list"
)
main_osd_mapping_audit <- main_make_id_audit(
  main_osd_source_ids, main_osd_entrez_candidates
)

# Summarize source-to-target matches and repeated targets before filtering any genes.
# Rows_sharing_target overlaps One_match; it is a separate duplication flag.
main_osd_mapping_summary <- data.frame(
  Study = "OSD758",
  Source_rows = nrow(main_osd_mapping_audit),
  One_match = sum(main_osd_mapping_audit$N_matches == 1L),
  Multiple_matches = sum(main_osd_mapping_audit$N_matches > 1L),
  No_match = sum(main_osd_mapping_audit$N_matches == 0L),
  Rows_sharing_target = sum(main_osd_mapping_audit$Target_repeated)
)
print(main_osd_mapping_summary, row.names = FALSE)
print(head(main_osd_mapping_audit, n = 5), row.names = FALSE)

# Save the complete row-level mapping audit and summary beside the input records.
main_osd_mapping_files <- file.path(
  main_analysis_dir,
  c("OSD758_Entrez_mapping_audit.csv", "OSD758_Entrez_mapping_summary.csv")
)
write.csv(main_osd_mapping_audit, main_osd_mapping_files[1], row.names = FALSE)
write.csv(main_osd_mapping_summary, main_osd_mapping_files[2], row.names = FALSE)
print(file.exists(main_osd_mapping_files))
