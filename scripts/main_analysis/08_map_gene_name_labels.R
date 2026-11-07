# Audit the original LHON and uveitis labels using the same mouse annotation snapshot.
# Prefer exact official symbols; use aliases only when no official-symbol match exists.
# https://bioconductor.org/packages/release/data/annotation/manuals/org.Mm.eg.db/man/org.Mm.eg.db.pdf

# Retrieve all candidates for supported keys, retaining unmatched labels as empty lists.
main_lookup_candidate_list <- function(source_ids, keytype, known_keys, database) {
  source_ids <- unique(as.character(source_ids))
  source_ids <- source_ids[!is.na(source_ids) & nzchar(source_ids)]
  candidates <- stats::setNames(vector("list", length(source_ids)), source_ids)
  query_ids <- source_ids[source_ids %in% known_keys]
  if (length(query_ids) > 0L) {
    query_matches <- AnnotationDbi::mapIds(
      database, keys = query_ids, column = "ENTREZID",
      keytype = keytype, multiVals = "list"
    )
    candidates[query_ids] <- query_matches[query_ids]
  }
  lapply(candidates, function(x) unique(as.character(x[!is.na(x) & nzchar(x)])))
}

# Keep every candidate from the preferred lookup route and record that route explicitly.
main_choose_label_candidates <- function(source_ids, ensembl_ids, ensembl_map, symbol_map, alias_map) {
  candidates <- stats::setNames(vector("list", length(source_ids)), source_ids)
  routes <- stats::setNames(rep("Unmatched", length(source_ids)), source_ids)
  ensembl_ids <- intersect(source_ids, ensembl_ids)
  if (length(ensembl_ids) > 0L) {
    candidates[ensembl_ids] <- ensembl_map[ensembl_ids]
    routes[ensembl_ids] <- "Ensembl"
  }
  symbol_hits <- intersect(
    setdiff(source_ids, ensembl_ids), names(symbol_map)[lengths(symbol_map) > 0L]
  )
  if (length(symbol_hits) > 0L) {
    candidates[symbol_hits] <- symbol_map[symbol_hits]
    routes[symbol_hits] <- "Official_symbol"
  }
  alias_hits <- intersect(
    setdiff(source_ids, c(ensembl_ids, symbol_hits)), names(alias_map)[lengths(alias_map) > 0L]
  )
  if (length(alias_hits) > 0L) {
    candidates[alias_hits] <- alias_map[alias_hits]
    routes[alias_hits] <- "Alias"
  }
  list(Candidates = candidates, Route = routes)
}

# Confirm that the full result tables, reviewed database and row-audit helpers are available.
stopifnot(
  exists("main_ranked_tables"), exists("main_mouse_db"),
  exists("main_make_id_audit"), exists("main_make_row_id_audit"),
  dir.exists(main_analysis_dir)
)

# Use the original labels established in the user's previews for each contrast.
main_symbol_mapping_plan <- data.frame(
  Study = c(rep("LHON", 3), "Autoimmune_uveitis"),
  Comparison = c("YARS2_vs_WT", "COI_vs_WT", "YARS2_COI_vs_WT", "EAU_vs_CON"),
  Source_column = c(rep("gene_name", 3), "AccID"),
  stringsAsFactors = FALSE
)
main_symbol_mapping_plan$Table <- paste(
  main_symbol_mapping_plan$Study, main_symbol_mapping_plan$Comparison, sep = "__"
)
stopifnot(all(main_symbol_mapping_plan$Table %in% names(main_ranked_tables)))
main_label_keytypes <- c("ENSEMBL", "SYMBOL", "ALIAS")
main_label_db_keys <- stats::setNames(lapply(main_label_keytypes, function(k) {
  AnnotationDbi::keys(main_mouse_db, keytype = k)
}), main_label_keytypes)
main_symbol_mapping_audits <- vector("list", nrow(main_symbol_mapping_plan))
names(main_symbol_mapping_audits) <- main_symbol_mapping_plan$Table
main_symbol_summary_rows <- vector("list", nrow(main_symbol_mapping_plan))
main_symbol_route_rows <- vector("list", nrow(main_symbol_mapping_plan))
main_symbol_audit_files <- character(nrow(main_symbol_mapping_plan))

# Apply the same lookup policy to every source row regardless of its expression statistics.
for (i in seq_len(nrow(main_symbol_mapping_plan))) {
  table_key <- main_symbol_mapping_plan$Table[i]
  source_column <- main_symbol_mapping_plan$Source_column[i]
  result_table <- main_ranked_tables[[table_key]]
  stopifnot(source_column %in% names(result_table))
  source_ids <- as.character(result_table[[source_column]])
  usable_ids <- unique(source_ids[!is.na(source_ids) & nzchar(source_ids)])

  # Keep Ensembl-form labels on the Ensembl route, preserving their original text.
  ensembl_ids <- usable_ids[grepl("^ENSMUSG[0-9]{11}(\\.[0-9]+)?$", usable_ids)]
  symbol_ids <- setdiff(usable_ids, ensembl_ids)
  ensembl_map <- main_lookup_candidate_list(
    ensembl_ids, "ENSEMBL", main_label_db_keys$ENSEMBL, main_mouse_db
  )
  symbol_map <- main_lookup_candidate_list(
    symbol_ids, "SYMBOL", main_label_db_keys$SYMBOL, main_mouse_db
  )

  # Try exact recorded aliases only for names without any official-symbol candidate.
  alias_ids <- names(symbol_map)[lengths(symbol_map) == 0L]
  alias_map <- main_lookup_candidate_list(
    alias_ids, "ALIAS", main_label_db_keys$ALIAS, main_mouse_db
  )
  chosen_maps <- main_choose_label_candidates(
    usable_ids, ensembl_ids, ensembl_map, symbol_map, alias_map
  )
  audit <- main_make_row_id_audit(source_ids, chosen_maps$Candidates)
  audit$Lookup_route <- unname(chosen_maps$Route[match(source_ids, usable_ids)])
  audit$Lookup_route[is.na(source_ids) | !nzchar(source_ids)] <- "Missing_source_ID"
  audit$Original_record_ID <- source_ids
  audit$Source_column <- source_column
  stopifnot(
    nrow(audit) == nrow(result_table), identical(audit$Source_ID, source_ids),
    !anyNA(audit$Lookup_route)
  )
  main_symbol_mapping_audits[[table_key]] <- audit

  # Summarize candidate multiplicity separately from the route that found the match.
  main_symbol_summary_rows[[i]] <- data.frame(
    Study = main_symbol_mapping_plan$Study[i],
    Comparison = main_symbol_mapping_plan$Comparison[i],
    Source_rows = nrow(audit),
    One_match = sum(audit$Mapping_status == "One_match"),
    Multiple_matches = sum(audit$Mapping_status == "Multiple_matches"),
    No_match = sum(audit$Mapping_status == "No_match"),
    Missing_ID = sum(audit$Mapping_status == "Missing_source_ID"),
    Rows_sharing_target = sum(audit$Target_repeated)
  )
  main_symbol_route_rows[[i]] <- data.frame(
    Study = main_symbol_mapping_plan$Study[i],
    Comparison = main_symbol_mapping_plan$Comparison[i],
    Official_symbol_matches = sum(audit$Lookup_route == "Official_symbol" & audit$N_matches > 0L),
    Alias_matches = sum(audit$Lookup_route == "Alias" & audit$N_matches > 0L),
    Ensembl_matches = sum(audit$Lookup_route == "Ensembl" & audit$N_matches > 0L)
  )
  main_symbol_audit_files[i] <- file.path(
    main_analysis_dir, paste0(table_key, "_Entrez_mapping_audit.csv")
  )
  write.csv(audit, main_symbol_audit_files[i], row.names = FALSE)
}

# Display and save the four-contrast mapping and provenance summaries for review.
main_symbol_mapping_summary <- do.call(rbind, main_symbol_summary_rows)
main_symbol_route_summary <- do.call(rbind, main_symbol_route_rows)
stopifnot(
  with(main_symbol_mapping_summary,
    all(One_match + Multiple_matches + No_match + Missing_ID == Source_rows)),
  with(main_symbol_route_summary,
    all(Official_symbol_matches + Alias_matches + Ensembl_matches ==
      main_symbol_mapping_summary$One_match + main_symbol_mapping_summary$Multiple_matches))
)
print(main_symbol_mapping_summary, row.names = FALSE)
print(main_symbol_route_summary, row.names = FALSE)
main_symbol_summary_files <- file.path(
  main_analysis_dir,
  c("Gene_name_Entrez_mapping_summary.csv", "Gene_name_mapping_route_summary.csv")
)
write.csv(main_symbol_mapping_summary, main_symbol_summary_files[1], row.names = FALSE)
write.csv(main_symbol_route_summary, main_symbol_summary_files[2], row.names = FALSE)
print(file.exists(c(main_symbol_audit_files, main_symbol_summary_files)))
