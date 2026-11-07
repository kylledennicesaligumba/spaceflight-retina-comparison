# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Align the reviewed gene labels to the accepted common-response gene order.
# Keep Entrez IDs as the unique keys and save labels for reuse in pathway outputs.

# Restore the saved response context and the raw annotation lookup reviewed by the user.
main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_annotation_response_file <- file.path(
  main_analysis_dir, "Main_analysis_common_responses_checkpoint_2026-11-06.rds"
)
main_annotation_lookup_file <- file.path(
  main_analysis_dir, "Main_analysis_common_gene_annotation_raw.csv"
)
stopifnot(file.exists(main_annotation_response_file), file.exists(main_annotation_lookup_file))
main_annotation_context <- readRDS(main_annotation_response_file)
main_annotation_lookup <- utils::read.csv(
  main_annotation_lookup_file, colClasses = "character", check.names = FALSE
)
main_annotation_common_ids <- main_annotation_context$common_ids
stopifnot(
  length(main_annotation_context$checks) > 0L, all(main_annotation_context$checks),
  identical(main_annotation_common_ids, rownames(main_annotation_context$rank_scores)),
  all(c("ENTREZID", "SYMBOL", "GENENAME") %in% names(main_annotation_lookup)),
  !anyNA(main_annotation_lookup$ENTREZID),
  anyDuplicated(main_annotation_lookup$ENTREZID) == 0L,
  setequal(main_annotation_lookup$ENTREZID, main_annotation_common_ids)
)

# Match by Entrez ID so labels follow the response matrices' exact gene order.
main_annotation_index <- match(main_annotation_common_ids, main_annotation_lookup$ENTREZID)
stopifnot(!anyNA(main_annotation_index))
main_common_gene_annotation <- data.frame(
  Entrez_ID = main_annotation_lookup$ENTREZID[main_annotation_index],
  Gene_symbol = main_annotation_lookup$SYMBOL[main_annotation_index],
  Gene_name = main_annotation_lookup$GENENAME[main_annotation_index],
  stringsAsFactors = FALSE, row.names = NULL
)

# Confirm complete labels and one ordered record per accepted common gene.
main_common_annotation_checks <- c(
  Gene_order = identical(main_common_gene_annotation$Entrez_ID, main_annotation_common_ids),
  Unique_IDs = anyDuplicated(main_common_gene_annotation$Entrez_ID) == 0L,
  Complete_symbols = all(!is.na(main_common_gene_annotation$Gene_symbol) &
    nzchar(main_common_gene_annotation$Gene_symbol)),
  Complete_names = all(!is.na(main_common_gene_annotation$Gene_name) &
    nzchar(main_common_gene_annotation$Gene_name))
)
print(dim(main_common_gene_annotation))
print(main_common_annotation_checks)
stopifnot(all(main_common_annotation_checks))

# Save the ordered annotation key and a plain-data checkpoint with its source records.
main_common_annotation_files <- file.path(
  main_analysis_dir,
  c(
    "Main_analysis_common_gene_annotation.csv",
    "Main_analysis_common_gene_annotation_checkpoint_2026-11-06.rds"
  )
)
utils::write.csv(main_common_gene_annotation, main_common_annotation_files[1], row.names = FALSE)
main_common_annotation_checkpoint <- list(
  annotation = main_common_gene_annotation,
  checks = main_common_annotation_checks,
  annotation_versions = main_annotation_context$annotation_versions,
  annotation_metadata = main_annotation_context$annotation_metadata,
  response_checkpoint = main_annotation_response_file,
  response_checkpoint_md5 = tools::md5sum(main_annotation_response_file),
  raw_lookup = main_annotation_lookup_file,
  raw_lookup_md5 = tools::md5sum(main_annotation_lookup_file),
  purpose = "Readable gene labels for the fixed common Entrez IDs; no response scores are changed"
)
saveRDS(main_common_annotation_checkpoint, main_common_annotation_files[2])
print(file.exists(main_common_annotation_files))
