# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Retrieve readable labels for the accepted common genes before pathway analysis.
# Preserve all returned labels and use Entrez IDs as the analysis keys.
# https://bioconductor.org/packages/release/bioc/manuals/AnnotationDbi/man/AnnotationDbi.pdf

# Restore the reviewed common-gene inputs without rerunning the earlier analyses.
main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_annotation_input_file <- file.path(
  main_analysis_dir, "Main_analysis_common_responses_checkpoint_2026-11-06.rds"
)
stopifnot(
  file.exists(main_annotation_input_file),
  requireNamespace("AnnotationDbi", quietly = TRUE),
  requireNamespace("org.Mm.eg.db", quietly = TRUE)
)
main_annotation_inputs <- readRDS(main_annotation_input_file)
main_annotation_ids <- main_annotation_inputs$common_ids
stopifnot(
  length(main_annotation_inputs$checks) > 0L,
  all(main_annotation_inputs$checks),
  is.character(main_annotation_ids), !anyNA(main_annotation_ids),
  length(main_annotation_ids) > 0L, anyDuplicated(main_annotation_ids) == 0L,
  identical(main_annotation_ids, rownames(main_annotation_inputs$rank_scores))
)

# Use the same mouse annotation snapshot that supported the reviewed ID mappings.
main_annotation_db <- org.Mm.eg.db::org.Mm.eg.db
main_annotation_saved_version <- main_annotation_inputs$annotation_versions$Version[
  main_annotation_inputs$annotation_versions$Component == "org.Mm.eg.db"
]
main_annotation_version <- as.character(utils::packageVersion("org.Mm.eg.db"))
stopifnot(identical(main_annotation_version, main_annotation_saved_version))

# Retrieve gene symbols and full names, retaining every annotation returned per ID.
main_common_gene_annotation_raw <- AnnotationDbi::select(
  main_annotation_db,
  keys = main_annotation_ids,
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)
stopifnot(
  !anyNA(main_common_gene_annotation_raw$ENTREZID),
  all(main_common_gene_annotation_raw$ENTREZID %in% main_annotation_ids),
  setequal(main_common_gene_annotation_raw$ENTREZID, main_annotation_ids)
)

# Count genes with labels and flag IDs that have more than one returned annotation row.
main_annotation_rows_per_id <- tabulate(
  match(main_common_gene_annotation_raw$ENTREZID, main_annotation_ids),
  nbins = length(main_annotation_ids)
)
main_annotation_has_symbol <- !is.na(main_common_gene_annotation_raw$SYMBOL) &
  nzchar(main_common_gene_annotation_raw$SYMBOL)
main_annotation_has_name <- !is.na(main_common_gene_annotation_raw$GENENAME) &
  nzchar(main_common_gene_annotation_raw$GENENAME)
print(c(
  Common_genes = length(main_annotation_ids),
  Genes_with_symbol = length(unique(main_common_gene_annotation_raw$ENTREZID[main_annotation_has_symbol])),
  Genes_with_full_name = length(unique(main_common_gene_annotation_raw$ENTREZID[main_annotation_has_name])),
  Genes_with_multiple_annotation_rows = sum(main_annotation_rows_per_id > 1L)
))
print(utils::head(main_common_gene_annotation_raw), row.names = FALSE)

# Save the complete lookup for review before attaching labels to response tables.
main_annotation_raw_file <- file.path(
  main_analysis_dir, "Main_analysis_common_gene_annotation_raw.csv"
)
utils::write.csv(main_common_gene_annotation_raw, main_annotation_raw_file, row.names = FALSE)
print(file.exists(main_annotation_raw_file))
