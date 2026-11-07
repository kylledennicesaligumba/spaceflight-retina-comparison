# Prepare mouse gene-ID mapping to a shared Entrez Gene identifier system.
# Install the missing database and record its version before querying any study genes.
# https://bioconductor.org/packages/org.Mm.eg.db/

# Confirm that the reviewed result tables and required package manager are available.
stopifnot(
  exists("main_ranked_tables"), dir.exists(main_analysis_dir),
  requireNamespace("BiocManager", quietly = TRUE),
  requireNamespace("AnnotationDbi", quietly = TRUE)
)

# Install the mouse database through the Bioconductor release used by this R setup.
# The update argument leaves other installed packages out of the update step.
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db", ask = FALSE, update = FALSE)
}
stopifnot(requireNamespace("org.Mm.eg.db", quietly = TRUE))

# Load the mouse database and declare the common target identifier type.
main_mouse_db <- org.Mm.eg.db::org.Mm.eg.db
main_target_keytype <- "ENTREZID"

# Verify support for the source-label types seen in the user's previews.
main_needed_keytypes <- c("ENSEMBL", "SYMBOL", "ALIAS", "ENTREZID")
main_mapping_key_checks <- stats::setNames(
  main_needed_keytypes %in% AnnotationDbi::keytypes(main_mouse_db),
  main_needed_keytypes
)
print(main_mapping_key_checks)
stopifnot(all(main_mapping_key_checks))

# Record the actual software versions and annotation-source metadata.
main_mapping_versions <- data.frame(
  Component = c("R", "Bioconductor", "AnnotationDbi", "org.Mm.eg.db"),
  Version = c(
    as.character(getRversion()), as.character(BiocManager::version()),
    as.character(utils::packageVersion("AnnotationDbi")),
    as.character(utils::packageVersion("org.Mm.eg.db"))
  )
)
print(main_mapping_versions, row.names = FALSE)
main_mapping_metadata <- AnnotationDbi::metadata(main_mouse_db)

# Save the version and database-source records in the main-analysis folder.
main_mapping_setup_files <- file.path(
  main_analysis_dir,
  c("Main_analysis_mapping_versions.csv", "Main_analysis_mouse_annotation_metadata.csv")
)
write.csv(main_mapping_versions, main_mapping_setup_files[1], row.names = FALSE)
write.csv(main_mapping_metadata, main_mapping_setup_files[2], row.names = FALSE)
print(file.exists(main_mapping_setup_files))
