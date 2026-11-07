# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare one recorded mouse Reactome resource for all later pathway comparisons.
# Retrieve pathway definitions here; review their coverage before enrichment testing.
# https://bioconductor.org/packages/release/data/annotation/manuals/reactome.db/man/reactome.db.pdf
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

# Install any missing pathway packages through this R setup's Bioconductor release.
main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
stopifnot(
  dir.exists(main_analysis_dir),
  requireNamespace("BiocManager", quietly = TRUE),
  requireNamespace("AnnotationDbi", quietly = TRUE)
)
main_pathway_packages <- c("fgsea", "reactome.db")
main_missing_pathway_packages <- main_pathway_packages[!vapply(
  main_pathway_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(main_missing_pathway_packages) > 0L) {
  BiocManager::install(main_missing_pathway_packages, ask = FALSE, update = FALSE)
}
stopifnot(all(vapply(main_pathway_packages, requireNamespace, logical(1), quietly = TRUE)))

# Record the actual installed package versions and the Reactome source metadata.
main_reactome_versions <- data.frame(
  Component = c("R", "Bioconductor", "AnnotationDbi", "fgsea", "reactome.db"),
  Version = c(
    as.character(getRversion()), as.character(BiocManager::version()),
    as.character(utils::packageVersion("AnnotationDbi")),
    as.character(utils::packageVersion("fgsea")),
    as.character(utils::packageVersion("reactome.db"))
  ),
  stringsAsFactors = FALSE
)
main_reactome_metadata <- AnnotationDbi::metadata(reactome.db::reactome.db)
print(main_reactome_versions, row.names = FALSE)
print(main_reactome_metadata, row.names = FALSE)

# Select mouse pathway names explicitly from the multi-species database.
main_reactome_name_map <- as.list(reactome.db::reactomePATHID2NAME)
main_reactome_is_mouse <- vapply(main_reactome_name_map, function(x) {
  length(x) == 1L && !is.na(x) && startsWith(x, "Mus musculus: ")
}, logical(1))
main_reactome_mouse_ids <- names(main_reactome_name_map)[main_reactome_is_mouse]
stopifnot(length(main_reactome_mouse_ids) > 0L, anyDuplicated(main_reactome_mouse_ids) == 0L)
main_reactome_mouse_names <- vapply(
  main_reactome_name_map[main_reactome_mouse_ids], as.character, character(1)
)

# Retain each mouse pathway's distinct Entrez members, keyed by its pathway ID.
main_reactome_member_map <- as.list(reactome.db::reactomePATHID2EXTID)
main_reactome_mouse_pathways <- stats::setNames(
  lapply(main_reactome_mouse_ids, function(pathway_id) {
    members <- as.character(main_reactome_member_map[[pathway_id]])
    sort(unique(members[!is.na(members) & nzchar(members)]))
  }),
  main_reactome_mouse_ids
)
main_reactome_mouse_members <- unique(unlist(main_reactome_mouse_pathways, use.names = FALSE))
stopifnot(
  length(main_reactome_mouse_members) > 0L,
  all(grepl("^[1-9][0-9]*$", main_reactome_mouse_members))
)

# Preview the resource's pathway names and full gene counts before common-gene restriction.
main_reactome_mouse_catalog <- data.frame(
  Pathway_ID = main_reactome_mouse_ids,
  Pathway_name = sub("^Mus musculus: ", "", main_reactome_mouse_names),
  Species = "Mus musculus",
  Entrez_genes_in_resource = lengths(main_reactome_mouse_pathways),
  stringsAsFactors = FALSE, row.names = NULL
)
print(c(
  Mouse_pathways = nrow(main_reactome_mouse_catalog),
  Pathways_with_Entrez_members = sum(main_reactome_mouse_catalog$Entrez_genes_in_resource > 0L),
  Distinct_mouse_Entrez_genes = length(main_reactome_mouse_members)
))
print(utils::head(main_reactome_mouse_catalog), row.names = FALSE)

# Save the exact resource and plain source records so later comparisons reuse this snapshot.
main_reactome_resource_files <- file.path(
  main_analysis_dir,
  c(
    "Main_analysis_Reactome_versions.csv",
    "Main_analysis_Reactome_source_metadata.csv",
    "Main_analysis_Reactome_mouse_catalog.csv",
    "Main_analysis_Reactome_mouse_resource_checkpoint_2026-11-06.rds"
  )
)
utils::write.csv(main_reactome_versions, main_reactome_resource_files[1], row.names = FALSE)
utils::write.csv(main_reactome_metadata, main_reactome_resource_files[2], row.names = FALSE)
utils::write.csv(main_reactome_mouse_catalog, main_reactome_resource_files[3], row.names = FALSE)
main_reactome_resource_checkpoint <- list(
  pathways = main_reactome_mouse_pathways,
  catalog = main_reactome_mouse_catalog,
  versions = main_reactome_versions,
  metadata = main_reactome_metadata,
  species = "Mus musculus",
  species_selection = "Exact Mus musculus: prefix in reactomePATHID2NAME",
  identifier_type = "Entrez Gene IDs as provided by reactomePATHID2EXTID",
  collection = "All mouse pathways represented in the recorded reactome.db snapshot",
  membership_rule = "Distinct members within each pathway; retain membership in multiple pathways",
  stage = "Resource preparation only; common-gene coverage and enrichment settings remain next"
)
saveRDS(main_reactome_resource_checkpoint, main_reactome_resource_files[4])
print(file.exists(main_reactome_resource_files))
