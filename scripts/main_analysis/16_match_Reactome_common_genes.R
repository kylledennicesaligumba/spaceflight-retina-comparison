# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Match mouse Reactome pathways to the fixed common-gene universe.
# Apply the same inclusive 15-500 matched-gene rule before any enrichment result.
# https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/#gene-sets
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

# Match memberships by identifier and classify pathways using their matched sizes.
main_match_common_pathways <- function(pathways, common_ids, min_size = 15L, max_size = 500L) {
  stopifnot(
    is.list(pathways), length(pathways) > 0L,
    !is.null(names(pathways)), !anyNA(names(pathways)), all(nzchar(names(pathways))),
    anyDuplicated(names(pathways)) == 0L,
    all(vapply(pathways, is.character, logical(1))),
    !anyNA(unlist(pathways, use.names = FALSE)),
    is.character(common_ids), length(common_ids) > 0L,
    !anyNA(common_ids), all(nzchar(common_ids)), anyDuplicated(common_ids) == 0L,
    length(min_size) == 1L, length(max_size) == 1L,
    is.finite(min_size), is.finite(max_size), min_size >= 1L, max_size >= min_size
  )
  matched <- lapply(pathways, function(members) common_ids[common_ids %in% members])
  full_sizes <- vapply(pathways, function(members) length(unique(members)), integer(1))
  matched_sizes <- lengths(matched)
  eligible <- matched_sizes >= min_size & matched_sizes <= max_size
  status <- ifelse(full_sizes == 0L, "No_Entrez_members",
    ifelse(matched_sizes == 0L, "No_common_members",
      ifelse(matched_sizes < min_size, "Below_minimum",
        ifelse(matched_sizes > max_size, "Above_maximum", "Eligible"))))
  list(
    matched = matched, eligible = matched[eligible],
    full_sizes = full_sizes, matched_sizes = matched_sizes,
    keep = eligible, status = status
  )
}

# Restore the recorded resource and reviewed common-response checkpoint.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_pathway_response_file <- file.path(
  main_analysis_dir, "Main_analysis_common_responses_checkpoint_2026-11-06.rds"
)
main_pathway_resource_file <- file.path(
  main_analysis_dir, "Main_analysis_Reactome_mouse_resource_checkpoint_2026-11-06.rds"
)
stopifnot(file.exists(main_pathway_response_file), file.exists(main_pathway_resource_file))
main_pathway_response_context <- readRDS(main_pathway_response_file)
main_pathway_resource_context <- readRDS(main_pathway_resource_file)
main_pathway_common_ids <- main_pathway_response_context$common_ids
main_pathway_resource_catalog <- main_pathway_resource_context$catalog
stopifnot(
  length(main_pathway_response_context$checks) > 0L, all(main_pathway_response_context$checks),
  identical(main_pathway_common_ids, rownames(main_pathway_response_context$rank_scores)),
  identical(main_pathway_resource_context$species, "Mus musculus"),
  all(main_pathway_resource_catalog$Species == "Mus musculus"),
  identical(main_pathway_resource_catalog$Pathway_ID, names(main_pathway_resource_context$pathways)),
  identical(as.integer(main_pathway_resource_catalog$Entrez_genes_in_resource),
    unname(lengths(main_pathway_resource_context$pathways)))
)

# Fix the size limits on common-gene memberships, without selecting genes by FDR.
main_pathway_min_size <- 15L
main_pathway_max_size <- 500L
main_pathway_match <- main_match_common_pathways(
  main_pathway_resource_context$pathways, main_pathway_common_ids,
  min_size = main_pathway_min_size, max_size = main_pathway_max_size
)
main_reactome_common_pathways <- main_pathway_match$eligible
main_reactome_common_coverage <- main_pathway_resource_catalog
main_reactome_common_coverage$Common_gene_members <- unname(main_pathway_match$matched_sizes)
main_reactome_common_coverage$Eligible_for_enrichment <- unname(main_pathway_match$keep)
main_reactome_common_coverage$Eligibility_status <- unname(main_pathway_match$status)

# Record common-gene coverage while retaining every common gene in the ranking universe.
main_pathway_any_members <- unique(unlist(main_pathway_match$matched, use.names = FALSE))
main_pathway_eligible_members <- unique(unlist(main_reactome_common_pathways, use.names = FALSE))
main_reactome_common_gene_coverage <- data.frame(
  Entrez_ID = main_pathway_common_ids,
  In_any_mouse_Reactome_pathway = main_pathway_common_ids %in% main_pathway_any_members,
  In_eligible_pathway = main_pathway_common_ids %in% main_pathway_eligible_members,
  stringsAsFactors = FALSE
)
main_pathway_gene_summary <- c(
  Common_genes_in_rankings = length(main_pathway_common_ids),
  With_any_Reactome_membership = length(main_pathway_any_members),
  Without_Reactome_membership = length(main_pathway_common_ids) - length(main_pathway_any_members),
  With_membership_in_eligible_pathways = length(main_pathway_eligible_members)
)
main_pathway_status_levels <- c(
  "No_Entrez_members", "No_common_members", "Below_minimum", "Above_maximum", "Eligible"
)
main_pathway_eligibility_summary <- data.frame(
  Status = main_pathway_status_levels,
  Pathways = as.integer(table(factor(main_pathway_match$status, levels = main_pathway_status_levels)))
)
print(main_pathway_gene_summary)
print(main_pathway_eligibility_summary, row.names = FALSE)

# Verify the shared ranking universe, pathway keys and inclusive size limits.
main_pathway_match_checks <- c(
  Common_gene_order = identical(main_reactome_common_gene_coverage$Entrez_ID, main_pathway_common_ids),
  Common_members_only = all(main_pathway_any_members %in% main_pathway_common_ids),
  Eligible_sizes = length(main_reactome_common_pathways) > 0L && all(
    lengths(main_reactome_common_pathways) >= main_pathway_min_size &
      lengths(main_reactome_common_pathways) <= main_pathway_max_size),
  Pathway_order = identical(names(main_reactome_common_pathways),
    main_reactome_common_coverage$Pathway_ID[main_reactome_common_coverage$Eligible_for_enrichment]),
  All_pathways_accounted_for = sum(main_pathway_eligibility_summary$Pathways) == nrow(main_pathway_resource_catalog)
)
print(main_pathway_match_checks)
stopifnot(all(main_pathway_match_checks))

# Save coverage decisions and the fixed pathway memberships for later enrichment.
main_pathway_match_files <- file.path(main_analysis_dir, c(
  "Main_analysis_Reactome_common_pathway_coverage.csv",
  "Main_analysis_Reactome_common_gene_coverage.csv",
  "Main_analysis_Reactome_pathway_eligibility_summary.csv",
  "Main_analysis_Reactome_common_pathways_checkpoint_2026-11-06.rds"
))
utils::write.csv(main_reactome_common_coverage, main_pathway_match_files[1], row.names = FALSE)
utils::write.csv(main_reactome_common_gene_coverage, main_pathway_match_files[2], row.names = FALSE)
utils::write.csv(main_pathway_eligibility_summary, main_pathway_match_files[3], row.names = FALSE)
main_pathway_match_checkpoint <- list(
  pathways = main_reactome_common_pathways,
  coverage = main_reactome_common_coverage,
  gene_coverage = main_reactome_common_gene_coverage,
  gene_summary = main_pathway_gene_summary,
  eligibility_summary = main_pathway_eligibility_summary,
  checks = main_pathway_match_checks,
  common_ids = main_pathway_common_ids,
  size_rule = list(minSize = main_pathway_min_size, maxSize = main_pathway_max_size,
    count_basis = "Distinct pathway members among the accepted all-13 common Entrez IDs; inclusive limits"),
  ranking_universe = "All accepted common genes, including those without Reactome memberships; no FDR filtering",
  response_checkpoint = main_pathway_response_file,
  response_checkpoint_md5 = tools::md5sum(main_pathway_response_file),
  resource_checkpoint = main_pathway_resource_file,
  resource_checkpoint_md5 = tools::md5sum(main_pathway_resource_file),
  resource_versions = main_pathway_resource_context$versions,
  resource_metadata = main_pathway_resource_context$metadata,
  gene_mapping_versions = main_pathway_response_context$annotation_versions,
  stage = "Common-gene membership matching completed; enrichment has not been run"
)
saveRDS(main_pathway_match_checkpoint, main_pathway_match_files[4])
print(file.exists(main_pathway_match_files))
