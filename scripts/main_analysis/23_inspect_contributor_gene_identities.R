# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Read the saved contributor identities, fold changes and gene-FDR support.
# Use two reviewed programs as examples; retain all12 contrasts in the cross-model view.

# Restore the completed contributor checkpoint without recalculating any expression or enrichment results.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_identity_input <- file.path(main_analysis_dir,
  "Main_analysis_leading_edge_full_panel_checkpoint_2026-11-07.rds")
stopifnot(file.exists(main_identity_input))
main_identity_context <- readRDS(main_identity_input)
stopifnot(length(main_identity_context$checks) > 0L, all(main_identity_context$checks))
main_identity_genes <- main_identity_context$shared_genes
main_identity_membership <- main_identity_context$membership
main_identity_cross <- main_identity_context$shared_genes_by_pathway
main_identity_manifest <- main_identity_context$contrasts
main_identity_keys <- paste(main_identity_manifest$Study, main_identity_manifest$Comparison, sep = "__")
stopifnot(length(main_identity_keys) == 12L, anyDuplicated(main_identity_keys) == 0L,
  all(main_identity_keys %in% names(main_identity_cross)),
  identical(main_identity_genes$Gene_FDR005_in_both,
    main_identity_genes$OSD_gene_FDR < 0.05 & main_identity_genes$Terrestrial_gene_FDR < 0.05),
  anyDuplicated(paste(main_identity_genes$Study, main_identity_genes$Comparison,
    main_identity_genes$Entrez_ID, sep = "__")) == 0L)

# Display every shared contributor meeting both original gene-FDR cutoffs, across the complete panel.
main_identity_both <- main_identity_genes[main_identity_genes$Gene_FDR005_in_both, , drop = FALSE]
main_identity_both <- main_identity_both[order(match(paste(main_identity_both$Study,
  main_identity_both$Comparison, sep = "__"), main_identity_keys),
  main_identity_both$Gene_symbol, main_identity_both$Entrez_ID), , drop = FALSE]
rownames(main_identity_both) <- NULL
cat("Shared contributors meeting individual-gene FDR < 0.05 in both studies\n")
cat("Contrast-specific records:", nrow(main_identity_both),
  "| Distinct Entrez genes:", length(unique(main_identity_both$Entrez_ID)), "\n")
cat("Gene directions below are the original fold-change signs; names are alphabetical within contrast.\n")
for (key in main_identity_keys) {
  x <- main_identity_both[paste(main_identity_both$Study,
    main_identity_both$Comparison, sep = "__") == key, , drop = FALSE]
  if (nrow(x) > 0L) {
    cat("\n", x$Study[1], " | ", x$Comparison[1], "\n", sep = "")
    print(x[, c("Gene_symbol", "OSD_logFC", "Terrestrial_logFC", "OSD_gene_FDR",
      "Terrestrial_gene_FDR", "Gene_logFC_pattern")], digits = 4, row.names = FALSE)
  }
}
if (nrow(main_identity_both) == 0L) cat("No shared contributors meet both individual-gene cutoffs.\n")

# Locate the exact two programs already discussed, preserving their original pathway IDs.
main_identity_program_names <- c("Extracellular matrix organization", "Neuronal System")
main_identity_program_ids <- vapply(main_identity_program_names, function(name) {
  ids <- unique(main_identity_context$pathway_overlap$Pathway_ID[
    main_identity_context$pathway_overlap$Pathway_name == name])
  stopifnot(length(ids) == 1L, !is.na(ids))
  ids
}, character(1))
main_identity_focus_cross <- main_identity_cross[
  main_identity_cross$Pathway_ID %in% main_identity_program_ids, , drop = FALSE]
main_identity_focus_cross <- main_identity_focus_cross[order(
  match(main_identity_focus_cross$Pathway_ID, main_identity_program_ids),
  main_identity_focus_cross$Gene_symbol, main_identity_focus_cross$Entrez_ID), , drop = FALSE]
rownames(main_identity_focus_cross) <- NULL
main_identity_model_keys <- c(Glaucoma2W = "Glaucoma__2W_vs_Naive", RP = "RP__P23H_vs_WT",
  Uveitis = "Autoimmune_uveitis__EAU_vs_CON")
stopifnot(all(main_identity_model_keys %in% main_identity_keys))
main_identity_common_tables <- vector("list", length(main_identity_program_ids))
main_identity_program_summaries <- vector("list", length(main_identity_program_ids))
main_identity_original_checks <- logical(length(main_identity_program_ids))
main_identity_cross_checks <- logical(length(main_identity_program_ids))

# Intersect flight-model contributors within each SAME pathway, without a gene-FDR filter.
for (i in seq_along(main_identity_program_ids)) {
  pathway_id <- unname(main_identity_program_ids[i])
  members <- main_identity_membership[main_identity_membership$Pathway_ID == pathway_id, , drop = FALSE]
  shared_lists <- lapply(main_identity_model_keys, function(key) {
    m <- members[paste(members$Study, members$Comparison, sep = "__") == key &
      members$In_both_leading_edges, , drop = FALSE]
    stopifnot(anyDuplicated(m$Entrez_ID) == 0L)
    m$Entrez_ID
  })
  common_ids <- Reduce(intersect, shared_lists)
  main_identity_program_summaries[[i]] <- data.frame(Pathway_ID = pathway_id,
    Pathway_name = names(main_identity_program_ids)[i],
    Glaucoma2W_shared_LE = length(shared_lists$Glaucoma2W), RP_shared_LE = length(shared_lists$RP),
    Uveitis_shared_LE = length(shared_lists$Uveitis), Shared_in_all_three_pairs = length(common_ids),
    stringsAsFactors = FALSE)

  # Use the existing all-panel table as a separate check of the same-pathway intersection.
  cross <- main_identity_focus_cross[main_identity_focus_cross$Pathway_ID == pathway_id, , drop = FALSE]
  present <- matrix(vapply(unname(main_identity_model_keys), function(key)
    grepl("^Shared;", cross[[key]]), logical(nrow(cross))),
    nrow = nrow(cross), ncol = length(main_identity_model_keys))
  cross_common <- cross$Entrez_ID[rowSums(present) == 3L]
  main_identity_cross_checks[i] <- setequal(common_ids, cross_common)
  cross_index <- match(common_ids, cross$Entrez_ID)
  stopifnot(!anyNA(cross_index))
  wide <- cross[cross_index, c("Pathway_ID", "Pathway_name", "Entrez_ID", "Gene_symbol",
    "Gene_name", "OSD_logFC", "OSD_gene_FDR"), drop = FALSE]
  if (nrow(wide) > 0L) wide <- wide[order(wide$Gene_symbol, wide$Entrez_ID), , drop = FALSE]
  checks <- logical(length(main_identity_model_keys))
  for (j in seq_along(main_identity_model_keys)) {
    key <- main_identity_model_keys[j]
    prefix <- names(main_identity_model_keys)[j]
    m <- members[paste(members$Study, members$Comparison, sep = "__") == key, , drop = FALSE]
    m <- m[match(wide$Entrez_ID, m$Entrez_ID), , drop = FALSE]
    stopifnot(!anyNA(m$Entrez_ID), all(m$In_both_leading_edges))
    wide[[paste0(prefix, "_logFC")]] <- m$Terrestrial_logFC
    wide[[paste0(prefix, "_gene_FDR")]] <- m$Terrestrial_gene_FDR
    wide[[paste0(prefix, "_gene_direction")]] <- m$Gene_logFC_pattern
    wide[[paste0(prefix, "_gene_FDR_support")]] <- m$Gene_FDR_support
    wide[[paste0(prefix, "_pathway_direction")]] <- m$NES_direction
    wide[[paste0(prefix, "_pathway_NES")]] <- m$Terrestrial_pathway_NES
    wide[[paste0(prefix, "_pathway_FDR")]] <- m$Terrestrial_pathway_FDR
    checks[j] <- identical(wide$OSD_logFC, m$OSD_logFC) &&
      identical(wide$OSD_gene_FDR, m$OSD_gene_FDR) &&
      identical(wide$Gene_symbol, m$Gene_symbol) && identical(wide$Gene_name, m$Gene_name)
  }
  main_identity_original_checks[i] <- all(checks)
  main_identity_common_tables[[i]] <- wide
}
main_identity_program_summary <- do.call(rbind, main_identity_program_summaries)
main_identity_common <- do.call(rbind, main_identity_common_tables)
rownames(main_identity_program_summary) <- rownames(main_identity_common) <- NULL

# Check the display against the saved complete-panel gene-support counts and memberships.
main_identity_expected_support <- main_identity_context$gene_FDR_support[
  main_identity_context$gene_FDR_support$Pathway_group == "All selected pathways", , drop = FALSE]
support_index <- match(main_identity_keys, paste(main_identity_expected_support$Study,
  main_identity_expected_support$Comparison, sep = "__"))
stopifnot(!anyNA(support_index))
main_identity_both_counts <- vapply(main_identity_keys, function(key)
  sum(paste(main_identity_both$Study, main_identity_both$Comparison, sep = "__") == key), integer(1))
main_identity_checks <- c(
  Original_both_FDR_counts = all(main_identity_both_counts ==
    main_identity_expected_support$Both_gene_FDR005[support_index]),
  Same_pathway_intersections = all(main_identity_cross_checks) &&
    anyDuplicated(paste(main_identity_common$Pathway_ID, main_identity_common$Entrez_ID, sep = "__")) == 0L,
  Original_gene_values = all(main_identity_original_checks),
  All_panel_columns = all(main_identity_keys %in% names(main_identity_focus_cross)),
  Complete_accounting = sum(main_identity_program_summary$Shared_in_all_three_pairs) == nrow(main_identity_common)
)
cat("\nSame-pathway contributors shared across flight-glaucoma2W, flight-RP and flight-uveitis\n")
print(main_identity_program_summary[, c("Pathway_name", "Glaucoma2W_shared_LE", "RP_shared_LE",
  "Uveitis_shared_LE", "Shared_in_all_three_pairs")], row.names = FALSE)
cat("\nCommon-contributor previews retain genes regardless of individual gene FDR.\n")
for (pathway_id in unname(main_identity_program_ids)) {
  x <- main_identity_common[main_identity_common$Pathway_ID == pathway_id, , drop = FALSE]
  cat("\n", names(main_identity_program_ids)[match(pathway_id, main_identity_program_ids)], "\n", sep = "")
  if (nrow(x) == 0L) {
    cat("No gene is shared across all three flight-model pairs within this pathway.\n")
  } else {
    cat("First ten alphabetically, not ranked by importance; all genes are saved.\n")
    print(utils::head(x[, c("Gene_symbol", "OSD_logFC", "Glaucoma2W_logFC", "RP_logFC", "Uveitis_logFC")],
      10L), digits = 4, row.names = FALSE)
  }
}
print(main_identity_checks)
stopifnot(all(main_identity_checks))

# Save the inspected identities and their original support without changing prior selections or rankings.
main_identity_method <- list(
  purpose = "Inspect identities and original support already generated in the full-panel contributor review",
  both_gene_cutoffs = "View every selected shared contributor with original individual-gene FDR < .05 in each study; all12 contrasts considered",
  illustrative_programs = unname(main_identity_program_ids), program_names = main_identity_program_names,
  focus_rationale = "Two broad programs already reviewed in the heatmaps; examples of matrix and neuronal responses, not a new discovery test",
  three_model_examples = main_identity_model_keys,
  shared_across_pairs = "Intersect Entrez IDs within the SAME pathway across the three existing flight-model leading-edge overlaps",
  common_gene_filter = "No gene-FDR filtering of common contributors; preserve gene fold changes/FDR and pathway NES directions separately",
  full_panel = "Two-program wide table retains all12 contrast membership labels; complete all-pathway tables remain in module22",
  preview = "First ten common gene symbols alphabetically per pathway, not an importance ranking; full values saved",
  counting = "Both-cutoff records can repeat a gene across contrasts; common-program records can repeat a gene across pathways; counts are not independent replications",
  inference = "Descriptive inspection, not independent validation, joint FDR, causality, disease specificity or a new disease ranking",
  zero_overlap = "No common contributor across all three selected pairs does not erase pairwise sharing or pathway enrichment",
  provenance = "Only the saved module22 checkpoint is read; no model, enrichment, annotation or adjusted P value is recomputed"
)
main_identity_files <- file.path(main_analysis_dir, c(
  "Main_analysis_OSD758_shared_contributors_both_gene_FDR005.csv",
  "Main_analysis_OSD758_ECM_neuronal_cross_contrast_contributors.csv",
  "Main_analysis_OSD758_Glaucoma2W_RP_uveitis_common_program_genes.csv",
  "Main_analysis_OSD758_Glaucoma2W_RP_uveitis_common_program_summary.csv",
  "Main_analysis_contributor_identity_review_checkpoint_2026-11-07.rds",
  "Main_analysis_contributor_identity_review_method_2026-11-07.txt"
))
utils::write.csv(main_identity_both, main_identity_files[1], row.names = FALSE)
utils::write.csv(main_identity_focus_cross, main_identity_files[2], row.names = FALSE)
utils::write.csv(main_identity_common, main_identity_files[3], row.names = FALSE)
utils::write.csv(main_identity_program_summary, main_identity_files[4], row.names = FALSE)
saveRDS(list(both_gene_FDR = main_identity_both, focus_cross_contrasts = main_identity_focus_cross,
  common_program_genes = main_identity_common, program_summary = main_identity_program_summary,
  checks = main_identity_checks, method = main_identity_method, input_file = main_identity_input,
  input_md5 = tools::md5sum(main_identity_input)), main_identity_files[5])
utils::capture.output(print(main_identity_method), file = main_identity_files[6])
print(file.exists(main_identity_files))
