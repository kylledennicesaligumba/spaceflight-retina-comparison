# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Compare leading-edge genes for every reviewed shared 2-week glaucoma pathway.
# Keep all contributors and the original gene statistics, without gene-FDR filtering.
# https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/#running-a-leading-edge-analysis

# Compare two leading edges within the same pathway, using the frozen common Entrez IDs.
main_le_compare_members <- function(osd_ids, glaucoma_ids, pathway_members, common_ids) {
  stopifnot(
    is.character(osd_ids), is.character(glaucoma_ids),
    length(osd_ids) > 0L, length(glaucoma_ids) > 0L,
    !anyNA(osd_ids), !anyNA(glaucoma_ids),
    anyDuplicated(osd_ids) == 0L, anyDuplicated(glaucoma_ids) == 0L,
    all(osd_ids %in% pathway_members), all(glaucoma_ids %in% pathway_members),
    all(osd_ids %in% common_ids), all(glaucoma_ids %in% common_ids)
  )
  ids <- common_ids[common_ids %in% union(osd_ids, glaucoma_ids)]
  in_osd <- ids %in% osd_ids
  in_glaucoma <- ids %in% glaucoma_ids
  data.frame(
    Entrez_ID = ids, In_OSD_leading_edge = in_osd,
    In_Glaucoma_leading_edge = in_glaucoma,
    In_both_leading_edges = in_osd & in_glaucoma,
    Leading_edge_membership = ifelse(in_osd & in_glaucoma, "Both",
      ifelse(in_osd, "OSD_only", "Glaucoma_only")),
    stringsAsFactors = FALSE
  )
}

# Attach names and original gene-level responses by Entrez ID, preserving non-DEGs.
main_le_gene_details <- function(gene_ids, annotation, response, osd_key, glaucoma_key) {
  annotation_index <- match(gene_ids, annotation$Entrez_ID)
  gene_index <- match(gene_ids, response$common_ids)
  stopifnot(!anyNA(annotation_index), !anyNA(gene_index))
  result <- data.frame(
    Entrez_ID = gene_ids, Gene_symbol = annotation$Gene_symbol[annotation_index],
    Gene_name = annotation$Gene_name[annotation_index],
    OSD_logFC = unname(response$logFC[gene_index, osd_key]),
    Glaucoma_logFC = unname(response$logFC[gene_index, glaucoma_key]),
    OSD_rank_score = unname(response$rank_scores[gene_index, osd_key]),
    Glaucoma_rank_score = unname(response$rank_scores[gene_index, glaucoma_key]),
    OSD_gene_FDR = unname(response$FDR[gene_index, osd_key]),
    Glaucoma_gene_FDR = unname(response$FDR[gene_index, glaucoma_key]),
    stringsAsFactors = FALSE
  )
  result$Gene_logFC_pattern <- with(result, as.character(
    ifelse(OSD_logFC > 0 & Glaucoma_logFC > 0, "Higher in both",
      ifelse(OSD_logFC < 0 & Glaucoma_logFC < 0, "Lower in both",
        ifelse(OSD_logFC > 0 & Glaucoma_logFC < 0, "Higher OSD / lower glaucoma",
          ifelse(OSD_logFC < 0 & Glaucoma_logFC > 0, "Lower OSD / higher glaucoma", "Zero logFC"))))))
  result$Gene_FDR005_in_both <- result$OSD_gene_FDR < 0.05 & result$Glaucoma_gene_FDR < 0.05
  result
}

# Restore the frozen comparison, enrichment, response, annotation and pathway records.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_le_input_files <- file.path(main_analysis_dir, c(
  "Main_analysis_pathway_comparison_checkpoint_2026-11-06.rds",
  "Main_analysis_Reactome_fgsea_checkpoint_2026-11-06.rds",
  "Main_analysis_common_responses_checkpoint_2026-11-06.rds",
  "Main_analysis_common_gene_annotation_checkpoint_2026-11-06.rds",
  "Main_analysis_Reactome_common_pathways_checkpoint_2026-11-06.rds"
))
stopifnot(all(file.exists(main_le_input_files)))
main_le_comparison <- readRDS(main_le_input_files[1])
main_le_enrichment <- readRDS(main_le_input_files[2])
main_le_response <- readRDS(main_le_input_files[3])
main_le_annotation_context <- readRDS(main_le_input_files[4])
main_le_pathway_context <- readRDS(main_le_input_files[5])
main_le_annotation <- main_le_annotation_context$annotation
main_le_common_ids <- main_le_response$common_ids
main_le_osd_key <- "OSD758__FLT_0G_vs_HGC"
main_le_glaucoma_key <- "Glaucoma__2W_vs_Naive"
main_le_input_md5 <- tools::md5sum(main_le_input_files)
stopifnot(
  all(main_le_comparison$checks), all(main_le_response$checks),
  all(main_le_annotation_context$checks), all(main_le_pathway_context$checks),
  identical(unname(main_le_input_md5[2]), unname(main_le_comparison$enrichment_checkpoint_md5)),
  identical(unname(main_le_input_md5[3]), unname(main_le_annotation_context$response_checkpoint_md5)),
  identical(unname(main_le_input_md5[3]),
    unname(main_le_enrichment$input_md5[main_le_enrichment$response_checkpoint])),
  identical(unname(main_le_input_md5[5]),
    unname(main_le_enrichment$input_md5[main_le_enrichment$pathway_checkpoint])),
  identical(main_le_common_ids, main_le_comparison$common_gene_ids),
  identical(main_le_common_ids, main_le_enrichment$common_ids),
  identical(main_le_common_ids, main_le_pathway_context$common_ids),
  !anyNA(main_le_annotation$Entrez_ID), anyDuplicated(main_le_annotation$Entrez_ID) == 0L,
  setequal(main_le_annotation$Entrez_ID, main_le_common_ids),
  all(c(main_le_osd_key, main_le_glaucoma_key) %in% names(main_le_enrichment$results))
)
for (matrix_name in c("logFC", "rank_scores", "FDR")) {
  response_matrix <- main_le_response[[matrix_name]]
  stopifnot(
    identical(rownames(response_matrix), main_le_common_ids),
    all(c(main_le_osd_key, main_le_glaucoma_key) %in% colnames(response_matrix)),
    all(is.finite(response_matrix[, c(main_le_osd_key, main_le_glaucoma_key), drop = FALSE]))
  )
}

# Select all reviewed both-enriched pathways for the focused contrast, without choosing themes.
main_le_focus <- main_le_comparison$pathway_pairs[
  main_le_comparison$pathway_pairs$Study == "Glaucoma" &
    main_le_comparison$pathway_pairs$Comparison == "2W_vs_Naive" &
    main_le_comparison$pathway_pairs$Both_enriched, , drop = FALSE
]
main_le_focus <- main_le_focus[order(main_le_focus$OSD_NES < 0,
  main_le_focus$OSD_FDR, main_le_focus$Terrestrial_FDR, main_le_focus$Pathway_ID), , drop = FALSE]
main_le_expected <- main_le_comparison$overlap_summary[
  main_le_comparison$overlap_summary$Study == "Glaucoma" &
    main_le_comparison$overlap_summary$Comparison == "2W_vs_Naive", , drop = FALSE
]
stopifnot(
  nrow(main_le_expected) == 1L, nrow(main_le_focus) > 0L,
  nrow(main_le_focus) == main_le_expected$Both_FDR005,
  anyDuplicated(main_le_focus$Pathway_ID) == 0L,
  all(main_le_focus$NES_direction == "Same"),
  all(main_le_focus$OSD_FDR < 0.05 & main_le_focus$Terrestrial_FDR < 0.05)
)
main_le_osd_results <- main_le_enrichment$results[[main_le_osd_key]]
main_le_glaucoma_results <- main_le_enrichment$results[[main_le_glaucoma_key]]
main_le_osd_index <- match(main_le_focus$Pathway_ID, main_le_osd_results$Pathway_ID)
main_le_glaucoma_index <- match(main_le_focus$Pathway_ID, main_le_glaucoma_results$Pathway_ID)
stopifnot(
  !anyNA(main_le_osd_index), !anyNA(main_le_glaucoma_index),
  all(main_le_focus$Pathway_ID %in% names(main_le_pathway_context$pathways)),
  identical(main_le_focus$OSD_NES, main_le_osd_results$NES[main_le_osd_index]),
  identical(main_le_focus$Terrestrial_NES, main_le_glaucoma_results$NES[main_le_glaucoma_index])
)

# Compare leading edges within each pathway and retain the complete union of contributors.
main_le_membership_tables <- vector("list", nrow(main_le_focus))
main_le_overlap_rows <- vector("list", nrow(main_le_focus))
for (i in seq_len(nrow(main_le_focus))) {
  pathway_id <- main_le_focus$Pathway_ID[i]
  edges <- main_le_compare_members(
    main_le_osd_results$leadingEdge[[main_le_osd_index[i]]],
    main_le_glaucoma_results$leadingEdge[[main_le_glaucoma_index[i]]],
    main_le_pathway_context$pathways[[pathway_id]], main_le_common_ids
  )
  details <- main_le_gene_details(edges$Entrez_ID, main_le_annotation, main_le_response,
    main_le_osd_key, main_le_glaucoma_key)
  membership <- data.frame(
    Pathway_ID = pathway_id, Pathway_name = main_le_focus$Pathway_name[i],
    Pathway_NES_pattern = ifelse(main_le_focus$OSD_NES[i] > 0, "Positive in both", "Negative in both"),
    OSD_pathway_NES = main_le_focus$OSD_NES[i], Glaucoma_pathway_NES = main_le_focus$Terrestrial_NES[i],
    OSD_pathway_FDR = main_le_focus$OSD_FDR[i], Glaucoma_pathway_FDR = main_le_focus$Terrestrial_FDR[i],
    details, edges[, setdiff(names(edges), "Entrez_ID"), drop = FALSE],
    stringsAsFactors = FALSE, row.names = NULL
  )
  main_le_membership_tables[[i]] <- membership
  shared <- membership[membership$In_both_leading_edges, , drop = FALSE]
  main_le_overlap_rows[[i]] <- data.frame(
    Pathway_ID = pathway_id, Pathway_name = main_le_focus$Pathway_name[i],
    Pathway_NES_pattern = membership$Pathway_NES_pattern[1],
    OSD_leading_genes = sum(edges$In_OSD_leading_edge),
    Glaucoma_leading_genes = sum(edges$In_Glaucoma_leading_edge),
    Shared_leading_genes = nrow(shared),
    OSD_only_leading = sum(edges$Leading_edge_membership == "OSD_only"),
    Glaucoma_only_leading = sum(edges$Leading_edge_membership == "Glaucoma_only"),
    Shared_gene_symbols = paste(sort(shared$Gene_symbol), collapse = ";"),
    stringsAsFactors = FALSE
  )
}
main_le_membership <- do.call(rbind, main_le_membership_tables)
main_le_overlap_summary <- do.call(rbind, main_le_overlap_rows)
main_le_shared_rows <- main_le_membership[main_le_membership$In_both_leading_edges, , drop = FALSE]

# Count a gene as shared only if it belongs to both leading edges of the same pathway.
main_le_shared_gene_ids <- main_le_common_ids[main_le_common_ids %in% main_le_shared_rows$Entrez_ID]
main_le_shared_genes <- main_le_gene_details(main_le_shared_gene_ids, main_le_annotation,
  main_le_response, main_le_osd_key, main_le_glaucoma_key)
main_le_shared_genes$Shared_leading_pathways <- vapply(main_le_shared_gene_ids, function(id) {
  length(unique(main_le_shared_rows$Pathway_ID[main_le_shared_rows$Entrez_ID == id]))
}, integer(1))
main_le_shared_genes$Shared_pathway_IDs <- vapply(main_le_shared_gene_ids, function(id) {
  paste(unique(main_le_shared_rows$Pathway_ID[main_le_shared_rows$Entrez_ID == id]), collapse = ";")
}, character(1))
main_le_shared_genes$Shared_pathway_names <- vapply(main_le_shared_gene_ids, function(id) {
  paste(unique(main_le_shared_rows$Pathway_name[main_le_shared_rows$Entrez_ID == id]), collapse = ";")
}, character(1))
main_le_shared_genes <- main_le_shared_genes[
  order(main_le_shared_genes$Gene_symbol, main_le_shared_genes$Entrez_ID), , drop = FALSE
]
rownames(main_le_shared_genes) <- NULL

# Verify memberships, names, original values and overlap accounting before interpretation.
main_le_checks <- c(
  Reviewed_pathways = identical(main_le_overlap_summary$Pathway_ID, main_le_focus$Pathway_ID),
  Common_gene_members = all(main_le_membership$Entrez_ID %in% main_le_common_ids),
  Complete_gene_annotations = all(!is.na(main_le_membership$Gene_symbol) &
    nzchar(main_le_membership$Gene_symbol) & !is.na(main_le_membership$Gene_name) &
    nzchar(main_le_membership$Gene_name)),
  Finite_gene_values = all(is.finite(as.matrix(main_le_membership[, c(
    "OSD_logFC", "Glaucoma_logFC", "OSD_rank_score", "Glaucoma_rank_score",
    "OSD_gene_FDR", "Glaucoma_gene_FDR"), drop = FALSE]))),
  Overlap_accounting = all(main_le_overlap_summary$OSD_leading_genes ==
    main_le_overlap_summary$Shared_leading_genes + main_le_overlap_summary$OSD_only_leading) &&
    all(main_le_overlap_summary$Glaucoma_leading_genes ==
      main_le_overlap_summary$Shared_leading_genes + main_le_overlap_summary$Glaucoma_only_leading)
)
print(main_le_overlap_summary[, c("Pathway_name", "OSD_leading_genes",
  "Glaucoma_leading_genes", "Shared_leading_genes")], row.names = FALSE)
print(c(
  Shared_pathways_reviewed = nrow(main_le_overlap_summary),
  Unique_shared_leading_genes = nrow(main_le_shared_genes),
  Shared_genes_higher_in_both = sum(main_le_shared_genes$Gene_logFC_pattern == "Higher in both"),
  Shared_genes_lower_in_both = sum(main_le_shared_genes$Gene_logFC_pattern == "Lower in both"),
  Shared_genes_opposite_logFC = sum(main_le_shared_genes$Gene_logFC_pattern %in% c(
    "Higher OSD / lower glaucoma", "Lower OSD / higher glaucoma")),
  Shared_genes_zero_logFC = sum(main_le_shared_genes$Gene_logFC_pattern == "Zero logFC")
))
print(main_le_checks)
stopifnot(all(main_le_checks))

# Preview genes alphabetically within each direction group; these are not importance rankings.
for (pattern in c("Higher in both", "Lower in both", "Higher OSD / lower glaucoma",
    "Lower OSD / higher glaucoma", "Zero logFC")) {
  preview <- main_le_shared_genes[main_le_shared_genes$Gene_logFC_pattern == pattern,
    c("Gene_symbol", "OSD_logFC", "Glaucoma_logFC", "OSD_gene_FDR",
      "Glaucoma_gene_FDR", "Shared_leading_pathways"), drop = FALSE]
  if (nrow(preview) > 0L) {
    cat("\n", pattern, " — first 10 genes alphabetically\n", sep = "")
    print(utils::head(preview, 10L), digits = 4, row.names = FALSE)
  }
}

# Save every contributor, the shared-gene table and the complete matching provenance.
main_le_output_files <- file.path(main_analysis_dir, c(
  "Main_analysis_OSD758_Glaucoma_2W_leading_edge_overlap.csv",
  "Main_analysis_OSD758_Glaucoma_2W_leading_edge_membership.csv",
  "Main_analysis_OSD758_Glaucoma_2W_shared_leading_edge_genes.csv",
  "Main_analysis_OSD758_Glaucoma_2W_leading_edge_checkpoint_2026-11-06.rds"
))
utils::write.csv(main_le_overlap_summary, main_le_output_files[1], row.names = FALSE)
utils::write.csv(main_le_membership, main_le_output_files[2], row.names = FALSE)
utils::write.csv(main_le_shared_genes, main_le_output_files[3], row.names = FALSE)
main_le_checkpoint <- list(
  pathway_overlap = main_le_overlap_summary, membership = main_le_membership,
  shared_genes = main_le_shared_genes, focus_pathways = main_le_focus,
  checks = main_le_checks, input_files = main_le_input_files, input_md5 = main_le_input_md5,
  reference = main_le_osd_key, terrestrial = main_le_glaucoma_key,
  annotation_versions = main_le_annotation_context$annotation_versions,
  method = list(
    selection = "Every reviewed both-enriched 2W glaucoma pathway; no further theme selection",
    shared_gene = "Present in the leading edges of both responses for the same pathway",
    gene_statistics = "Original within-study logFC, signed rank score and gene FDR from the frozen common-response checkpoint",
    gene_filter = "No gene-FDR filter; gene-level significance retained as a descriptive flag",
    exclusivity = "OSD_only/Glaucoma_only refers to leading-edge membership in a pathway, not absence or no expression response",
    repeated_membership = "Pathways overlap; counts across sets are not independent evidence or causal gene rankings",
    preview = "First ten gene symbols alphabetically within each original logFC direction group",
    inference = "Descriptive leading-edge overlap; no overlap P values, joint FDR or causal mechanism test"
  )
)
saveRDS(main_le_checkpoint, main_le_output_files[4])
print(file.exists(main_le_output_files))
