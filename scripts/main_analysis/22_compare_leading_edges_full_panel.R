# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Compare contributors within every pathway enriched in flight and a terrestrial contrast.
# Keep all contributors, original gene statistics and separate same/opposite pathway directions.
# https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/#running-a-leading-edge-analysis

# Match leading-edge Entrez IDs within the same pathway, retaining the full union for review.
main_all_le_members <- function(osd_ids, terrestrial_ids, pathway_members, common_ids) {
  stopifnot(is.character(osd_ids), is.character(terrestrial_ids),
    length(osd_ids) > 0L, length(terrestrial_ids) > 0L,
    !anyNA(osd_ids), !anyNA(terrestrial_ids),
    anyDuplicated(osd_ids) == 0L, anyDuplicated(terrestrial_ids) == 0L,
    all(osd_ids %in% pathway_members), all(terrestrial_ids %in% pathway_members),
    all(osd_ids %in% common_ids), all(terrestrial_ids %in% common_ids))
  ids <- common_ids[common_ids %in% union(osd_ids, terrestrial_ids)]
  in_osd <- ids %in% osd_ids
  in_terrestrial <- ids %in% terrestrial_ids
  data.frame(Entrez_ID = ids, In_OSD_leading_edge = in_osd,
    In_Terrestrial_leading_edge = in_terrestrial,
    In_both_leading_edges = in_osd & in_terrestrial,
    Leading_edge_membership = ifelse(in_osd & in_terrestrial, "Both",
      ifelse(in_osd, "OSD_only", "Terrestrial_only")), stringsAsFactors = FALSE)
}

# Attach existing names, fold changes, rank scores and original gene-level FDRs by ID.
main_all_le_details <- function(ids, annotation, response, reference, terrestrial) {
  a <- match(ids, annotation$Entrez_ID)
  g <- match(ids, response$common_ids)
  stopifnot(!anyNA(a), !anyNA(g))
  x <- data.frame(Entrez_ID = ids, Gene_symbol = annotation$Gene_symbol[a],
    Gene_name = annotation$Gene_name[a], OSD_logFC = unname(response$logFC[g, reference]),
    Terrestrial_logFC = unname(response$logFC[g, terrestrial]),
    OSD_rank_score = unname(response$rank_scores[g, reference]),
    Terrestrial_rank_score = unname(response$rank_scores[g, terrestrial]),
    OSD_gene_FDR = unname(response$FDR[g, reference]),
    Terrestrial_gene_FDR = unname(response$FDR[g, terrestrial]), stringsAsFactors = FALSE)
  x$Gene_logFC_pattern <- with(x, as.character(ifelse(OSD_logFC > 0 & Terrestrial_logFC > 0,
    "Higher in both", ifelse(OSD_logFC < 0 & Terrestrial_logFC < 0, "Lower in both",
      ifelse(OSD_logFC > 0 & Terrestrial_logFC < 0, "Higher flight / lower terrestrial",
        ifelse(OSD_logFC < 0 & Terrestrial_logFC > 0,
          "Lower flight / higher terrestrial", "Zero logFC"))))))
  osd_sig <- x$OSD_gene_FDR < 0.05
  terrestrial_sig <- x$Terrestrial_gene_FDR < 0.05
  x$Gene_FDR_support <- as.character(ifelse(osd_sig & terrestrial_sig, "Both",
    ifelse(osd_sig, "Flight only", ifelse(terrestrial_sig, "Terrestrial only", "Neither"))))
  x$Gene_FDR005_in_both <- osd_sig & terrestrial_sig
  x
}

# Preserve typed empty tables when a contrast has no eligible pathways or no shared contributors.
main_all_le_bind <- function(tables, template) {
  tables <- Filter(function(x) nrow(x) > 0L, tables)
  if (length(tables) == 0L) return(template)
  x <- do.call(rbind, tables)
  rownames(x) <- NULL
  x
}

# Restore the frozen records and the previously reviewed glaucoma result, without rerunning models.
if (!exists("main_analysis_dir")) {
  main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
}
main_all_le_files_in <- file.path(main_analysis_dir, c(
  "Main_analysis_pathway_comparison_checkpoint_2026-11-06.rds",
  "Main_analysis_Reactome_fgsea_checkpoint_2026-11-06.rds",
  "Main_analysis_common_responses_checkpoint_2026-11-06.rds",
  "Main_analysis_common_gene_annotation_checkpoint_2026-11-06.rds",
  "Main_analysis_Reactome_common_pathways_checkpoint_2026-11-06.rds",
  "Main_analysis_OSD758_Glaucoma_2W_leading_edge_checkpoint_2026-11-06.rds"
))
stopifnot(all(file.exists(main_all_le_files_in)))
main_all_le_input_md5 <- tools::md5sum(main_all_le_files_in)
main_all_le_comparison <- readRDS(main_all_le_files_in[1])
main_all_le_enrichment <- readRDS(main_all_le_files_in[2])
main_all_le_response <- readRDS(main_all_le_files_in[3])
main_all_le_annotation_context <- readRDS(main_all_le_files_in[4])
main_all_le_pathway_context <- readRDS(main_all_le_files_in[5])
main_all_le_previous_glaucoma <- readRDS(main_all_le_files_in[6])
main_all_le_reference <- "OSD758__FLT_0G_vs_HGC"
main_all_le_annotation <- main_all_le_annotation_context$annotation
main_all_le_ids <- main_all_le_response$common_ids
main_all_le_manifest <- main_all_le_comparison$contrasts
main_all_le_manifest_keys <- paste(main_all_le_manifest$Study, main_all_le_manifest$Comparison, sep = "__")
main_all_le_terrestrial <- main_all_le_manifest[main_all_le_manifest$Role == "Terrestrial", , drop = FALSE]
main_all_le_keys <- paste(main_all_le_terrestrial$Study, main_all_le_terrestrial$Comparison, sep = "__")
stopifnot(length(main_all_le_keys) == 12L, anyDuplicated(main_all_le_manifest_keys) == 0L,
  sum(main_all_le_manifest$Role == "Reference") == 1L,
  main_all_le_manifest$Role[match(main_all_le_reference, main_all_le_manifest_keys)] == "Reference",
  setequal(c(main_all_le_reference, main_all_le_keys), main_all_le_manifest_keys),
  setequal(names(main_all_le_enrichment$results), main_all_le_manifest_keys),
  all(main_all_le_comparison$checks), all(main_all_le_response$checks),
  all(main_all_le_annotation_context$checks), all(main_all_le_pathway_context$checks),
  all(main_all_le_previous_glaucoma$checks),
  identical(unname(main_all_le_input_md5[2]), unname(main_all_le_comparison$enrichment_checkpoint_md5)),
  identical(unname(main_all_le_input_md5[3]), unname(main_all_le_annotation_context$response_checkpoint_md5)),
  identical(unname(main_all_le_input_md5[3]),
    unname(main_all_le_enrichment$input_md5[main_all_le_enrichment$response_checkpoint])),
  identical(unname(main_all_le_input_md5[5]),
    unname(main_all_le_enrichment$input_md5[main_all_le_enrichment$pathway_checkpoint])),
  identical(unname(main_all_le_input_md5[1:5]),
    unname(main_all_le_previous_glaucoma$input_md5[main_all_le_files_in[1:5]])),
  identical(main_all_le_ids, main_all_le_comparison$common_gene_ids),
  identical(main_all_le_ids, main_all_le_enrichment$common_ids),
  identical(main_all_le_ids, main_all_le_pathway_context$common_ids),
  !anyNA(main_all_le_ids), anyDuplicated(main_all_le_ids) == 0L,
  !anyNA(main_all_le_annotation$Entrez_ID), anyDuplicated(main_all_le_annotation$Entrez_ID) == 0L,
  setequal(main_all_le_annotation$Entrez_ID, main_all_le_ids))
for (field in c("logFC", "rank_scores", "FDR")) {
  x <- main_all_le_response[[field]]
  stopifnot(is.matrix(x), identical(rownames(x), main_all_le_ids),
    identical(colnames(x), main_all_le_manifest_keys), all(is.finite(x)))
}
stopifnot(all(main_all_le_response$FDR >= 0 & main_all_le_response$FDR <= 1))

# Retain all pathways passing the original FDR cutoff in both responses, including opposing signs.
main_all_le_pairs <- main_all_le_comparison$pathway_pairs
stopifnot(identical(main_all_le_pairs$Both_enriched,
  main_all_le_pairs$OSD_FDR < 0.05 & main_all_le_pairs$Terrestrial_FDR < 0.05))
main_all_le_focus <- main_all_le_pairs[main_all_le_pairs$Both_enriched, , drop = FALSE]
main_all_le_focus_keys <- paste(main_all_le_focus$Study, main_all_le_focus$Comparison, sep = "__")
stopifnot(all(main_all_le_focus_keys %in% main_all_le_keys))
main_all_le_focus <- main_all_le_focus[order(match(main_all_le_focus_keys, main_all_le_keys),
  match(main_all_le_focus$NES_direction, c("Same", "Opposite", "Zero_NES")),
  main_all_le_focus$OSD_FDR, main_all_le_focus$Terrestrial_FDR, main_all_le_focus$Pathway_ID), , drop = FALSE]
rownames(main_all_le_focus) <- NULL
main_all_le_empty_details <- main_all_le_details(character(), main_all_le_annotation,
  main_all_le_response, main_all_le_reference, main_all_le_keys[1])
main_all_le_membership_template <- data.frame(Study = character(), Comparison = character(),
  Pathway_ID = character(), Pathway_name = character(), NES_direction = character(),
  OSD_pathway_NES = double(), Terrestrial_pathway_NES = double(),
  OSD_pathway_FDR = double(), Terrestrial_pathway_FDR = double(), main_all_le_empty_details,
  In_OSD_leading_edge = logical(), In_Terrestrial_leading_edge = logical(),
  In_both_leading_edges = logical(), Leading_edge_membership = character(), stringsAsFactors = FALSE)
main_all_le_overlap_template <- data.frame(Study = character(), Comparison = character(),
  Pathway_ID = character(), Pathway_name = character(), NES_direction = character(),
  OSD_NES = double(), Terrestrial_NES = double(), OSD_FDR = double(), Terrestrial_FDR = double(),
  OSD_leading_genes = integer(), Terrestrial_leading_genes = integer(), Shared_leading_genes = integer(),
  OSD_only_leading = integer(), Terrestrial_only_leading = integer(),
  Shared_gene_symbols = character(), stringsAsFactors = FALSE)
main_all_le_membership_tables <- vector("list", nrow(main_all_le_focus))
main_all_le_overlap_tables <- vector("list", nrow(main_all_le_focus))
main_all_le_reference_results <- main_all_le_enrichment$results[[main_all_le_reference]]

# Compare each pair within its own pathway before attaching that contrast's gene responses.
for (i in seq_len(nrow(main_all_le_focus))) {
  pair <- main_all_le_focus[i, , drop = FALSE]
  key <- paste(pair$Study, pair$Comparison, sep = "__")
  result <- main_all_le_enrichment$results[[key]]
  o <- match(pair$Pathway_ID, main_all_le_reference_results$Pathway_ID)
  t <- match(pair$Pathway_ID, result$Pathway_ID)
  stopifnot(!anyNA(c(o, t)), pair$Pathway_ID %in% names(main_all_le_pathway_context$pathways),
    pair$OSD_NES == main_all_le_reference_results$NES[o], pair$OSD_FDR == main_all_le_reference_results$FDR[o],
    pair$Terrestrial_NES == result$NES[t], pair$Terrestrial_FDR == result$FDR[t])
  expected_direction <- ifelse(pair$OSD_NES * pair$Terrestrial_NES > 0, "Same",
    ifelse(pair$OSD_NES * pair$Terrestrial_NES < 0, "Opposite", "Zero_NES"))
  stopifnot(pair$NES_direction == expected_direction)
  edges <- main_all_le_members(main_all_le_reference_results$leadingEdge[[o]], result$leadingEdge[[t]],
    main_all_le_pathway_context$pathways[[pair$Pathway_ID]], main_all_le_ids)
  details <- main_all_le_details(edges$Entrez_ID, main_all_le_annotation,
    main_all_le_response, main_all_le_reference, key)
  main_all_le_membership_tables[[i]] <- data.frame(Study = pair$Study, Comparison = pair$Comparison,
    Pathway_ID = pair$Pathway_ID, Pathway_name = pair$Pathway_name, NES_direction = pair$NES_direction,
    OSD_pathway_NES = pair$OSD_NES, Terrestrial_pathway_NES = pair$Terrestrial_NES,
    OSD_pathway_FDR = pair$OSD_FDR, Terrestrial_pathway_FDR = pair$Terrestrial_FDR,
    details, edges[, setdiff(names(edges), "Entrez_ID"), drop = FALSE], stringsAsFactors = FALSE)
  main_all_le_overlap_tables[[i]] <- data.frame(Study = pair$Study, Comparison = pair$Comparison,
    Pathway_ID = pair$Pathway_ID, Pathway_name = pair$Pathway_name, NES_direction = pair$NES_direction,
    OSD_NES = pair$OSD_NES, Terrestrial_NES = pair$Terrestrial_NES,
    OSD_FDR = pair$OSD_FDR, Terrestrial_FDR = pair$Terrestrial_FDR,
    OSD_leading_genes = sum(edges$In_OSD_leading_edge),
    Terrestrial_leading_genes = sum(edges$In_Terrestrial_leading_edge),
    Shared_leading_genes = sum(edges$In_both_leading_edges),
    OSD_only_leading = sum(edges$Leading_edge_membership == "OSD_only"),
    Terrestrial_only_leading = sum(edges$Leading_edge_membership == "Terrestrial_only"),
    Shared_gene_symbols = paste(sort(details$Gene_symbol[edges$In_both_leading_edges]), collapse = ";"),
    stringsAsFactors = FALSE)
}
main_all_le_membership <- main_all_le_bind(main_all_le_membership_tables, main_all_le_membership_template)
main_all_le_overlap <- main_all_le_bind(main_all_le_overlap_tables, main_all_le_overlap_template)
main_all_le_shared_rows <- main_all_le_membership[main_all_le_membership$In_both_leading_edges, , drop = FALSE]

# Deduplicate shared genes within each contrast while retaining their same/opposite pathway links.
main_all_le_gene_template <- data.frame(Study = character(), Comparison = character(),
  main_all_le_empty_details, Shared_pathways_same_NES = integer(), Shared_pathways_opposite_NES = integer(),
  Shared_pathways_zero_NES = integer(), Shared_pathways_total = integer(),
  Same_NES_pathway_IDs = character(), Opposite_NES_pathway_IDs = character(),
  Zero_NES_pathway_IDs = character(), Shared_pathway_names = character(), stringsAsFactors = FALSE)
main_all_le_gene_tables <- vector("list", length(main_all_le_keys))
main_all_le_summary_tables <- vector("list", length(main_all_le_keys))
main_all_le_support_tables <- list()
for (j in seq_along(main_all_le_keys)) {
  study <- main_all_le_terrestrial$Study[j]
  comparison <- main_all_le_terrestrial$Comparison[j]
  rows <- main_all_le_shared_rows[main_all_le_shared_rows$Study == study &
    main_all_le_shared_rows$Comparison == comparison, , drop = FALSE]
  genes <- main_all_le_ids[main_all_le_ids %in% rows$Entrez_ID]
  details <- main_all_le_details(genes, main_all_le_annotation, main_all_le_response,
    main_all_le_reference, main_all_le_keys[j])
  gene_table <- data.frame(Study = rep(study, length(genes)),
    Comparison = rep(comparison, length(genes)), details, stringsAsFactors = FALSE)
  for (direction in c("Same", "Opposite", "Zero_NES")) {
    suffix <- switch(direction, Same = "same_NES", Opposite = "opposite_NES", Zero_NES = "zero_NES")
    gene_table[[paste0("Shared_pathways_", suffix)]] <- vapply(genes, function(id)
      length(unique(rows$Pathway_ID[rows$Entrez_ID == id & rows$NES_direction == direction])), integer(1))
  }
  gene_table$Shared_pathways_total <- vapply(genes, function(id)
    length(unique(rows$Pathway_ID[rows$Entrez_ID == id])), integer(1))
  for (direction in c("Same", "Opposite", "Zero_NES")) {
    label <- switch(direction, Same = "Same", Opposite = "Opposite", Zero_NES = "Zero")
    gene_table[[paste0(label, "_NES_pathway_IDs")]] <- vapply(genes, function(id)
      paste(unique(rows$Pathway_ID[rows$Entrez_ID == id & rows$NES_direction == direction]), collapse = ";"), character(1))
  }
  gene_table$Shared_pathway_names <- vapply(genes, function(id)
    paste(unique(rows$Pathway_name[rows$Entrez_ID == id]), collapse = ";"), character(1))
  gene_table <- gene_table[order(gene_table$Gene_symbol, gene_table$Entrez_ID),
    names(main_all_le_gene_template), drop = FALSE]
  main_all_le_gene_tables[[j]] <- gene_table
  pathways <- main_all_le_overlap[main_all_le_overlap$Study == study &
    main_all_le_overlap$Comparison == comparison, , drop = FALSE]
  main_all_le_summary_tables[[j]] <- data.frame(Study = study, Comparison = comparison,
    Both_enriched_pathways = nrow(pathways), Same_NES_pathways = sum(pathways$NES_direction == "Same"),
    Opposite_NES_pathways = sum(pathways$NES_direction == "Opposite"),
    Zero_NES_pathways = sum(pathways$NES_direction == "Zero_NES"),
    Unique_shared_LE_genes = nrow(gene_table),
    Shared_genes_same_NES = sum(gene_table$Shared_pathways_same_NES > 0),
    Shared_genes_opposite_NES = sum(gene_table$Shared_pathways_opposite_NES > 0), stringsAsFactors = FALSE)

  # Count original gene FDR support without removing contributors or readjusting P values.
  for (group in c("All selected pathways", "Same NES pathways", "Opposite NES pathways", "Zero NES pathways")) {
    keep <- switch(group, "All selected pathways" = rep(TRUE, nrow(gene_table)),
      "Same NES pathways" = gene_table$Shared_pathways_same_NES > 0,
      "Opposite NES pathways" = gene_table$Shared_pathways_opposite_NES > 0,
      "Zero NES pathways" = gene_table$Shared_pathways_zero_NES > 0)
    x <- gene_table[keep, , drop = FALSE]
    main_all_le_support_tables[[length(main_all_le_support_tables) + 1L]] <- data.frame(
      Study = study, Comparison = comparison, Pathway_group = group, Shared_genes = nrow(x),
      Both_gene_FDR005 = sum(x$Gene_FDR_support == "Both"),
      Flight_only_gene_FDR005 = sum(x$Gene_FDR_support == "Flight only"),
      Terrestrial_only_gene_FDR005 = sum(x$Gene_FDR_support == "Terrestrial only"),
      Neither_gene_FDR005 = sum(x$Gene_FDR_support == "Neither"),
      Higher_in_both = sum(x$Gene_logFC_pattern == "Higher in both"),
      Lower_in_both = sum(x$Gene_logFC_pattern == "Lower in both"),
      Opposite_logFC = sum(x$Gene_logFC_pattern %in% c("Higher flight / lower terrestrial",
        "Lower flight / higher terrestrial")), Zero_logFC = sum(x$Gene_logFC_pattern == "Zero logFC"),
      stringsAsFactors = FALSE)
  }
}
main_all_le_genes <- main_all_le_bind(main_all_le_gene_tables, main_all_le_gene_template)
main_all_le_summary <- do.call(rbind, main_all_le_summary_tables)
main_all_le_support <- do.call(rbind, main_all_le_support_tables)
rownames(main_all_le_summary) <- rownames(main_all_le_support) <- NULL

# Compare identities across models only for the same pathway and the same Entrez gene.
main_all_le_cross_key <- paste(main_all_le_shared_rows$Pathway_ID, main_all_le_shared_rows$Entrez_ID, sep = "__")
main_all_le_cross <- main_all_le_shared_rows[!duplicated(main_all_le_cross_key),
  c("Pathway_ID", "Pathway_name", "Entrez_ID", "Gene_symbol", "Gene_name", "OSD_logFC", "OSD_gene_FDR"), drop = FALSE]
main_all_le_cross <- main_all_le_cross[order(match(main_all_le_cross$Pathway_ID,
  names(main_all_le_pathway_context$pathways)), main_all_le_cross$Gene_symbol, main_all_le_cross$Entrez_ID), , drop = FALSE]
rownames(main_all_le_cross) <- NULL
main_all_le_cross_ids <- paste(main_all_le_cross$Pathway_ID, main_all_le_cross$Entrez_ID, sep = "__")
main_all_le_cross_matrix <- matrix("Pathway not both-enriched", nrow(main_all_le_cross), length(main_all_le_keys),
  dimnames = list(main_all_le_cross_ids, main_all_le_keys))
for (j in seq_along(main_all_le_keys)) {
  study <- main_all_le_terrestrial$Study[j]
  comparison <- main_all_le_terrestrial$Comparison[j]
  pathways <- main_all_le_focus[main_all_le_focus$Study == study &
    main_all_le_focus$Comparison == comparison, , drop = FALSE]
  p <- match(main_all_le_cross$Pathway_ID, pathways$Pathway_ID)
  eligible <- !is.na(p)
  if (any(eligible)) {
    main_all_le_cross_matrix[eligible, j] <- paste("Not shared;", pathways$NES_direction[p[eligible]], "NES")
  }
  shared <- main_all_le_shared_rows[main_all_le_shared_rows$Study == study &
    main_all_le_shared_rows$Comparison == comparison, , drop = FALSE]
  matched <- match(paste(shared$Pathway_ID, shared$Entrez_ID, sep = "__"), main_all_le_cross_ids)
  stopifnot(!anyNA(matched))
  if (length(matched) > 0L) {
    main_all_le_cross_matrix[matched, j] <- paste("Shared;", shared$NES_direction, "NES")
  }
}
main_all_le_cross$Shared_same_NES_contrasts <- rowSums(main_all_le_cross_matrix == "Shared; Same NES")
main_all_le_cross$Shared_opposite_NES_contrasts <- rowSums(main_all_le_cross_matrix == "Shared; Opposite NES")
main_all_le_cross <- data.frame(main_all_le_cross, main_all_le_cross_matrix, check.names = FALSE)

# Verify contrast coverage, original values, gene accounting and the earlier glaucoma comparison.
main_all_le_expected <- main_all_le_comparison$overlap_summary
expected_index <- match(main_all_le_keys, paste(main_all_le_expected$Study,
  main_all_le_expected$Comparison, sep = "__"))
stopifnot(!anyNA(expected_index))
main_all_le_expected <- main_all_le_expected[expected_index, , drop = FALSE]
gene_index <- match(main_all_le_membership$Entrez_ID, main_all_le_ids)
contrast_index <- match(paste(main_all_le_membership$Study,
  main_all_le_membership$Comparison, sep = "__"), main_all_le_manifest_keys)
annotation_index <- match(main_all_le_membership$Entrez_ID, main_all_le_annotation$Entrez_ID)
original_checks <- vapply(c("logFC", "rank_scores", "FDR"), function(field) {
  suffix <- switch(field, logFC = "logFC", rank_scores = "rank_score", FDR = "gene_FDR")
  identical(main_all_le_membership[[paste0("OSD_", suffix)]],
    unname(main_all_le_response[[field]][gene_index, main_all_le_reference])) &&
    identical(main_all_le_membership[[paste0("Terrestrial_", suffix)]],
      unname(main_all_le_response[[field]][cbind(gene_index, contrast_index)]))
}, logical(1))
old_genes <- main_all_le_previous_glaucoma$shared_genes
new_genes <- main_all_le_genes[main_all_le_genes$Study == "Glaucoma" &
  main_all_le_genes$Comparison == "2W_vs_Naive", , drop = FALSE]
g <- match(old_genes$Entrez_ID, new_genes$Entrez_ID)
old_pathways <- main_all_le_previous_glaucoma$pathway_overlap
new_pathways <- main_all_le_overlap[main_all_le_overlap$Study == "Glaucoma" &
  main_all_le_overlap$Comparison == "2W_vs_Naive", , drop = FALSE]
p <- match(old_pathways$Pathway_ID, new_pathways$Pathway_ID)
main_all_le_glaucoma_checks <- c(!anyNA(g), !anyNA(p), nrow(old_genes) == nrow(new_genes),
  nrow(old_pathways) == nrow(new_pathways),
  identical(old_genes$OSD_logFC, new_genes$OSD_logFC[g]),
  identical(old_genes$Glaucoma_logFC, new_genes$Terrestrial_logFC[g]),
  identical(old_genes$OSD_gene_FDR, new_genes$OSD_gene_FDR[g]),
  identical(old_genes$Glaucoma_gene_FDR, new_genes$Terrestrial_gene_FDR[g]),
  identical(old_genes$Gene_FDR005_in_both, new_genes$Gene_FDR005_in_both[g]),
  identical(old_pathways$Shared_leading_genes, new_pathways$Shared_leading_genes[p]))
cross_counts <- vapply(main_all_le_keys, function(key) sum(grepl("^Shared;", main_all_le_cross_matrix[, key])), integer(1))
shared_row_counts <- vapply(main_all_le_keys, function(key) sum(paste(main_all_le_shared_rows$Study,
  main_all_le_shared_rows$Comparison, sep = "__") == key), integer(1))
main_all_le_checks <- c(
  All_contrasts = identical(paste(main_all_le_summary$Study, main_all_le_summary$Comparison, sep = "__"), main_all_le_keys),
  Reviewed_pathways = all(main_all_le_summary$Both_enriched_pathways == main_all_le_expected$Both_FDR005) &&
    all(main_all_le_summary$Same_NES_pathways == main_all_le_expected$Same_direction_both) &&
    all(main_all_le_summary$Opposite_NES_pathways == main_all_le_expected$Opposite_direction_both) &&
    all(main_all_le_summary$Zero_NES_pathways == main_all_le_expected$Zero_NES_both),
  Same_pathway_members = all(vapply(unique(main_all_le_membership$Pathway_ID), function(pathway_id)
    all(main_all_le_membership$Entrez_ID[main_all_le_membership$Pathway_ID == pathway_id] %in%
      main_all_le_pathway_context$pathways[[pathway_id]]), logical(1))),
  Original_gene_values = all(original_checks),
  Complete_annotations = identical(main_all_le_membership$Gene_symbol, main_all_le_annotation$Gene_symbol[annotation_index]) &&
    identical(main_all_le_membership$Gene_name, main_all_le_annotation$Gene_name[annotation_index]) &&
    all(!is.na(main_all_le_membership$Gene_symbol) & nzchar(main_all_le_membership$Gene_symbol)) &&
    all(!is.na(main_all_le_membership$Gene_name) & nzchar(main_all_le_membership$Gene_name)),
  Support_accounting = all(rowSums(main_all_le_support[, c("Both_gene_FDR005", "Flight_only_gene_FDR005",
    "Terrestrial_only_gene_FDR005", "Neither_gene_FDR005")]) == main_all_le_support$Shared_genes) &&
    all(rowSums(main_all_le_support[, c("Higher_in_both", "Lower_in_both", "Opposite_logFC", "Zero_logFC")]) ==
      main_all_le_support$Shared_genes) &&
    all(main_all_le_overlap$Shared_leading_genes + main_all_le_overlap$OSD_only_leading == main_all_le_overlap$OSD_leading_genes) &&
    all(main_all_le_overlap$Shared_leading_genes + main_all_le_overlap$Terrestrial_only_leading == main_all_le_overlap$Terrestrial_leading_genes),
  Cross_contrast_membership = identical(cross_counts, shared_row_counts),
  Glaucoma_reproduced = all(main_all_le_glaucoma_checks)
)
print(main_all_le_summary, row.names = FALSE)
cat("\nOriginal gene-FDR support among unique shared contributors per contrast\n")
print(main_all_le_support[main_all_le_support$Pathway_group == "All selected pathways",
  c("Study", "Comparison", "Shared_genes", "Both_gene_FDR005", "Flight_only_gene_FDR005",
    "Terrestrial_only_gene_FDR005", "Neither_gene_FDR005")], row.names = FALSE)
print(main_all_le_checks)
stopifnot(all(main_all_le_checks))

# Save complete contributors, per-pathway overlap and original gene support without reranking diseases.
main_all_le_method <- list(
  selection = "Every pathway with original FDR < 0.05 in both flight and each terrestrial contrast; all 12 contrasts retained",
  shared_gene = "The same Entrez gene belongs to both leading edges of the same pathway",
  direction = "Pathway NES Same/Opposite/Zero_NES is recorded separately from the actual gene logFC signs",
  gene_filter = "None; retain contributors with any original gene FDR and label Both/Flight only/Terrestrial only/Neither at .05",
  gene_statistics = "Original full-tested-family within-study logFC, signed score and gene FDR; no recalculation",
  counting = "Unique genes counted once within each contrast and pathway-direction group; groups/contrasts can overlap and are not additive independent evidence",
  cross_model_table = "One row per pathway plus Entrez ID shared in at least one contrast; cells indicate shared/not-shared leading edges or pathway not both-enriched",
  missing_overlap = "Zero shared leading-edge genes in an opposing pathway can reflect different contributing subsets; not absence of pathway-level opposition",
  exclusions = "A pathway not both-enriched or a gene outside a leading edge is not evidence of absent expression or no response",
  inference = "Descriptive contributor overlap and original FDR flags only; no joint FDR, overlap P value, causal ranking or new disease similarity measure",
  provenance = "Frozen comparison/enrichment/common-response/annotation/pathway records and exact reproduction of reviewed glaucoma2W shared genes",
  documentation = "https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/#running-a-leading-edge-analysis"
)
main_all_le_files_out <- file.path(main_analysis_dir, c(
  "Main_analysis_OSD758_leading_edge_overlap_full_panel.csv",
  "Main_analysis_OSD758_leading_edge_membership_full_panel.csv",
  "Main_analysis_OSD758_shared_leading_edge_genes_full_panel.csv",
  "Main_analysis_OSD758_leading_edge_summary_full_panel.csv",
  "Main_analysis_OSD758_leading_edge_gene_FDR_support_full_panel.csv",
  "Main_analysis_OSD758_shared_leading_genes_by_pathway_full_panel.csv",
  "Main_analysis_leading_edge_full_panel_checkpoint_2026-11-07.rds",
  "Main_analysis_OSD758_leading_edge_full_panel_method_2026-11-07.txt"
))
tables <- list(main_all_le_overlap, main_all_le_membership, main_all_le_genes,
  main_all_le_summary, main_all_le_support, main_all_le_cross)
for (i in seq_along(tables)) utils::write.csv(tables[[i]], main_all_le_files_out[i], row.names = FALSE)
saveRDS(list(pathway_overlap = main_all_le_overlap, membership = main_all_le_membership,
  shared_genes = main_all_le_genes, summary = main_all_le_summary, gene_FDR_support = main_all_le_support,
  shared_genes_by_pathway = main_all_le_cross, focus_pathways = main_all_le_focus,
  checks = main_all_le_checks, method = main_all_le_method, reference = main_all_le_reference,
  contrasts = main_all_le_terrestrial, input_files = main_all_le_files_in, input_md5 = main_all_le_input_md5,
  annotation_versions = main_all_le_annotation_context$annotation_versions), main_all_le_files_out[7])
utils::capture.output(print(main_all_le_method), file = main_all_le_files_out[8])
print(file.exists(main_all_le_files_out))
