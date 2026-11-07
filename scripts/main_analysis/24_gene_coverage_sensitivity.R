# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Check whether gene-level resemblance depends on the all-study common-gene restriction.
# Preserve the original signed scores, mapping exclusions and primary analysis.

# Extract eligible genes by their verified links to original ranked-table rows.
main_coverage_extract_scores <- function(table, map, ids, source_column) {
  stopifnot(anyDuplicated(map$Entrez_ID) == 0L, anyDuplicated(ids) == 0L,
    all(c(source_column, "rank_score", "logFC") %in% names(table)))
  m <- match(ids, as.character(map$Entrez_ID))
  stopifnot(!anyNA(m))
  selected <- map[m, , drop = FALSE]
  rows <- selected$Source_row
  stopifnot(all(is.finite(rows)), all(rows == floor(rows)), all(rows >= 1L & rows <= nrow(table)),
    identical(as.character(table[[source_column]][rows]), as.character(selected$Source_ID)))
  scores <- table$rank_score[rows]
  changes <- table$logFC[rows]
  stopifnot(is.numeric(scores), all(is.finite(scores)), all(is.finite(changes)),
    all(scores == 0 | sign(scores) == sign(changes)), length(unique(scores)) > 1L)
  list(scores = scores, source_rows = rows, source_ids = as.character(selected$Source_ID))
}

# Draw both gene-coverage choices on the full correlation scale.
main_coverage_plot <- function(result, labels, file, type) {
  if (type == "png") grDevices::png(file, width = 2800, height = 1850, res = 180, pointsize = 15)
  else grDevices::pdf(file, width = 15.55, height = 10.28, pointsize = 15)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- result[order(result$Primary_rank), , drop = FALSE]
  y <- rev(seq_len(nrow(result)))
  keys <- paste(result$Study, result$Comparison, sep = "__")
  graphics::par(mar = c(6, 18.8, 5.8, 2))
  graphics::plot(NA, xlim = c(-1, 1), ylim = c(0.5, nrow(result) + 0.6),
    xlab = "Spearman rho of signed gene-response scores", ylab = "", yaxt = "n", bty = "l", xaxs = "i")
  graphics::abline(v = 0, lty = 2, col = "#777777")
  graphics::axis(2, at = y, labels = unname(labels[keys]), las = 1, tick = FALSE, cex.axis = .83)
  graphics::segments(result$Primary_rho, y + .12, result$Pairwise_rho, y - .12, col = "#999999", lwd = 1.4)
  graphics::points(result$Primary_rho, y + .12, pch = 16, cex = 1.1, col = "#0072B2")
  graphics::points(result$Pairwise_rho, y - .12, pch = 17, cex = 1.1, col = "#D55E00")
  graphics::title("Gene-coverage sensitivity of resemblance to OSD-758", line = 3, cex.main = 1.15)
  graphics::mtext("Primary: same genes in all 13 contrasts | Sensitivity: all eligible genes shared by each flight-model pair",
    side = 3, line = 1.5, cex = .76)
  graphics::legend("topleft", legend = c("All-study common genes", "Available pairwise genes"),
    pch = c(16, 17), col = c("#0072B2", "#D55E00"), bty = "n", cex = .78)
  graphics::mtext("Original scores and ID eligibility retained; all genes included regardless of FDR.",
    side = 1, line = 3.8, cex = .78)
  graphics::mtext("This is a coverage sensitivity analysis, not a replacement primary ranking or a confidence interval.",
    side = 1, line = 5, cex = .75)
}

# Restore the reviewed mapping, common responses and primary gene correlations.
if (!exists("main_analysis_dir")) main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_cov_input <- file.path(main_analysis_dir, c("Main_analysis_gene_coverage_checkpoint_2026-11-06.rds",
  "Main_analysis_common_responses_checkpoint_2026-11-06.rds", "Main_analysis_gene_Spearman_checkpoint_2026-11-06.rds"))
stopifnot(all(file.exists(main_cov_input)))
main_cov_mapping <- readRDS(main_cov_input[1])
main_cov_response <- readRDS(main_cov_input[2])
main_cov_primary <- readRDS(main_cov_input[3])
main_cov_manifest <- main_cov_mapping$inputs
main_cov_keys <- paste(main_cov_manifest$Study, main_cov_manifest$Comparison, sep = "__")
main_cov_reference <- "OSD758__FLT_0G_vs_HGC"
main_cov_terrestrial <- which(main_cov_manifest$Role == "Terrestrial")
main_cov_common_ids <- main_cov_response$common_ids
stopifnot(length(main_cov_terrestrial) == 12L, anyDuplicated(main_cov_keys) == 0L,
  sum(main_cov_manifest$Role == "Reference") == 1L,
  main_cov_manifest$Role[match(main_cov_reference, main_cov_keys)] == "Reference",
  all(main_cov_response$checks), identical(main_cov_common_ids, main_cov_mapping$candidate_common_ids),
  identical(colnames(main_cov_response$rank_scores), main_cov_keys),
  identical(rownames(main_cov_response$rank_scores), main_cov_common_ids),
  identical(unname(tools::md5sum(main_cov_input[2])), unname(main_cov_primary$response_checkpoint_md5)))
main_cov_source_md5 <- tools::md5sum(main_cov_manifest$File)
main_cov_integrity <- identical(main_cov_source_md5, main_cov_mapping$input_md5) &&
  identical(main_cov_source_md5, main_cov_response$input_md5)
stopifnot(main_cov_integrity)

# Fix this follow-up's rules before calculating any new pairwise correlations.
main_cov_method <- list(
  status = "Exploratory planned coverage sensitivity; rules fixed before these new results",
  measure = "stats::cor(method='spearman') with default average tie ranks",
  primary = "Same accepted all13 common genes in every contrast; primary results unchanged",
  sensitivity = "For each terrestrial contrast, use all eligible Entrez IDs also eligible in the fixed OSD reference",
  eligibility = "Exactly the frozen module09 row maps; no rescue of ambiguous/repeated/unmapped targets",
  score = "Original saved signed within-study scores, unchanged; no refit, normalization or FDR selection",
  comparison = "Reproduce primary scores/rho before expanding coverage; retain all12 separate contrasts",
  inference = "Changes assess dependence on coverage; different pairwise gene universes are not a fairer replacement ranking or biological confidence intervals",
  documentation = "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html"
)
main_cov_files <- file.path(main_analysis_dir, c("Main_analysis_OSD758_gene_coverage_sensitivity.csv",
  "Main_analysis_pairwise_gene_coverage_row_key.csv", "Main_analysis_gene_coverage_sensitivity_checkpoint_2026-11-07.rds",
  "Main_analysis_gene_coverage_sensitivity_method_2026-11-07.txt",
  "Main_analysis_OSD758_gene_coverage_sensitivity.png", "Main_analysis_OSD758_gene_coverage_sensitivity.pdf"))
utils::capture.output(print(main_cov_method), file = main_cov_files[4])
main_cov_tables <- lapply(main_cov_manifest$File, utils::read.csv, check.names = FALSE, stringsAsFactors = FALSE)
names(main_cov_tables) <- main_cov_keys
main_cov_source_columns <- vapply(main_cov_keys, function(key) {
  info <- main_cov_response$source_row_key
  field <- unique(info$Source_column[paste(info$Study, info$Comparison, sep = "__") == key])
  stopifnot(length(field) == 1L)
  field
}, character(1))
main_cov_reference_map <- main_cov_mapping$eligible_row_maps[[main_cov_reference]]
main_cov_reference_common <- main_coverage_extract_scores(main_cov_tables[[main_cov_reference]],
  main_cov_reference_map, main_cov_common_ids, main_cov_source_columns[[main_cov_reference]])
main_cov_rows <- main_cov_pair_tables <- main_cov_pair_scores <- vector("list", length(main_cov_terrestrial))
main_cov_score_checks <- main_cov_count_checks <- main_cov_subset_checks <- logical(length(main_cov_terrestrial))
for (j in seq_along(main_cov_terrestrial)) {
  i <- main_cov_terrestrial[j]
  key <- main_cov_keys[i]
  map <- main_cov_mapping$eligible_row_maps[[key]]
  ids <- sort(intersect(as.character(main_cov_reference_map$Entrez_ID), as.character(map$Entrez_ID)), method = "radix")
  main_cov_subset_checks[j] <- all(main_cov_common_ids %in% ids)
  a <- main_coverage_extract_scores(main_cov_tables[[main_cov_reference]], main_cov_reference_map,
    ids, main_cov_source_columns[[main_cov_reference]])
  b <- main_coverage_extract_scores(main_cov_tables[[key]], map, ids, main_cov_source_columns[[key]])
  primary_b <- main_coverage_extract_scores(main_cov_tables[[key]], map,
    main_cov_common_ids, main_cov_source_columns[[key]])
  main_cov_score_checks[j] <- identical(main_cov_reference_common$scores,
    unname(main_cov_response$rank_scores[, main_cov_reference])) &&
    identical(primary_b$scores, unname(main_cov_response$rank_scores[, key])) &&
    identical(a$scores[match(main_cov_common_ids, ids)], main_cov_reference_common$scores) &&
    identical(b$scores[match(main_cov_common_ids, ids)], primary_b$scores)
  primary_rho <- unname(stats::cor(main_cov_reference_common$scores, primary_b$scores, method = "spearman"))
  pairwise_rho <- unname(stats::cor(a$scores, b$scores, method = "spearman"))
  main_cov_count_checks[j] <- length(ids) == main_cov_mapping$coverage$Shared_with_OSD[
    match(key, paste(main_cov_mapping$coverage$Study, main_cov_mapping$coverage$Comparison, sep = "__"))]
  main_cov_rows[[j]] <- data.frame(Study = main_cov_manifest$Study[i], Comparison = main_cov_manifest$Comparison[i],
    Primary_genes = length(main_cov_common_ids), Pairwise_genes = length(ids), Added_genes = length(ids) - length(main_cov_common_ids),
    Primary_rho = primary_rho, Pairwise_rho = pairwise_rho, Rho_change = pairwise_rho - primary_rho)
  main_cov_pair_tables[[j]] <- data.frame(Study = main_cov_manifest$Study[i], Comparison = main_cov_manifest$Comparison[i],
    Entrez_ID = ids, OSD_source_row = a$source_rows, Terrestrial_source_row = b$source_rows,
    OSD_source_ID = a$source_ids, Terrestrial_source_ID = b$source_ids)
  main_cov_pair_scores[[j]] <- list(Entrez_ID = ids, OSD_rank_score = a$scores, Terrestrial_rank_score = b$scores)
}
main_cov_result <- do.call(rbind, main_cov_rows)
main_cov_result$Primary_rank <- rank(-main_cov_result$Primary_rho, ties.method = "min")
main_cov_result$Pairwise_rank <- rank(-main_cov_result$Pairwise_rho, ties.method = "min")
main_cov_result$Rank_change <- main_cov_result$Pairwise_rank - main_cov_result$Primary_rank
main_cov_primary_index <- match(paste(main_cov_result$Study, main_cov_result$Comparison, sep = "__"),
  paste(main_cov_primary$results$Study, main_cov_primary$results$Comparison, sep = "__"))
main_cov_checks <- c(Input_integrity = main_cov_integrity,
  All_contrasts = nrow(main_cov_result) == 12L && !anyNA(main_cov_primary_index),
  Preserved_scores = all(main_cov_score_checks), Common_subset = all(main_cov_subset_checks),
  Pairwise_counts = all(main_cov_count_checks),
  Primary_reproduced = max(abs(main_cov_result$Primary_rho - main_cov_primary$results$Spearman_rho[main_cov_primary_index])) < 1e-12,
  Finite_rho = all(is.finite(main_cov_result$Pairwise_rho) & abs(main_cov_result$Pairwise_rho) <= 1))
print(main_cov_result[order(main_cov_result$Primary_rank), ], digits = 6, row.names = FALSE)
print(main_cov_checks)
stopifnot(all(main_cov_checks))
main_cov_labels <- c(Glaucoma__3D_vs_Naive = "Glaucoma: 3 days vs naive",
  Glaucoma__2W_vs_Naive = "Glaucoma: 2 weeks vs naive", Glaucoma__6W_vs_Naive = "Glaucoma: 6 weeks vs naive",
  RP__P23H_vs_WT = "RP model: P23H vs WT", LHON__YARS2_vs_WT = "LHON-linked: YARS2 vs WT",
  LHON__COI_vs_WT = "LHON-linked: COI vs WT", LHON__YARS2_COI_vs_WT = "LHON-linked: YARS2 + COI vs WT",
  DR__Diabetic_vs_Normal = "Diabetic retinal injury: diabetic vs normal",
  AMD__Young_CfhKO_vs_WT = "AMD-related: young CfhKO vs WT", AMD__Aged_CfhKO_vs_WT = "AMD-related: aged CfhKO vs WT",
  Autoimmune_uveitis__EAU_vs_CON = "Autoimmune uveitis: EAU vs control", Light_damage__LD_vs_CTRL = "Light damage: LD vs control")
stopifnot(all(main_cov_keys[main_cov_terrestrial] %in% names(main_cov_labels)))
utils::write.csv(main_cov_result, main_cov_files[1], row.names = FALSE)
utils::write.csv(do.call(rbind, main_cov_pair_tables), main_cov_files[2], row.names = FALSE)
names(main_cov_pair_scores) <- main_cov_keys[main_cov_terrestrial]
saveRDS(list(results = main_cov_result, pairwise_scores = main_cov_pair_scores, checks = main_cov_checks,
  method = main_cov_method, labels = main_cov_labels, common_ids = main_cov_common_ids,
  input_files = main_cov_input, input_md5 = tools::md5sum(main_cov_input), source_input_md5 = main_cov_source_md5), main_cov_files[3])
main_coverage_plot(main_cov_result, main_cov_labels, main_cov_files[5], "png")
main_coverage_plot(main_cov_result, main_cov_labels, main_cov_files[6], "pdf")
print(file.exists(main_cov_files))
