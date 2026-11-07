# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Assess conditional ranking stability by deleting related pathway neighborhoods together.
# These deterministic sensitivity ranges are NOT biological confidence intervals.

# Measure overlap relative to the smaller matched pathway, retaining nested relationships.
main_stability_overlap <- function(pathways) {
  stopifnot(length(pathways) >= 4L, !is.null(names(pathways)), anyDuplicated(names(pathways)) == 0L,
    all(lengths(pathways) > 0L), all(vapply(pathways, function(x) !anyNA(x) && anyDuplicated(x) == 0L, logical(1))),
    requireNamespace("Matrix", quietly = TRUE))
  genes <- sort(unique(unlist(pathways, use.names = FALSE)), method = "radix")
  incidence <- Matrix::sparseMatrix(i = rep(seq_along(pathways), lengths(pathways)),
    j = match(unlist(pathways, use.names = FALSE), genes), x = 1,
    dims = c(length(pathways), length(genes)))
  shared <- as.matrix(Matrix::tcrossprod(incidence))
  coefficient <- shared / outer(lengths(pathways), lengths(pathways), pmin)
  dimnames(coefficient) <- list(names(pathways), names(pathways))
  stopifnot(all(is.finite(coefficient)), all(coefficient >= 0 & coefficient <= 1),
    all(diag(coefficient) == 1), max(abs(coefficient - t(coefficient))) < 1e-12)
  coefficient
}

# Recalculate all reference correlations on the identical retained pathways in a deletion scenario.
main_stability_rho <- function(nes, reference, terrestrial_keys) {
  stopifnot(nrow(nes) >= 3L, all(is.finite(nes)),
    all(apply(nes[, c(reference, terrestrial_keys), drop = FALSE], 2, function(x) length(unique(x)) > 1L)))
  stats::setNames(vapply(terrestrial_keys, function(key) unname(stats::cor(nes[, reference],
    nes[, key], method = "spearman")), numeric(1)), terrestrial_keys)
}

# Draw deletion-scenario ranges and baseline points, explicitly separating these from confidence intervals.
main_stability_plot <- function(summary, labels, file, type) {
  if (type == "png") grDevices::png(file, width = 3100, height = 1900, res = 180, pointsize = 15)
  else grDevices::pdf(file, width = 17.22, height = 10.56, pointsize = 15)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- summary[order(summary$Baseline_rank), , drop = FALSE]
  y <- rev(seq_len(nrow(x)))
  graphics::par(oma = c(5.1, 0, 4.8, 0))
  graphics::layout(matrix(1:3, nrow = 1), widths = c(1.2, 1.2, .85))
  graphics::par(mar = c(4, .5, 1.5, .3))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(.5, nrow(x) + .5), xaxs = "i", yaxs = "i")
  graphics::text(.98, y, labels = unname(labels[paste(x$Study, x$Comparison, sep = "__")]),
    adj = c(1, .5), cex = .78)
  graphics::par(mar = c(4, .2, 1.5, 1.5))
  graphics::plot(NA, xlim = c(-1, 1), ylim = c(.5, nrow(x) + .5), xaxs = "i", yaxs = "i",
    yaxt = "n", xlab = "Pathway-profile Spearman rho", ylab = "", bty = "l")
  graphics::abline(v = 0, lty = 2, col = "#999999")
  graphics::segments(x$Minimum_rho, y, x$Maximum_rho, y, lwd = 3, col = "#56B4E9")
  graphics::points(x$Baseline_rho, y, pch = 16, cex = .95)
  graphics::par(mar = c(4, .6, 1.5, 1.8))
  graphics::plot(NA, xlim = c(.5, nrow(x) + .5), ylim = c(.5, nrow(x) + .5),
    xaxs = "i", yaxs = "i", yaxt = "n", xaxt = "n", xlab = "Observed rank (1 = highest)", ylab = "", bty = "l")
  graphics::axis(1, at = seq_len(nrow(x)), cex.axis = .8)
  graphics::segments(x$Best_rank, y, x$Worst_rank, y, lwd = 3, col = "#E69F00")
  graphics::points(x$Baseline_rank, y, pch = 16, cex = .95)
  graphics::mtext("Pathway ranking stability after related-program deletions", side = 3, outer = TRUE,
    line = 2.8, cex = 1.15, font = 2)
  graphics::mtext("Black points: original 799-pathway results | Lines: minimum-to-maximum deletion-scenario values",
    side = 3, outer = TRUE, line = 1.1, cex = .82)
  graphics::mtext("Neighborhood: shared genes / size of smaller pathway >= 0.5; remove the anchor and all its direct neighbors.",
    side = 1, outer = TRUE, line = 1.1, cex = .78)
  graphics::mtext("Ranges shown retain at least half of the pathways; all deletion scenarios are saved. Duplicate neighborhoods count once.",
    side = 1, outer = TRUE, line = 2.6, cex = .76)
  graphics::mtext("Sensitivity ranges are not confidence intervals, significance tests or probabilities of being the closest disease.",
    side = 1, outer = TRUE, line = 4.1, cex = .78, font = 2)
}

# Restore the original complete pathway comparison and frozen matched gene memberships.
if (!exists("main_analysis_dir")) main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_stab_inputs <- file.path(main_analysis_dir, c("Main_analysis_pathway_comparison_checkpoint_2026-11-06.rds",
  "Main_analysis_Reactome_common_pathways_checkpoint_2026-11-06.rds",
  "Main_analysis_gene_coverage_sensitivity_checkpoint_2026-11-07.rds"))
stopifnot(all(file.exists(main_stab_inputs)))
main_stab_context <- readRDS(main_stab_inputs[1])
main_stab_path_context <- readRDS(main_stab_inputs[2])
main_stab_coverage <- readRDS(main_stab_inputs[3])
main_stab_nes <- main_stab_context$NES
main_stab_path_ids <- rownames(main_stab_nes)
main_stab_reference <- "OSD758__FLT_0G_vs_HGC"
main_stab_info <- main_stab_context$contrasts[main_stab_context$contrasts$Role == "Terrestrial", , drop = FALSE]
main_stab_keys <- paste(main_stab_info$Study, main_stab_info$Comparison, sep = "__")
main_stab_paths <- main_stab_path_context$pathways[main_stab_path_ids]
stopifnot(all(main_stab_context$checks), all(main_stab_path_context$checks), all(main_stab_coverage$checks),
  length(main_stab_keys) == 12L, setequal(colnames(main_stab_nes), c(main_stab_reference, main_stab_keys)),
  setequal(main_stab_path_ids, names(main_stab_path_context$pathways)),
  identical(main_stab_path_ids, names(main_stab_paths)),
  identical(unname(tools::md5sum(main_stab_inputs[2])),
    unname(main_stab_context$original_input_md5[main_stab_inputs[2]])),
  all(main_stab_keys %in% names(main_stab_coverage$labels)))

# Record this project-specific sensitivity protocol before examining new deletion outcomes.
main_stab_method <- list(
  status = "Exploratory post-primary robustness follow-up; numerical rules fixed before new deletion outcomes",
  purpose = "Check whether the observed pathway ranking depends on particular groups of overlapping gene sets",
  neighborhoods = "For each retained pathway, delete it and all direct neighbors with intersection size / smaller matched-set size >= 0.5",
  no_transitive_expansion = "Direct neighborhoods only; no connected-component closure or outcome-based manual grouping",
  duplicate_handling = "Identical removed-pathway sets evaluated once; all anchor identities retained",
  scope = "All retained pathways, regardless of FDR; same retained set for all13 contrast columns within each scenario",
  score = "Recompute descriptive Spearman correlations of the saved NES values after deletion; do not rerun enrichment",
  minimum = "At least3 retained pathways and nonconstant profiles are required to estimate correlations",
  principal_display = "Summaries/figure use valid deletions retaining at least half of the original pathways; heavier deletions remain separately saved",
  inference = "Finite deterministic sensitivity ranges; not biological confidence intervals, P values, independence-based bootstrap or probabilities of a winning disease",
  counts = "Top counts refer to distinct predefined deletion scenarios; neighborhoods are not independent trials and top counts do not have a sampling probability interpretation",
  biological_uncertainty = "Not quantified: original biological units would require design-aware resampling/refitting; glaucoma has one pooled library per sex/group, shared controls and nominal residual df3",
  decision = "Keep conclusions descriptive as highest observed resemblance, not a statistically demonstrated closest model",
  literature_component = "Merico2010 supports overlap-coefficient representation of related sets, not this project's deletion protocol or confidence intervals",
  references = c("https://pmc.ncbi.nlm.nih.gov/articles/PMC2981572/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/", "https://pubmed.ncbi.nlm.nih.gov/17303618/")
)
main_stab_files <- file.path(main_analysis_dir, c("Main_analysis_pathway_rank_stability_summary.csv",
  "Main_analysis_pathway_rank_stability_scenarios.csv", "Main_analysis_pathway_rank_stability_profiles.csv",
  "Main_analysis_pathway_rank_stability_checkpoint_2026-11-07.rds",
  "Main_analysis_pathway_rank_stability_method_2026-11-07.txt",
  "Main_analysis_OSD758_pathway_rank_stability.png", "Main_analysis_OSD758_pathway_rank_stability.pdf"))
utils::capture.output(print(main_stab_method), file = main_stab_files[5])
main_stab_baseline <- main_stability_rho(main_stab_nes, main_stab_reference, main_stab_keys)
main_stab_original_index <- match(main_stab_keys, paste(main_stab_context$similarity$Study,
  main_stab_context$similarity$Comparison, sep = "__"))
stopifnot(!anyNA(main_stab_original_index),
  max(abs(main_stab_baseline - main_stab_context$similarity$Spearman_rho[main_stab_original_index])) < 1e-12)
main_stab_overlap <- main_stability_overlap(main_stab_paths)
main_stab_anchor_removals <- lapply(seq_along(main_stab_path_ids), function(i) which(main_stab_overlap[i, ] >= .5))
main_stab_signatures <- vapply(main_stab_anchor_removals, function(i) paste(i, collapse = ","), character(1))
main_stab_unique_signatures <- unique(main_stab_signatures)
main_stab_first_anchor <- match(main_stab_unique_signatures, main_stab_signatures)
main_stab_removals <- main_stab_anchor_removals[main_stab_first_anchor]
main_stab_scenarios <- main_stab_profiles <- vector("list", length(main_stab_removals))
for (s in seq_along(main_stab_removals)) {
  excluded <- main_stab_removals[[s]]
  kept <- setdiff(seq_along(main_stab_path_ids), excluded)
  anchors <- which(main_stab_signatures == main_stab_unique_signatures[s])
  estimable <- length(kept) >= 3L && all(apply(main_stab_nes[kept, , drop = FALSE], 2,
    function(x) length(unique(x)) > 1L))
  rho <- if (estimable) main_stability_rho(main_stab_nes[kept, , drop = FALSE], main_stab_reference, main_stab_keys) else
    stats::setNames(rep(NA_real_, length(main_stab_keys)), main_stab_keys)
  ranks <- if (estimable) rank(-rho, ties.method = "min") else rep(NA_real_, length(rho))
  main_stab_scenarios[[s]] <- data.frame(Scenario = s,
    Anchor_pathway_IDs = paste(main_stab_path_ids[anchors], collapse = ";"),
    Anchor_pathway_names = paste(main_stab_context$pathway_key$Pathway_name[match(main_stab_path_ids[anchors],
      main_stab_context$pathway_key$Pathway_ID)], collapse = ";"),
    Removed_pathways = length(excluded), Retained_pathways = length(kept),
    Half_or_more_retained = length(kept) >= ceiling(nrow(main_stab_nes) / 2),
    Status = if (estimable) "Estimated" else "Not estimable",
    Highest_contrast = if (estimable) paste(main_stab_keys[ranks == 1], collapse = ";") else NA_character_,
    stringsAsFactors = FALSE)
  main_stab_profiles[[s]] <- data.frame(Scenario = s, Study = main_stab_info$Study,
    Comparison = main_stab_info$Comparison, Retained_pathways = length(kept),
    Spearman_rho = unname(rho), Observed_rank = unname(ranks), stringsAsFactors = FALSE)
}
main_stab_scenario_table <- do.call(rbind, main_stab_scenarios)
main_stab_profile_table <- do.call(rbind, main_stab_profiles)
main_stab_selected <- main_stab_scenario_table$Scenario[
  main_stab_scenario_table$Status == "Estimated" & main_stab_scenario_table$Half_or_more_retained]
stopifnot(length(main_stab_selected) > 0L)
main_stab_summary <- do.call(rbind, lapply(seq_along(main_stab_keys), function(j) {
  x <- main_stab_profile_table[main_stab_profile_table$Scenario %in% main_stab_selected &
    main_stab_profile_table$Study == main_stab_info$Study[j] &
    main_stab_profile_table$Comparison == main_stab_info$Comparison[j], ]
  data.frame(Study = main_stab_info$Study[j], Comparison = main_stab_info$Comparison[j],
    Baseline_rho = unname(main_stab_baseline[j]), Baseline_rank = unname(rank(-main_stab_baseline, ties.method = "min")[j]),
    Valid_half_retained_scenarios = nrow(x), Minimum_rho = min(x$Spearman_rho), Maximum_rho = max(x$Spearman_rho),
    Best_rank = min(x$Observed_rank), Worst_rank = max(x$Observed_rank),
    Top_rank_scenarios = sum(x$Observed_rank == 1), stringsAsFactors = FALSE)
}))
rownames(main_stab_summary) <- NULL
main_stab_checks <- c(Primary_reproduced = max(abs(main_stab_baseline -
    main_stab_context$similarity$Spearman_rho[main_stab_original_index])) < 1e-12,
  Frozen_memberships = identical(main_stab_path_ids, names(main_stab_paths)),
  Neighborhoods_valid = all(vapply(seq_along(main_stab_path_ids), function(i)
    i %in% main_stab_anchor_removals[[i]] && all(main_stab_overlap[i, main_stab_anchor_removals[[i]]] >= .5), logical(1))),
  Unique_deletions = anyDuplicated(vapply(main_stab_removals, function(i) paste(i, collapse = ","), character(1))) == 0L,
  All_contrasts = nrow(main_stab_profile_table) == length(main_stab_removals) * 12L,
  Finite_selected = all(is.finite(main_stab_profile_table$Spearman_rho[main_stab_profile_table$Scenario %in% main_stab_selected])),
  Rank_accounting = all(main_stab_summary$Valid_half_retained_scenarios == length(main_stab_selected)) &&
    all(main_stab_summary$Best_rank >= 1 & main_stab_summary$Worst_rank <= 12))
print(main_stab_summary[order(main_stab_summary$Baseline_rank), ], digits = 6, row.names = FALSE)
print(c(Anchor_pathways = length(main_stab_path_ids), Unique_deletion_scenarios = length(main_stab_removals),
  Estimated_scenarios = sum(main_stab_scenario_table$Status == "Estimated"), Principal_scenarios = length(main_stab_selected)))
print(main_stab_checks)
stopifnot(all(main_stab_checks))
utils::write.csv(main_stab_summary, main_stab_files[1], row.names = FALSE)
utils::write.csv(main_stab_scenario_table, main_stab_files[2], row.names = FALSE)
utils::write.csv(main_stab_profile_table, main_stab_files[3], row.names = FALSE)
saveRDS(list(summary = main_stab_summary, scenarios = main_stab_scenario_table, profiles = main_stab_profile_table,
  removed_pathway_indices = main_stab_removals, pathway_ids = main_stab_path_ids,
  overlap_coefficient = main_stab_overlap, anchor_scenario = match(main_stab_signatures, main_stab_unique_signatures),
  principal_scenarios = main_stab_selected, checks = main_stab_checks, method = main_stab_method,
  baseline_rho = main_stab_baseline, input_files = main_stab_inputs, input_md5 = tools::md5sum(main_stab_inputs)), main_stab_files[4])
main_stability_plot(main_stab_summary, main_stab_coverage$labels, main_stab_files[6], "png")
main_stability_plot(main_stab_summary, main_stab_coverage$labels, main_stab_files[7], "pdf")
print(file.exists(main_stab_files))
