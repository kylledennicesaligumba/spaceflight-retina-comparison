# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Independently verify the final checkpoint, representative selection, statistics and plotted values.
# Run in a fresh R process after the full analysis runner; this does not refit the model.
amd_check_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "amd")
amd_verified <- new.env()
load(file.path(amd_check_dir, "AMD_GSE38671_gene_level_withinstudy_final_checkpoint_2026-11-05.RData"), envir = amd_verified)
amd_checks <- character()
amd_check <- function(condition, label) {
  if (!isTRUE(condition)) stop(paste("FAILED:", label), call. = FALSE)
  amd_checks <<- c(amd_checks, paste("PASS:", label))
}
amd_check(identical(dim(amd_verified$amd_gene_expr), c(19825L, 16L)), "19,825 gene rows and 16 arrays")
amd_check(identical(colnames(amd_verified$amd_gene_expr), amd_verified$amd_within_metadata$geo_accession), "Sample order matches metadata")
amd_check(anyDuplicated(rownames(amd_verified$amd_gene_expr)) == 0L, "Gene IDs unique")
amd_check(identical(rownames(amd_verified$amd_gene_expr), amd_verified$amd_gene_annotation$entrez_id), "Gene annotation aligned")
amd_check(qr(amd_verified$amd_design)$rank == 4L && all(amd_verified$amd_gene_fit$df.residual == 12L), "Full-rank design and 12 residual degrees of freedom")

# Compare selected expression directly with the untouched precheck, without reprocessing it.
amd_original <- new.env()
load(amd_verified$amd_precheck_file, envir = amd_original)
amd_rep <- amd_verified$amd_gene_representatives
amd_check(identical(unname(amd_verified$amd_gene_expr), unname(amd_original$amd_expr_log2[amd_rep$ID, , drop = FALSE])),
          "Selected expression values exactly match precheck")
amd_check(identical(unname(tools::md5sum(amd_verified$amd_precheck_file)), unname(amd_verified$amd_input_checksum)), "Input checkpoint checksum unchanged")
amd_selection <- amd_verified$amd_feature_selection
amd_check(isTRUE(all.equal(amd_selection$mean_log2_expression,
  unname(rowMeans(amd_original$amd_expr_log2[amd_selection$ID, , drop = FALSE])), tolerance = 1e-12)), "All candidate means match original values")
amd_maxima <- tapply(amd_selection$mean_log2_expression, amd_selection$entrez_id, max)
amd_check(all(amd_rep$mean_log2_expression == amd_maxima[amd_rep$entrez_id]), "Maximum mean selected for every gene")
amd_tie_choices <- vapply(split(amd_selection, amd_selection$entrez_id), function(x) {
  sort(x$ID[x$mean_log2_expression == max(x$mean_log2_expression)], method = "radix")[1]
}, character(1))
amd_check(all(amd_rep$ID == amd_tie_choices[amd_rep$entrez_id]), "Deterministic feature-ID tie rule")
amd_check(nrow(amd_selection) == 21016L && sum(amd_selection$selected) == 19825L &&
          all(table(amd_selection$entrez_id[amd_selection$selected]) == 1L), "One selected representative per eligible gene")
amd_audit <- amd_verified$amd_mapping_audit
amd_check(sum(amd_audit$eligible) == 21016L && all(with(amd_audit[amd_audit$eligible, ],
  recorded_gene_IDs == 1 & !truncated_annotation & !missing_gene_ID_assignment)), "Conservative annotation eligibility")

# Check effect estimates against direct differences of the four observed group means.
amd_group_means <- sapply(levels(amd_verified$amd_within_metadata$group), function(g) {
  rowMeans(amd_verified$amd_gene_expr[, amd_verified$amd_within_metadata$group == g, drop = FALSE])
})
amd_expected_effects <- cbind(
  Young_CfhKO_vs_WT = amd_group_means[, "Young_CfhKO"] - amd_group_means[, "Young_WT"],
  Aged_CfhKO_vs_WT = amd_group_means[, "Aged_CfhKO"] - amd_group_means[, "Aged_WT"])
amd_check(max(abs(amd_expected_effects - amd_verified$amd_gene_contrast_fit$coefficients)) < 1e-10, "Age-matched contrast effect estimates")
for (comparison in names(amd_verified$amd_gene_results)) {
  x <- amd_verified$amd_gene_results[[comparison]]
  hits <- amd_verified$amd_gene_degs[[comparison]]
  ranked <- amd_verified$amd_gene_ranked[[comparison]]
  amd_check(nrow(x) == 19825L && all(is.finite(x$t)) && all(is.finite(x$P.Value)), paste(comparison, "complete finite results"))
  amd_check(max(abs(x$adj.P.Val - stats::p.adjust(x$P.Value, "BH"))) < 1e-12, paste(comparison, "BH FDR across all tested genes"))
  amd_check(identical(hits$entrez_id, x$entrez_id[x$adj.P.Val < 0.05]), paste(comparison, "DEG membership"))
  amd_check(nrow(ranked) == 19825L && all(diff(ranked$rank_score) <= 0) &&
    all(ranked$rank_score == ranked$t) && setequal(ranked$entrez_id, x$entrez_id), paste(comparison, "full signed moderated-t ranking"))
  csv <- read.csv(file.path(amd_check_dir, paste0("AMD_GSE38671_gene_level_", comparison, "_all_genes.csv")),
                  colClasses = c(entrez_id = "character", feature_ID = "character"))
  amd_check(identical(csv$entrez_id, x$entrez_id) && max(abs(csv$logFC - x$logFC)) < 1e-12 &&
    max(abs(csv$adj.P.Val - x$adj.P.Val)) < 1e-12, paste(comparison, "saved CSV matches checkpoint"))
  heat <- amd_verified$amd_heatmaps[[comparison]]
  amd_check(setequal(heat$genes$entrez_id, hits$entrez_id), paste(comparison, "heatmap includes every DEG"))
  amd_check(isTRUE(all.equal(heat$expression, amd_verified$amd_gene_expr[heat$genes$entrez_id,
    amd_verified$amd_heat_samples$geo_accession, drop = FALSE])), paste(comparison, "heatmap expression matches input"))
  amd_check(isTRUE(all.equal(heat$zscores, t(scale(t(heat$expression))))), paste(comparison, "heatmap row z-scores"))
}

# Confirm that every dot is the exact measured value for its recorded gene and sample.
amd_dots <- amd_verified$amd_dot_data
amd_values <- amd_verified$amd_gene_expr[cbind(match(amd_dots$entrez_id, rownames(amd_verified$amd_gene_expr)),
  match(amd_dots$sample_id, colnames(amd_verified$amd_gene_expr)))]
amd_check(nrow(amd_dots) == 176L && length(unique(amd_dots$entrez_id)) == 11L &&
  anyDuplicated(amd_dots[, c("entrez_id", "sample_id")]) == 0L &&
  all(table(amd_dots$entrez_id, amd_dots$group) == 4L) && all(amd_dots$Expression == amd_values), "11 genes x 16 exact individual-array dot values")
amd_check(all(table(amd_verified$amd_md_data$Comparison) == 19825L), "All genes shown in both MD panels")
amd_check(length(amd_verified$amd_figure_files) == 5L && all(file.exists(amd_verified$amd_figure_files)) &&
  all(file.info(amd_verified$amd_figure_files)$size > 0), "Five nonempty saved PNG figures")
amd_legacy <- read.csv(file.path(amd_check_dir, "AMD_GSE38671_significant_feature_summary.csv"))
amd_check(identical(as.integer(amd_legacy$Significant_features), c(5L, 51L)), "Earlier feature-level summary preserved")

# Save the verification record beside the results for later review.
amd_verification_file <- file.path(amd_check_dir, "AMD_GSE38671_gene_level_verification_2026-11-05.txt")
writeLines(c("AMD GSE38671 final gene-level verification — 5 November 2026, Asia/Manila", amd_checks,
  "Visual inspection of all five figures is recorded separately in the analysis notes."), amd_verification_file)
cat(paste(amd_checks, collapse = "\n"), "\n")
