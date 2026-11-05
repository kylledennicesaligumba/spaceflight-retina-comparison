# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# RP: full signed gene ranking for the sex-adjusted P23H-versus-WT comparison.
# Run in the current within-study R session after reviewing the result figures.

# Confirm that the complete results and RP output folder are available.
stopifnot(exists("rp_results"), dir.exists(rp_within_dir))

# Report the original F values before applying the same numerical guard as glaucoma.
rp_rank_F_check <- c(
  Minimum_F = min(rp_results$F, na.rm = TRUE),
  Negative_F = sum(rp_results$F < 0, na.rm = TRUE),
  Nonfinite_F = sum(!is.finite(rp_results$F))
)
print(rp_rank_F_check)

# Allow only tiny negative F values within R's default numerical comparison tolerance.
# This is a numerical guard, not a differential-expression significance threshold.
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/all.equal.html
rp_rank_tolerance <- sqrt(.Machine$double.eps)
stopifnot(
  all(is.finite(rp_results$F)),
  all(rp_results$F >= -rp_rank_tolerance)
)

# Add a signed score to every tested gene, retaining the original result columns.
# logFC supplies direction; sqrt(F) supplies magnitude. Tiny negative inputs become zero.
rp_ranked <- transform(
  rp_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)
stopifnot(all(is.finite(rp_ranked$rank_score)))

# Sort from the largest positive score to the most negative score.
rp_ranked <- rp_ranked[
  order(rp_ranked$rank_score, decreasing = TRUE),
]

# Save the full ranked table, including gene identifiers, annotations and statistics.
rp_ranked_file <- file.path(
  rp_within_dir,
  "RP_GSE156533_P23H_vs_WT_ranked_all_genes.csv"
)
write.csv(rp_ranked, file = rp_ranked_file, row.names = FALSE)

# Confirm the number of ranked genes and whether the file was saved.
print(nrow(rp_ranked))
print(file.exists(rp_ranked_file))

# Preview the first five and last five genes using the same columns as glaucoma.
print(head(rp_ranked[, c("gene_name", "logFC", "F", "rank_score", "FDR")], 5))
print(tail(rp_ranked[, c("gene_name", "logFC", "F", "rank_score", "FDR")], 5))
