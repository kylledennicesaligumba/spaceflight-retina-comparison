# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# LHON-linked model: full signed gene rankings for all three genotype-versus-WT tests.
# Run in the current within-study session after reviewing the expression figures.
# Reuse the agreed signed-sqrt-F score from glaucoma and RP, retaining every tested gene.

# Confirm that the complete results and LHON output folder are available.
stopifnot(
  exists("lhon_YARS2_results"), exists("lhon_COI_results"),
  exists("lhon_YARS2_COI_results"), dir.exists(lhon_within_dir)
)

# Report the original F statistics before applying the familiar numerical guard.
lhon_rank_F_check <- rbind(
  YARS2 = c(
    Minimum_F = min(lhon_YARS2_results$F, na.rm = TRUE),
    Negative_F = sum(lhon_YARS2_results$F < 0, na.rm = TRUE),
    Nonfinite_F = sum(!is.finite(lhon_YARS2_results$F))
  ),
  COI = c(
    Minimum_F = min(lhon_COI_results$F, na.rm = TRUE),
    Negative_F = sum(lhon_COI_results$F < 0, na.rm = TRUE),
    Nonfinite_F = sum(!is.finite(lhon_COI_results$F))
  ),
  YARS2_COI = c(
    Minimum_F = min(lhon_YARS2_COI_results$F, na.rm = TRUE),
    Negative_F = sum(lhon_YARS2_COI_results$F < 0, na.rm = TRUE),
    Nonfinite_F = sum(!is.finite(lhon_YARS2_COI_results$F))
  )
)
print(lhon_rank_F_check)

# Allow only tiny negative F values within R's default numerical comparison tolerance.
# This is a numerical guard, not a differential-expression significance threshold.
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/all.equal.html
lhon_rank_tolerance <- sqrt(.Machine$double.eps)
stopifnot(
  all(is.finite(lhon_YARS2_results$F)),
  all(lhon_YARS2_results$F >= -lhon_rank_tolerance),
  all(is.finite(lhon_COI_results$F)),
  all(lhon_COI_results$F >= -lhon_rank_tolerance),
  all(is.finite(lhon_YARS2_COI_results$F)),
  all(lhon_YARS2_COI_results$F >= -lhon_rank_tolerance)
)

# Add signed scores to separate copies of the three complete result tables.
# logFC supplies direction; sqrt(F) supplies the test-statistic magnitude.
# Only within-tolerance negative inputs become zero for ranking; original F values remain.
lhon_YARS2_ranked <- transform(
  lhon_YARS2_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)
lhon_COI_ranked <- transform(
  lhon_COI_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)
lhon_YARS2_COI_ranked <- transform(
  lhon_YARS2_COI_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)

# Check that every score is finite before sorting or exporting.
stopifnot(
  all(is.finite(lhon_YARS2_ranked$rank_score)),
  all(is.finite(lhon_COI_ranked$rank_score)),
  all(is.finite(lhon_YARS2_COI_ranked$rank_score))
)

# Sort each table from the largest positive score to the most negative score.
lhon_YARS2_ranked <- lhon_YARS2_ranked[
  order(lhon_YARS2_ranked$rank_score, decreasing = TRUE),
]
lhon_COI_ranked <- lhon_COI_ranked[
  order(lhon_COI_ranked$rank_score, decreasing = TRUE),
]
lhon_YARS2_COI_ranked <- lhon_YARS2_COI_ranked[
  order(lhon_YARS2_COI_ranked$rank_score, decreasing = TRUE),
]

# Set one ranked-table filename per genotype comparison inside the LHON output folder.
lhon_ranked_files <- file.path(
  lhon_within_dir,
  c(
    "LHON_GSE304707_YARS2_vs_WT_ranked_all_genes.csv",
    "LHON_GSE304707_COI_vs_WT_ranked_all_genes.csv",
    "LHON_GSE304707_YARS2_COI_vs_WT_ranked_all_genes.csv"
  )
)

# Save gene names and every original result column alongside the ranking score.
write.csv(lhon_YARS2_ranked, file = lhon_ranked_files[1], row.names = FALSE)
write.csv(lhon_COI_ranked, file = lhon_ranked_files[2], row.names = FALSE)
write.csv(lhon_YARS2_COI_ranked, file = lhon_ranked_files[3], row.names = FALSE)

# Confirm the number of ranked genes per comparison and whether all files were saved.
print(c(
  YARS2 = nrow(lhon_YARS2_ranked),
  COI = nrow(lhon_COI_ranked),
  YARS2_COI = nrow(lhon_YARS2_COI_ranked)
))
print(file.exists(lhon_ranked_files))

# Preview the first five and last five genes in the combined-genotype ranked list.
# All three lists are saved; this console preview shows YARS2/COI versus WT only.
print(head(
  lhon_YARS2_COI_ranked[, c("gene_name", "logFC", "F", "rank_score", "FDR")], 5
))
print(tail(
  lhon_YARS2_COI_ranked[, c("gene_name", "logFC", "F", "rank_score", "FDR")], 5
))
