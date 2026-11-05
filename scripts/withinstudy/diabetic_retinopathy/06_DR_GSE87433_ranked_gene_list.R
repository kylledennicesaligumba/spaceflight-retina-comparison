# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# DR: full signed gene ranking for Diabetic versus Normal retina.
# Keep every tested gene, using the existing limma moderated t statistic as the score.
# The score already contains direction, so no signed-square-root-F conversion is needed.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

# Confirm that the complete results, primary expression matrix and output folder are available.
stopifnot(
  exists("dr_results"), exists("dr_primary_expr"), dir.exists(dr_within_dir)
)

# Report the number of tested genes and any scores that cannot be used for ranking.
dr_rank_check <- c(
  Tested = nrow(dr_results),
  Nonfinite_scores = sum(!is.finite(dr_results$t))
)
print(dr_rank_check)
stopifnot(
  nrow(dr_results) == nrow(dr_primary_expr),
  all(is.finite(dr_results$t)),
  all(is.finite(dr_results$logFC)),
  anyDuplicated(dr_results$ID) == 0L
)

# Add the signed moderated t score to a separate copy of the complete results.
# Positive scores point toward higher Diabetic expression; negative scores toward lower.
dr_ranked <- transform(dr_results, rank_score = t)

# Order from the largest positive score to the most negative score.
# Keep the original annotation, effect estimates, P values and FDR values alongside the score.
dr_ranked <- dr_ranked[
  order(dr_ranked$rank_score, decreasing = TRUE), , drop = FALSE
]

# Save the complete annotated ranking in the DR within-study folder.
dr_ranked_file <- file.path(
  dr_within_dir, "DR_GSE87433_Diabetic_vs_Normal_ranked_all_genes.csv"
)
write.csv(dr_ranked, dr_ranked_file, row.names = FALSE)

# Confirm the full ranked gene count and that the CSV exists.
print(nrow(dr_ranked))
print(file.exists(dr_ranked_file))

# Preview the first five and last five genes in the signed ranking.
# FDR remains the DEG criterion; a position in this full ranking does not itself define a DEG.
print(head(
  dr_ranked[, c("Symbol", "logFC", "t", "rank_score", "adj.P.Val")], n = 5
))
print(tail(
  dr_ranked[, c("Symbol", "logFC", "t", "rank_score", "adj.P.Val")], n = 5
))
