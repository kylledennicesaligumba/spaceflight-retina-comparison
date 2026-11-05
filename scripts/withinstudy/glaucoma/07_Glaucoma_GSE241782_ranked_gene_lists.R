# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the saved glaucoma results to continue from the completed analyses.
load(file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma", "Glaucoma_GSE241782_withinstudy_results_checkpoint_2026-11-04.RData"))

# Use local repository output storage after restoring saved path objects.
glaucoma_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma")

# Allow only tiny negative F values within R's default numerical comparison tolerance.
# This is a numerical guard, not a differential-expression significance threshold.
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/all.equal.html
glaucoma_rank_tolerance <- sqrt(.Machine$double.eps)
stopifnot(
  all(is.finite(glaucoma_3D_results$F)),
  all(glaucoma_3D_results$F >= -glaucoma_rank_tolerance),
  all(is.finite(glaucoma_2W_results$F)),
  all(glaucoma_2W_results$F >= -glaucoma_rank_tolerance),
  all(is.finite(glaucoma_6W_results$F)),
  all(glaucoma_6W_results$F >= -glaucoma_rank_tolerance)
)

# Create separate tables with signed ranking scores for every tested gene.
# logFC supplies the direction, while F supplies the strength of the test statistic.
# pmax(F, 0) treats the reviewed tiny negative inputs as zero only for ranking.
# Original F, logFC, PValue and FDR columns are retained without modification.
glaucoma_3D_ranked <- transform(
  glaucoma_3D_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)
glaucoma_2W_ranked <- transform(
  glaucoma_2W_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)
glaucoma_6W_ranked <- transform(
  glaucoma_6W_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)

# Check that all scores are usable numbers before sorting or exporting.
stopifnot(
  all(is.finite(glaucoma_3D_ranked$rank_score)),
  all(is.finite(glaucoma_2W_ranked$rank_score)),
  all(is.finite(glaucoma_6W_ranked$rank_score))
)

# Sort each table from the strongest positive score to the strongest negative score.
glaucoma_3D_ranked <- glaucoma_3D_ranked[
  order(glaucoma_3D_ranked$rank_score, decreasing = TRUE),
]
glaucoma_2W_ranked <- glaucoma_2W_ranked[
  order(glaucoma_2W_ranked$rank_score, decreasing = TRUE),
]
glaucoma_6W_ranked <- glaucoma_6W_ranked[
  order(glaucoma_6W_ranked$rank_score, decreasing = TRUE),
]

# Set the three ranked-table filenames inside the glaucoma output folder.
glaucoma_ranked_files <- file.path(
  glaucoma_within_dir,
  c(
    "Glaucoma_GSE241782_3D_vs_Naive_ranked_all_genes.csv",
    "Glaucoma_GSE241782_2W_vs_Naive_ranked_all_genes.csv",
    "Glaucoma_GSE241782_6W_vs_Naive_ranked_all_genes.csv"
  )
)

# Save all original result columns plus the new ranking score for each stage.
write.csv(glaucoma_3D_ranked, file = glaucoma_ranked_files[1], row.names = FALSE)
write.csv(glaucoma_2W_ranked, file = glaucoma_ranked_files[2], row.names = FALSE)
write.csv(glaucoma_6W_ranked, file = glaucoma_ranked_files[3], row.names = FALSE)

# Check how many genes remain in each ranked table and whether all files were saved.
print(c(
  `3D` = nrow(glaucoma_3D_ranked),
  `2W` = nrow(glaucoma_2W_ranked),
  `6W` = nrow(glaucoma_6W_ranked)
))
print(file.exists(glaucoma_ranked_files))

# Preview the first five and last five genes in the 3-day ranked list.
print(head(glaucoma_3D_ranked[, c("gene_symbol", "logFC", "F", "rank_score", "FDR")], 5))
print(tail(glaucoma_3D_ranked[, c("gene_symbol", "logFC", "F", "rank_score", "FDR")], 5))
