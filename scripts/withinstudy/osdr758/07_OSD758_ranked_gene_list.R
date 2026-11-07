# OSD-758: full signed gene ranking for flight 0G versus HGC retina.
# Use the same project ranking rule and numerical guard as RP and glaucoma.

# Confirm that the complete results and output folder are available.
stopifnot(exists("osd758_results"), dir.exists(osd758_within_dir))

# Report the original F statistics before applying the existing numerical guard.
osd758_rank_F_check <- c(
  Minimum_F = min(osd758_results$F, na.rm = TRUE),
  Negative_F = sum(osd758_results$F < 0, na.rm = TRUE),
  Nonfinite_F = sum(!is.finite(osd758_results$F))
)
print(osd758_rank_F_check)

# Allow only tiny negative F values within R's default comparison tolerance.
# This guard handles numerical precision and does not change significance testing.
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/all.equal.html
osd758_rank_tolerance <- sqrt(.Machine$double.eps)
stopifnot(
  all(is.finite(osd758_results$F)),
  all(osd758_results$F >= -osd758_rank_tolerance),
  all(is.finite(osd758_results$logFC))
)

# Score every tested gene: logFC supplies direction and sqrt(F) supplies magnitude.
# Tiny negative inputs become zero for ranking; the original F column is retained.
osd758_ranked <- transform(
  osd758_results, rank_score = sign(logFC) * sqrt(pmax(F, 0))
)
stopifnot(
  all(is.finite(osd758_ranked$rank_score)),
  nrow(osd758_ranked) == nrow(osd758_results),
  anyDuplicated(osd758_ranked$gene_id) == 0L
)

# Sort from the strongest positive score to the most negative score.
osd758_ranked <- osd758_ranked[
  order(osd758_ranked$rank_score, decreasing = TRUE), , drop = FALSE
]

# Save all tested genes with their original identifiers and result columns.
osd758_ranked_file <- file.path(
  osd758_within_dir,
  "OSD758_FLT_0G_vs_HGC_ranked_all_genes.csv"
)
write.csv(osd758_ranked, file = osd758_ranked_file, row.names = FALSE)

# Confirm the number of ranked genes and that the file was saved.
print(nrow(osd758_ranked))
print(file.exists(osd758_ranked_file))

# Preview the first five and last five genes using the original Ensembl IDs.
print(head(osd758_ranked[, c("gene_id", "logFC", "F", "rank_score", "FDR")], 5))
print(tail(osd758_ranked[, c("gene_id", "logFC", "F", "rank_score", "FDR")], 5))
