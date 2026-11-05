# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Inspect the original F statistics after the ranking script's square-root warning.
# Run in the current R session, where the three complete result tables are loaded.
glaucoma_F_check <- rbind(
  `3D` = c(
    Minimum_F = min(glaucoma_3D_results$F, na.rm = TRUE),
    Negative_F = sum(glaucoma_3D_results$F < 0, na.rm = TRUE),
    Nonfinite_F = sum(!is.finite(glaucoma_3D_results$F))
  ),
  `2W` = c(
    Minimum_F = min(glaucoma_2W_results$F, na.rm = TRUE),
    Negative_F = sum(glaucoma_2W_results$F < 0, na.rm = TRUE),
    Nonfinite_F = sum(!is.finite(glaucoma_2W_results$F))
  ),
  `6W` = c(
    Minimum_F = min(glaucoma_6W_results$F, na.rm = TRUE),
    Negative_F = sum(glaucoma_6W_results$F < 0, na.rm = TRUE),
    Nonfinite_F = sum(!is.finite(glaucoma_6W_results$F))
  )
)
print(glaucoma_F_check)

# Display the five lowest F values in the 3-day table with their test results.
print(head(
  glaucoma_3D_results[
    order(glaucoma_3D_results$F),
    c("gene_symbol", "logFC", "F", "PValue", "FDR")
  ],
  5
))
