# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Display and save the limma mean-variance diagnostic from the moderated AMD comparison fit.
# This module reuses existing estimates and does not repeat fitting or empirical Bayes moderation.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

# Confirm that the moderated comparison fit and AMD output folder are available.
stopifnot(
  exists("amd_contrast_fit"),
  !is.null(amd_contrast_fit$s2.prior),
  dir.exists(amd_within_dir)
)

# Keep the requested blue for unusual feature variances, with regular estimates in black.
# The fitted variance-trend line keeps limma's default darker blue.
amd_sa_colours <- c("black", "#5E97F3")

# Display the diagnostic using the saved log2 microarray expression scale.
# The Normal/Outlier legend refers to feature variability, not sample groups or DEG status.
limma::plotSA(
  amd_contrast_fit,
  main = "CFH-deficiency model: mean-variance diagnostic",
  xlab = "Average log2 expression",
  ylab = "Quarter-root residual variance",
  col = amd_sa_colours
)

# Set the diagnostic figure path inside the AMD within-study output folder.
amd_sa_plot_file <- file.path(
  amd_within_dir, "AMD_GSE38671_mean_variance_diagnostic.png"
)

# Save the same diagnostic with the established figure dimensions and resolution.
png(amd_sa_plot_file, width = 2400, height = 1800, res = 300)
limma::plotSA(
  amd_contrast_fit,
  main = "CFH-deficiency model: mean-variance diagnostic",
  xlab = "Average log2 expression",
  ylab = "Quarter-root residual variance",
  col = amd_sa_colours
)
dev.off()

# Confirm that R can find the saved diagnostic figure.
print(file.exists(amd_sa_plot_file))
