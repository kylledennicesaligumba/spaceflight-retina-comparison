# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Display and save the DR limma SA diagnostic from the existing moderated fit.
# This script changes the figure's presentation without repeating model fitting.
# Method: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

# Confirm that the current fitted model and DR output folder are available.
stopifnot(exists("dr_fit"), dir.exists(dr_within_dir))

# Use the blue sampled from the user's reference screenshot for unusual gene variances.
# plotSA's two colours refer to regular and unusual gene variances respectively.
dr_sa_colours <- c("black", "#5E97F3")

# Display the limma diagnostic with its method identified in the title.
# The legend's Normal/Outlier describes gene variability, not sample groups.
limma::plotSA(
  dr_fit, main = "DR: limma mean-variance diagnostic",
  xlab = "Average processed expression", ylab = "Quarter-root residual variance",
  col = dr_sa_colours
)

# Use the same diagnostic filename in the DR within-study folder.
dr_sa_plot_file <- file.path(
  dr_within_dir, "DR_GSE87433_mean_variance_diagnostic.png"
)

# Save the revised colours and title using the existing fitted model.
png(dr_sa_plot_file, width = 2400, height = 1800, res = 300)
limma::plotSA(
  dr_fit, main = "DR: limma mean-variance diagnostic",
  xlab = "Average processed expression", ylab = "Quarter-root residual variance",
  col = dr_sa_colours
)
dev.off()

# Confirm that the updated diagnostic image was saved.
print(file.exists(dr_sa_plot_file))
