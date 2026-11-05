# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Glaucoma GSE241782: restore the saved within-study model.
# Prepared: 4 November 2026. Run this block and review the output before continuing.

# Restore the checkpoint containing the prepared counts, metadata, design and fit.
load(file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma", "Glaucoma_GSE241782_withinstudy_checkpoint_2026-11-03.RData"))

# Use local repository output storage after restoring saved path objects.
glaucoma_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma")

# Set the saved glaucoma within-study folder as the working location.
setwd(glaucoma_within_dir)

# Check the restored model's gene rows and coefficient columns.
dim(glaucoma_fit$coefficients)

# Confirm that the design rows still match the count columns in order.
identical(rownames(glaucoma_design), colnames(glaucoma_within_dge$counts))
