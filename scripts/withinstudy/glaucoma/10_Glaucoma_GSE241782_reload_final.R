# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the completed gene-level checkpoint when returning to glaucoma later.
glaucoma_final_checkpoint <- file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma", "Glaucoma_GSE241782_withinstudy_final_checkpoint_2026-11-04.RData")
stopifnot(file.exists(glaucoma_final_checkpoint))
load(glaucoma_final_checkpoint)

# Use local repository output storage after restoring saved path objects.
glaucoma_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma")

# Return to the glaucoma output folder and display the saved stage summary.
setwd(glaucoma_within_dir)
print(getwd())
print(glaucoma_stage_summary)

# Confirm the number of genes restored in each complete ranked table.
print(c(
  `3D` = nrow(glaucoma_3D_ranked),
  `2W` = nrow(glaucoma_2W_ranked),
  `6W` = nrow(glaucoma_6W_ranked)
))
