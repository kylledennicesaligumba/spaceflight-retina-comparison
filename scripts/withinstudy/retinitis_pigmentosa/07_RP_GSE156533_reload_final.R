# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the saved final RP gene-level analysis without repeating calculations.

# Locate the checkpoint created by the final saving script.
rp_final_checkpoint <- file.path(portfolio_repo_dir, "results", "withinstudy", "retinitis_pigmentosa", "RP_GSE156533_withinstudy_final_checkpoint_2026-11-05.RData")
stopifnot(file.exists(rp_final_checkpoint))
load(rp_final_checkpoint)

# Use local repository output storage after restoring saved path objects.
rp_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "retinitis_pigmentosa")

# Return to the RP output folder and display the restored result summary.
setwd(rp_within_dir)
print(getwd())
print(rp_summary)
print(nrow(rp_ranked))

# Confirm that the restored metadata rows still match the count columns.
print(identical(
  rp_within_metadata$sample_id,
  colnames(rp_within_dge$counts)
))
