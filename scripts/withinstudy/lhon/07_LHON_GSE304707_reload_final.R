# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the saved final LHON gene-level analysis without repeating calculations.

# Locate the checkpoint created by the final saving script.
lhon_final_checkpoint <- file.path(portfolio_repo_dir, "results", "withinstudy", "lhon", "LHON_GSE304707_withinstudy_final_checkpoint_2026-11-05.RData")
stopifnot(file.exists(lhon_final_checkpoint))
load(lhon_final_checkpoint)

# Use local repository output storage after restoring saved path objects.
lhon_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "lhon")

# Return to the LHON output folder and display the restored comparison summary.
setwd(lhon_within_dir)
print(getwd())
print(lhon_genotype_summary)

# Confirm that every complete ranked list was restored.
print(c(
  YARS2 = nrow(lhon_YARS2_ranked),
  COI = nrow(lhon_COI_ranked),
  YARS2_COI = nrow(lhon_YARS2_COI_ranked)
))

# Confirm that restored metadata rows still match the count columns.
print(identical(
  lhon_within_metadata$sample_id,
  colnames(lhon_within_dge$counts)
))
