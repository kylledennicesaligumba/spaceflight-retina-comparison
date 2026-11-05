# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the completed Light damage GSE314043 analysis without repeating modelling.
ld_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "light_damage")
load(file.path(ld_within_dir, "Light_Damage_GSE314043_withinstudy_final_checkpoint_2026-11-05.RData"))

# Use local repository output storage after restoring saved path objects.
ld_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "light_damage")

# Check sample order and display the saved result summary and full ranking size.
stopifnot(identical(colnames(ld_within_dge), ld_within_metadata$sample_id))
print(ld_summary)
print(nrow(ld_ranked))
