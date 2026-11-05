# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the completed Autoimmune uveitis GSE284505 analysis without repeating modelling.
uveitis_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "autoimmune_uveitis")
load(file.path(uveitis_within_dir, "Autoimmune_Uveitis_GSE284505_withinstudy_final_checkpoint_2026-11-05.RData"))

# Use local repository output storage after restoring saved path objects.
uveitis_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "autoimmune_uveitis")

# Check sample order and display the saved result summary and full ranking size.
stopifnot(identical(colnames(uveitis_within_dge), uveitis_within_metadata$sample_id))
print(uveitis_summary)
print(nrow(uveitis_ranked))
