# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the saved final DR analysis without repeating any analysis calculations.

# Locate the checkpoint created by the final saving script.
dr_final_checkpoint <- file.path(portfolio_repo_dir, "results", "withinstudy", "diabetic_retinopathy", "DR_GSE87433_withinstudy_final_checkpoint_2026-11-05.RData")
stopifnot(file.exists(dr_final_checkpoint))
load(dr_final_checkpoint)

# Use local repository output storage after restoring saved path objects.
dr_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "diabetic_retinopathy")

# Return to the DR output folder and display the restored summary and ranked gene count.
setwd(dr_within_dir)
print(getwd())
print(dr_summary)
print(nrow(dr_ranked))

# Confirm that the primary sample order and gene-annotation order are still correct.
print(identical(
  dr_primary_metadata$geo_accession,
  colnames(dr_primary_expr)
))
print(identical(
  dr_within_annotation$ID,
  rownames(dr_primary_expr)
))
