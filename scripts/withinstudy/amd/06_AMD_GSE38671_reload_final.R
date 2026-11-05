# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the completed AMD gene-level analysis without rerunning models or plots.
amd_final_checkpoint <- file.path(portfolio_repo_dir, "results", "withinstudy", "amd", "AMD_GSE38671_gene_level_withinstudy_final_checkpoint_2026-11-05.RData")
stopifnot(file.exists(amd_final_checkpoint))
load(amd_final_checkpoint)
# Use local repository output storage after restoring saved path objects.
amd_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "amd")

print(dim(amd_gene_expr))
print(amd_gene_summary)
