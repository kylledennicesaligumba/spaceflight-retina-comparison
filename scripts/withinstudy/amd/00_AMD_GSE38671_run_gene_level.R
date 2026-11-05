# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Prepare local output storage before the original workflow checks it.
dir.create(file.path(portfolio_repo_dir, "results", "withinstudy", "amd"), recursive = TRUE, showWarnings = FALSE)

# Reproduce the authorized AMD gene-level analysis, five figures, ranking and final save.
# Run this script with Rscript, or source it in a new R session.
# Earlier manual feature-level scripts/results remain separate and are not overwritten.
amd_script_dir <- file.path(portfolio_repo_dir, "scripts", "withinstudy", "amd")
source(file.path(amd_script_dir, "03_AMD_GSE38671_gene_level_analysis.R"))
source(file.path(amd_script_dir, "04_AMD_GSE38671_gene_level_plots.R"))
source(file.path(amd_script_dir, "05_AMD_GSE38671_gene_level_rank_save.R"))
