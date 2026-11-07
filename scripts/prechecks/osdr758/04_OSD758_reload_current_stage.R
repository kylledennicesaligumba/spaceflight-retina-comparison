# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the saved OSD-758 precheck stage without rerunning the analysis.
osd758_output_dir <- file.path(portfolio_repo_dir, "results", "osdr758")
osd758_checkpoint_file <- file.path(
  osd758_output_dir, "OSD758_precheck_checkpoint_2026-11-06.RData"
)
stopifnot(file.exists(osd758_checkpoint_file))
load(osd758_checkpoint_file, envir = .GlobalEnv)
osd758_output_dir <- file.path(portfolio_repo_dir, "results", "osdr758")

# Reload plotting definitions and confirm the stored data and sample alignment.
source(file.path(portfolio_repo_dir, "scripts", "prechecks", "osdr758", "02_OSD758_projection_style.R"))
stopifnot(
  identical(dim(osd758_dge_qc), c(21335L, 17L)),
  identical(colnames(osd758_dge_qc), osd758_metadata_selected$sample),
  identical(dim(osd758_log_cpm_qc), dim(osd758_dge_qc)),
  identical(rownames(osd758_log_cpm_qc), rownames(osd758_dge_qc))
)

# Display the saved review status, including the retained HGC05 flag.
osd758_precheck_review
