# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the reviewed OSD-758 input after restarting R or the computer.
# This reloads saved filtering/TMM results without rerunning analysis or models.

# Recreate the within-study paths in the fresh session.
osd758_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758")
dir.create(osd758_within_dir, recursive = TRUE, showWarnings = FALSE)
osd758_precheck_file <- Sys.getenv("OSD758_PRECHECK_FILE", unset = file.path(portfolio_repo_dir, "results", "osdr758", "OSD758_precheck_checkpoint_2026-11-06.RData"))
stopifnot(file.exists(osd758_precheck_file))

# Load the saved precheck into a separate environment.
osd758_precheck <- new.env()
osd758_loaded_objects <- load(osd758_precheck_file, envir = osd758_precheck)

# Restore the working count object and its matching sample descriptions.
osd758_within_dge <- osd758_precheck$osd758_dge_qc
osd758_within_metadata <- osd758_precheck$osd758_metadata_selected

# Check the restored dimensions, condition counts and sample order.
print(dim(osd758_within_dge$counts))
print(table(osd758_within_metadata$condition))
print(identical(
  colnames(osd758_within_dge$counts), osd758_within_metadata$sample
))

# Verify that saved TMM factors are usable; retain all 17 reviewed samples.
stopifnot(
  identical(dim(osd758_within_dge$counts), c(21335L, 17L)),
  identical(colnames(osd758_within_dge$counts), osd758_within_metadata$sample),
  all(is.finite(osd758_within_dge$samples$norm.factors)),
  all(osd758_within_dge$samples$norm.factors > 0)
)
