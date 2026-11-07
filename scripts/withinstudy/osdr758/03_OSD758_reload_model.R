# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the fitted OSD-758 model after its checkpoint has been saved.
# This script loads saved objects without repeating filtering, TMM or model fitting.
osd758_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758")
osd758_model_checkpoint_file <- file.path(
  osd758_within_dir, "OSD758_withinstudy_model_checkpoint_2026-11-06.RData"
)
stopifnot(file.exists(osd758_model_checkpoint_file))
load(osd758_model_checkpoint_file, envir = .GlobalEnv)
osd758_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758")

# Check the restored model dimensions, terms and alignment with the count input.
print(dim(osd758_fit$coefficients))
print(colnames(osd758_fit$coefficients))
print(identical(rownames(osd758_design), colnames(osd758_within_dge$counts)))
stopifnot(
  identical(dim(osd758_fit$coefficients), c(21335L, 2L)),
  identical(colnames(osd758_fit$coefficients), c("(Intercept)", "groupFLT_0G")),
  identical(rownames(osd758_design), colnames(osd758_within_dge$counts))
)
