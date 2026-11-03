# Set the AMD output folder before saving.
if (!exists("precheck_project_dir", inherits = TRUE)) {
  precheck_project_dir <- normalizePath(getwd(), mustWork = TRUE)
}
amd_output_dir <- file.path(precheck_project_dir, "results", "amd")
dir.create(amd_output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(amd_output_dir)

# Give the saved files a new timestamp.
amd_save_stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
amd_script_file <- file.path(precheck_project_dir, "scripts", "16_AMD_GSE38671_precheck.R")
amd_script_backup_file <- file.path(
  amd_output_dir, paste0("AMD_GSE38671_precheck_script_", amd_save_stamp, ".R")
)
amd_checkpoint_file <- file.path(
  amd_output_dir, paste0("AMD_GSE38671_precheck_checkpoint_", amd_save_stamp, ".RData")
)
amd_history_file <- file.path(
  amd_output_dir, paste0("AMD_GSE38671_commands_", amd_save_stamp, ".Rhistory")
)
amd_session_info_file <- file.path(
  amd_output_dir, paste0("AMD_GSE38671_sessionInfo_", amd_save_stamp, ".txt")
)
amd_latest_checkpoint_file <- file.path(
  amd_output_dir, "AMD_GSE38671_latest_checkpoint.txt"
)

# Preserve a dated copy of the teaching script.
stopifnot(file.exists(amd_script_file), !file.exists(amd_checkpoint_file))
stopifnot(file.copy(amd_script_file, amd_script_backup_file, overwrite = FALSE))

# Save AMD objects; XML readers are reopened by the reload helper.
amd_save_objects <- setdiff(
  ls(pattern = "^amd_", envir = .GlobalEnv),
  c("amd_metadata_xml", "amd_sample_nodes", "amd_characteristic_nodes",
    "amd_platform_nodes", "amd_platform_column_nodes", "amd_save_objects")
)
save(
  list = amd_save_objects, file = amd_checkpoint_file,
  envir = .GlobalEnv, compress = TRUE
)
stopifnot(file.exists(amd_checkpoint_file))
writeLines(amd_checkpoint_file, amd_latest_checkpoint_file)

# Record package versions and the current console history separately.
writeLines(capture.output(sessionInfo()), amd_session_info_file)
amd_history_saved <- tryCatch(
  {
    utils::savehistory(amd_history_file)
    file.exists(amd_history_file)
  },
  error = function(e) {
    message("R objects and script were saved. Console history was unavailable: ",
            conditionMessage(e))
    FALSE
  }
)

# Display the saved locations and confirm which files exist.
amd_saved_files <- c(
  script = amd_script_backup_file,
  checkpoint = amd_checkpoint_file,
  history = amd_history_file,
  session_info = amd_session_info_file
)
print(amd_saved_files)
print(setNames(file.exists(amd_saved_files), names(amd_saved_files)))
