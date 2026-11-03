# Find the latest completed checkpoint, preserving the earlier fallback.
precheck_project_dir <- normalizePath(getwd(), mustWork = TRUE)
amd_latest_checkpoint_file <- file.path(precheck_project_dir, "results", "amd", "AMD_GSE38671_latest_checkpoint.txt")
amd_reload_checkpoint_file <- if (file.exists(amd_latest_checkpoint_file)) {
  readLines(amd_latest_checkpoint_file, n = 1L, warn = FALSE)
} else {
  file.path(precheck_project_dir, "results", "amd", "AMD_GSE38671_precheck_checkpoint_2026-11-02.RData")
}
stopifnot(length(amd_reload_checkpoint_file) == 1L,
          file.exists(amd_reload_checkpoint_file))
cat("Loading checkpoint:", amd_reload_checkpoint_file, "\n")
load(amd_reload_checkpoint_file, envir = .GlobalEnv)
setwd(amd_output_dir)

# Reopen the metadata readers from the preserved XML file.
stopifnot(file.exists(amd_metadata_file))
amd_metadata_xml <- xml2::read_xml(amd_metadata_file)
amd_sample_nodes <- xml2::xml_find_all(
  amd_metadata_xml, "/geo:MINiML/geo:Sample", ns = amd_ns
)
amd_characteristic_nodes <- xml2::xml_find_all(
  amd_sample_nodes, ".//geo:Characteristics", ns = amd_ns
)

amd_platform_nodes <- xml2::xml_find_all(
  amd_metadata_xml, "/geo:MINiML/geo:Platform[@iid='GPL6246']", ns = amd_ns
)
amd_platform_column_nodes <- xml2::xml_find_all(
  amd_platform_nodes, ".//geo:Data-Table/geo:Column", ns = amd_ns
)

# Check the restored dimensions against the user-reviewed results.
stopifnot(
  identical(dim(amd_expr_all), c(35556L, 32L)),
  nrow(amd_metadata_all) == 32L
)
print(getwd())
print(dim(amd_expr_all))
if (exists("amd_expr_log2", envir = .GlobalEnv, inherits = FALSE)) {
  print(dim(amd_expr_log2))
}
cat("AMD checkpoint restored. See docs/amd-precheck.md for the completed record and checkpoint scope.\n")
