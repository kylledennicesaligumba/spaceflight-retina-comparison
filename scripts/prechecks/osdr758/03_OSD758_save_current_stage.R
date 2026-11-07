# Save the reviewed OSD-758 computational precheck stage from the user's session.
# This script does not recompute counts, normalization, projections or models.

# Require the reviewed data and projections before making a checkpoint.
osd758_required_objects <- c(
  "osd758_output_dir", "osd758_metadata_selected", "osd758_counts_mouse",
  "osd758_counts_ercc", "osd758_dge_qc", "osd758_keep_genes",
  "osd758_log_cpm_qc", "osd758_mds", "osd758_pca"
)
osd758_missing_objects <- setdiff(
  osd758_required_objects, ls(envir = .GlobalEnv)
)
if (length(osd758_missing_objects) > 0) {
  stop("Missing current-stage objects: ", paste(osd758_missing_objects, collapse = ", "))
}
stopifnot(
  identical(dim(osd758_dge_qc), c(21335L, 17L)),
  identical(colnames(osd758_dge_qc), osd758_metadata_selected$sample),
  identical(rownames(osd758_log_cpm_qc), rownames(osd758_dge_qc)),
  all(is.finite(osd758_log_cpm_qc)),
  length(osd758_mds$x) == 17L,
  length(osd758_pca$x) == 17L
)

# Preserve the reviewed retention decision and the unresolved sample flag.
osd758_precheck_review <- list(
  stage = "Exploratory computational QC reviewed; study-design verification pending",
  retained_samples = as.character(osd758_metadata_selected$sample),
  excluded_samples = character(0),
  flagged_for_review = "MHU-8_RTN_HGC_HGC05",
  flag_reason = "Distinct second-axis position in both MDS and PCA; cause not established",
  retained_mouse_rows = 21335L,
  expression_filtered_mouse_rows = 35413L,
  separated_ercc_rows = 92L,
  projection_observation = paste(
    "Partial flight/ground overlap; FLT09 and FLT17 are furthest right.",
    "Flight samples do not form one compact cluster.",
    "Retain all 17 at this stage; review source annotations before modelling."
  )
)

# Set dated checkpoint and software-record filenames in the prequality folder.
osd758_checkpoint_file <- file.path(
  osd758_output_dir, "OSD758_precheck_checkpoint_2026-11-06.RData"
)
osd758_session_file <- file.path(
  osd758_output_dir, "OSD758_precheck_sessionInfo_2026-11-06.txt"
)

# Save only this dataset's prefixed objects from the current R session.
osd758_checkpoint_objects <- ls(envir = .GlobalEnv, pattern = "^osd758_")
save(
  list = osd758_checkpoint_objects,
  envir = .GlobalEnv,
  file = osd758_checkpoint_file
)

# Record the software versions actually loaded in the user's session.
writeLines(capture.output(sessionInfo()), osd758_session_file)

# Confirm that both files were written.
c(
  checkpoint = file.exists(osd758_checkpoint_file),
  session_info = file.exists(osd758_session_file)
)
