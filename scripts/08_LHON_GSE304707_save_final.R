# Portable repository copy; run from the repository root.
# Original working files are preserved. GitHub sync did not execute this code.
# Save the completed LHON processed-data precheck from the current R session.
# Run this file after completing TMM, MDS and PCA. It saves existing results;
# it does not repeat filtering, normalization, projections or figure exports.

# Require the essential completed-stage objects before writing a final save.
lhon_final_required_objects <- c(
  "lhon_output_dir", "lhon_raw", "lhon_counts_all",
  "lhon_metadata_all", "lhon_metadata_matched", "lhon_gene_annotation",
  "lhon_qc_group", "lhon_dge", "lhon_keep_genes", "lhon_dge_qc",
  "lhon_log_cpm", "lhon_log_cpm_qc", "lhon_total_counts",
  "lhon_detected_genes", "lhon_mds", "lhon_pca",
  "lhon_mds_plot", "lhon_pca_plot"
)
lhon_final_missing_objects <- lhon_final_required_objects[
  !vapply(
    lhon_final_required_objects, exists, logical(1),
    envir = .GlobalEnv, inherits = FALSE
  )
]
if (length(lhon_final_missing_objects)) {
  stop("Missing completed-stage objects: ", paste(lhon_final_missing_objects, collapse = ", "))
}

# Confirm the completed object sizes and their current sample/gene alignment.
stopifnot(
  dir.exists(lhon_output_dir),
  identical(dim(lhon_counts_all), c(49714L, 12L)),
  identical(dim(lhon_dge_qc), c(17860L, 12L)),
  identical(dim(lhon_log_cpm_qc), c(17860L, 12L)),
  identical(lhon_dge_qc$samples$sample_id, colnames(lhon_dge_qc)),
  identical(lhon_dge_qc$genes$gene_name, rownames(lhon_dge_qc)),
  identical(colnames(lhon_log_cpm_qc), colnames(lhon_dge_qc)),
  all(is.finite(lhon_log_cpm_qc)),
  all(is.finite(lhon_dge_qc$samples$norm.factors)),
  all(lhon_dge_qc$samples$norm.factors > 0)
)

# Preserve a readable copy of the already calculated normalization factors.
lhon_tmm_summary <- lhon_dge_qc$samples[
  , c("sample_id", "group", "sex", "lib.size", "norm.factors"), drop = FALSE
]

# Use new dated files, preserving the earlier 29 October checkpoint.
lhon_final_checkpoint_file <- file.path(
  lhon_output_dir, "LHON_GSE304707_precheck_checkpoint_2026-10-30.RData"
)
lhon_final_metadata_file <- file.path(
  lhon_output_dir, "LHON_GSE304707_metadata_matched_2026-10-30.csv"
)
lhon_final_tmm_file <- file.path(
  lhon_output_dir, "LHON_GSE304707_TMM_summary_2026-10-30.csv"
)
lhon_final_session_file <- file.path(
  lhon_output_dir, "LHON_GSE304707_sessionInfo_2026-10-30.txt"
)

# Save explicit LHON data/text/plot objects only. XML document/node pointers
# are reopened from the saved XML path when needed, rather than serialized.
lhon_final_checkpoint_objects <- unique(c(
  lhon_final_required_objects,
  "lhon_dir", "lhon_count_file", "lhon_metadata_archive",
  "lhon_metadata_dir", "lhon_metadata_file", "lhon_ns", "lhon_sample_ids",
  "lhon_gene_labels", "lhon_ensembl_prefix", "lhon_characteristic_tags",
  "lhon_extraction_protocol", "lhon_get_characteristic",
  "lhon_draw_boxplot", "lhon_boxplot_file",
  "lhon_density_curves", "lhon_density_colors", "lhon_draw_density", "lhon_density_file",
  "lhon_draw_mds", "lhon_draw_pca", "lhon_mds_file", "lhon_pca_file",
  "lhon_projection_plot", "lhon_tmm_summary",
  "lhon_final_checkpoint_file", "lhon_final_metadata_file",
  "lhon_final_tmm_file", "lhon_final_session_file",
  "lhon_final_checkpoint_objects"
))
lhon_final_checkpoint_objects <- intersect(
  lhon_final_checkpoint_objects, ls(envir = .GlobalEnv, all.names = TRUE)
)
save(list = lhon_final_checkpoint_objects, file = lhon_final_checkpoint_file)

# Export the current matched metadata, factors and actual session information.
write.csv(lhon_metadata_matched, lhon_final_metadata_file, row.names = FALSE)
write.csv(lhon_tmm_summary, lhon_final_tmm_file, row.names = FALSE)
writeLines(capture.output(sessionInfo()), lhon_final_session_file)

# Confirm successful file creation. Rerunning updates these final dated files.
lhon_final_save_status <- c(
  checkpoint = file.exists(lhon_final_checkpoint_file),
  metadata = file.exists(lhon_final_metadata_file),
  TMM_summary = file.exists(lhon_final_tmm_file),
  session_info = file.exists(lhon_final_session_file)
)
print(lhon_final_save_status)
cat("Final checkpoint:", lhon_final_checkpoint_file, "\n")
# Expected: TRUE for all four. This script does not close R or eject the drive.
