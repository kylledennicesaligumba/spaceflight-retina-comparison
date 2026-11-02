# Save the current uveitis precheck from the R session used for this walkthrough.
# Existing counts, annotations, normalization results and plots are preserved.

# Require the completed-stage objects before creating the new checkpoint.
uveitis_save_required <- c(
  "uveitis_output_dir", "uveitis_raw", "uveitis_counts_all",
  "uveitis_metadata_all", "uveitis_metadata_matched",
  "uveitis_counts_selected", "uveitis_metadata_selected",
  "uveitis_gene_annotation", "uveitis_ensembl_annotation",
  "uveitis_all_lookup_result", "uveitis_qc_group", "uveitis_dge",
  "uveitis_keep_genes", "uveitis_dge_qc", "uveitis_log_cpm",
  "uveitis_log_cpm_qc", "uveitis_mds", "uveitis_pca",
  "uveitis_mds_plot", "uveitis_pca_plot"
)
uveitis_save_missing <- uveitis_save_required[
  !vapply(uveitis_save_required, exists, logical(1),
          envir = .GlobalEnv, inherits = FALSE)
]
if (length(uveitis_save_missing)) {
  stop("Missing current-stage objects: ", paste(uveitis_save_missing, collapse = ", "))
}

# Confirm the current dimensions and sample/gene alignment before saving.
stopifnot(
  dir.exists(uveitis_output_dir),
  identical(dim(uveitis_counts_all), c(56859L, 12L)),
  identical(dim(uveitis_counts_selected), c(56859L, 6L)),
  identical(dim(uveitis_gene_annotation), c(56859L, 9L)),
  identical(dim(uveitis_dge_qc), c(18187L, 6L)),
  identical(dim(uveitis_log_cpm_qc), c(18187L, 6L)),
  identical(colnames(uveitis_counts_selected), uveitis_metadata_selected$sample_id),
  identical(uveitis_gene_annotation$AccID, rownames(uveitis_counts_selected)),
  identical(uveitis_dge_qc$samples$sample_id, colnames(uveitis_dge_qc)),
  identical(uveitis_dge_qc$genes$AccID, rownames(uveitis_dge_qc)),
  identical(colnames(uveitis_log_cpm_qc), colnames(uveitis_dge_qc))
)

# Set new dated filenames in the established output folder.
uveitis_checkpoint_file <- file.path(
  uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_precheck_checkpoint_2026-11-02.RData"
)
uveitis_saved_files <- c(
  checkpoint = uveitis_checkpoint_file,
  metadata = file.path(uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_metadata_selected_2026-11-02.csv"),
  annotation = file.path(uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_gene_annotation_2026-11-02.csv"),
  TMM_summary = file.path(uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_TMM_summary_2026-11-02.csv"),
  session_info = file.path(uveitis_output_dir, "Autoimmune_Uveitis_GSE284505_sessionInfo_2026-11-02.txt")
)
uveitis_tmm_summary <- uveitis_dge_qc$samples[
  , c("sample_id", "group", "lib.size", "norm.factors"), drop = FALSE
]
uveitis_session_info <- utils::sessionInfo()

# Save all uveitis objects except live XML readers and the HTTP response.
# Parsed metadata and Ensembl results are retained; XML readers can be reopened.
uveitis_checkpoint_objects <- setdiff(
  ls(pattern = "^uveitis_", envir = .GlobalEnv),
  c("uveitis_metadata_xml", "uveitis_sample_nodes",
    "uveitis_characteristic_nodes", "uveitis_protocol_nodes",
    "uveitis_lookup_response")
)
save(list = uveitis_checkpoint_objects, file = uveitis_checkpoint_file,
     envir = .GlobalEnv)

# Export readable tables and the package versions from this actual R session.
write.csv(uveitis_metadata_selected, uveitis_saved_files[["metadata"]], row.names = FALSE)
write.csv(uveitis_gene_annotation, uveitis_saved_files[["annotation"]], row.names = FALSE)
write.csv(uveitis_tmm_summary, uveitis_saved_files[["TMM_summary"]], row.names = FALSE)
writeLines(capture.output(print(uveitis_session_info)), uveitis_saved_files[["session_info"]])

# Confirm all five outputs exist; rerunning updates these 2 November files.
print(setNames(file.exists(uveitis_saved_files), names(uveitis_saved_files)))
cat("Checkpoint:", uveitis_checkpoint_file, "\n")
