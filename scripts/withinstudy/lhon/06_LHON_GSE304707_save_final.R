# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Save the completed LHON gene-level analysis from the current walkthrough session.
# This preserves existing objects without repeating fitting, tests, ranking or plotting.

# Confirm that the output folder, fitted model and three completed rankings are available.
lhon_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "lhon")
stopifnot(
  dir.exists(lhon_within_dir), exists("lhon_fit"),
  exists("lhon_YARS2_ranked"), exists("lhon_COI_ranked"),
  exists("lhon_YARS2_COI_ranked")
)

# Record the R and package versions used at the time of saving.
lhon_within_session <- sessionInfo()

# Give the completed gene-level analysis a dated checkpoint in the LHON output folder.
lhon_final_checkpoint <- file.path(
  lhon_within_dir,
  "LHON_GSE304707_withinstudy_final_checkpoint_2026-11-05.RData"
)

# Save the working data, model, all three tests/results, figures and full rankings.
# Keep the original precheck environment and unrelated study objects out of this file.
save(
  lhon_within_dir, lhon_precheck_file, lhon_final_checkpoint, lhon_within_session,
  lhon_within_dge, lhon_within_metadata, lhon_design, lhon_fit, lhon_ql_plot_file,
  lhon_YARS2_test, lhon_COI_test, lhon_YARS2_COI_test,
  lhon_YARS2_results, lhon_COI_results, lhon_YARS2_COI_results,
  lhon_YARS2_degs, lhon_COI_degs, lhon_YARS2_COI_degs,
  lhon_YARS2_summary, lhon_COI_summary, lhon_YARS2_COI_summary,
  lhon_YARS2_result_files, lhon_COI_result_files, lhon_YARS2_COI_result_files,
  lhon_genotype_summary, lhon_summary_file,
  lhon_md_data, lhon_md_hits, lhon_md_labels, lhon_md_plot, lhon_md_file,
  lhon_heat_genes, lhon_heat_log_cpm, lhon_heat_samples, lhon_heat_z, lhon_heat_limit,
  lhon_heat_annotations, lhon_heat_row_annotations, lhon_heat_annotation_colours,
  lhon_heatmap, lhon_heatmap_files,
  lhon_dot_higher, lhon_dot_lower, lhon_dot_names, lhon_dot_genes,
  lhon_dot_cpm, lhon_dot_data, lhon_dotplot, lhon_dotplot_files,
  lhon_YARS2_ranked, lhon_COI_ranked, lhon_YARS2_COI_ranked,
  lhon_rank_F_check, lhon_rank_tolerance, lhon_ranked_files,
  file = lhon_final_checkpoint
)

# Save a readable record of the R and package versions beside the checkpoint.
capture.output(
  lhon_within_session,
  file = file.path(
    lhon_within_dir,
    "LHON_GSE304707_sessionInfo_final_2026-11-05.txt"
  )
)

# Confirm that the checkpoint exists and report its size in bytes.
print(file.exists(lhon_final_checkpoint))
print(file.info(lhon_final_checkpoint)$size)
