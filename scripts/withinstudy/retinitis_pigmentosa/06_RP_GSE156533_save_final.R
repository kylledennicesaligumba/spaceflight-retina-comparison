# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Save the completed RP gene-level analysis from the R session used for the walkthrough.
# This preserves existing objects without repeating fitting, tests or plotting.

# Confirm that the output folder and completed model/ranking are available.
rp_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "retinitis_pigmentosa")
stopifnot(
  dir.exists(rp_within_dir),
  exists("rp_fit"),
  exists("rp_results"),
  exists("rp_ranked")
)

# Record the R and package versions used at the time of saving.
rp_within_session <- sessionInfo()

# Give the completed analysis a dated checkpoint in the RP output folder.
rp_final_checkpoint <- file.path(
  rp_within_dir,
  "RP_GSE156533_withinstudy_final_checkpoint_2026-11-05.RData"
)

# Save the working data, design, model, tests, tables, figures and full ranking.
# Keep the original precheck environment and unrelated datasets out of this file.
save(
  rp_within_dir, rp_precheck_file, rp_final_checkpoint, rp_within_session,
  rp_within_dge, rp_within_metadata, rp_design, rp_fit, rp_ql_plot_file,
  rp_test, rp_results, rp_degs, rp_summary, rp_result_files,
  rp_md_data, rp_md_hits, rp_md_labels, rp_md_plot, rp_md_file,
  rp_heat_genes, rp_heat_log_cpm, rp_heat_samples, rp_heat_z, rp_heat_limit,
  rp_heat_annotations, rp_heat_row_annotations, rp_heat_annotation_colours,
  rp_heatmap, rp_heatmap_files,
  rp_dot_genes, rp_dot_cpm, rp_dot_data, rp_dot_gene_labels, rp_dotplot,
  rp_dotplot_file, rp_dotplot_data_file,
  rp_ranked, rp_ranked_file, rp_rank_tolerance, rp_rank_F_check,
  file = rp_final_checkpoint
)

# Save a readable record of the R and package versions beside the checkpoint.
capture.output(
  rp_within_session,
  file = file.path(
    rp_within_dir,
    "RP_GSE156533_sessionInfo_final_2026-11-05.txt"
  )
)

# Confirm that the checkpoint exists and report its size in bytes.
print(file.exists(rp_final_checkpoint))
print(file.info(rp_final_checkpoint)$size)
