# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Save the completed OSD-758 within-study analysis from the walkthrough session.
# Preserve existing results without repeating model fitting, tests or plotting.

# Confirm that the completed model, results, ranking and sample alignment are available.
osd758_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758")
stopifnot(
  dir.exists(osd758_within_dir),
  exists("osd758_fit"), exists("osd758_results"), exists("osd758_ranked"),
  identical(osd758_within_metadata$sample, colnames(osd758_within_dge$counts)),
  nrow(osd758_ranked) == nrow(osd758_results)
)

# Confirm that the previously generated result tables and figures are present.
osd758_completed_output_files <- c(
  osd758_result_files, osd758_ql_plot_file, osd758_md_file,
  osd758_heatmap_files, osd758_dotplot_files, osd758_ranked_file
)
stopifnot(
  all(file.exists(osd758_completed_output_files)),
  all(file.info(osd758_completed_output_files)$size > 0)
)

# Preserve the reviewed comparison, significance rule and retained sample flag.
osd758_within_review <- list(
  comparison = "FLT_0G_vs_HGC",
  DEG_rule = "BH FDR < 0.05; no additional fold-change cutoff",
  ranking_rule = "sign(logFC) * sqrt(pmax(F, 0))",
  retained_projection_flag = "MHU-8_RTN_HGC_HGC05",
  excluded_samples = character(0),
  interpretation = "Flight-associated mouse retinal expression response"
)

# Record the R and package versions used at the time of saving.
osd758_within_session <- sessionInfo()
osd758_final_checkpoint <- file.path(
  osd758_within_dir, "OSD758_withinstudy_final_checkpoint_2026-11-06.RData"
)
osd758_final_session_file <- file.path(
  osd758_within_dir, "OSD758_sessionInfo_final_2026-11-06.txt"
)

# Save the working inputs, model, tests, tables, plotted data, figures and full ranking.
# Keep the original precheck environment and unrelated datasets out of this checkpoint.
save(
  osd758_within_dir, osd758_precheck_file, osd758_final_checkpoint,
  osd758_final_session_file, osd758_within_session, osd758_within_review,
  osd758_completed_output_files,
  osd758_within_dge, osd758_within_metadata, osd758_design, osd758_fit,
  osd758_ql_plot_file, osd758_test, osd758_results, osd758_degs,
  osd758_summary, osd758_result_files,
  osd758_md_data, osd758_md_hits, osd758_md_labels, osd758_md_plot, osd758_md_file,
  osd758_heat_genes, osd758_heat_log_cpm, osd758_heat_samples,
  osd758_heat_z, osd758_heat_limit, osd758_heat_annotations,
  osd758_heat_row_annotations, osd758_heat_annotation_colours,
  osd758_heatmap, osd758_heatmap_files,
  osd758_dot_genes, osd758_dot_cpm, osd758_dot_data,
  osd758_dotplot, osd758_dotplot_files,
  osd758_ranked, osd758_ranked_file, osd758_rank_tolerance, osd758_rank_F_check,
  file = osd758_final_checkpoint
)

# Save a readable record of the R and package versions beside the checkpoint.
capture.output(osd758_within_session, file = osd758_final_session_file)

# Confirm that both final files exist and report the checkpoint size in bytes.
print(file.exists(c(osd758_final_checkpoint, osd758_final_session_file)))
print(file.info(osd758_final_checkpoint)$size)
