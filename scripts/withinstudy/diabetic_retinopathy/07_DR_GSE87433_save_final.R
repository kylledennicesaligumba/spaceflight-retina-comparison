# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Save the completed DR gene-level analysis from the current walkthrough session.
# Preserve existing objects without repeating fitting, testing, ranking or plotting.

# Confirm that the output folder, fitted model and completed ranking are available.
dr_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "diabetic_retinopathy")
stopifnot(
  dir.exists(dr_within_dir), exists("dr_fit"),
  exists("dr_results"), exists("dr_ranked"),
  nrow(dr_ranked) == nrow(dr_results),
  all(is.finite(dr_ranked$rank_score))
)

# Record the R and package versions used at the time of saving.
dr_within_session <- sessionInfo()

# Give the completed analysis a dated checkpoint in the DR output folder.
dr_final_checkpoint <- file.path(
  dr_within_dir, "DR_GSE87433_withinstudy_final_checkpoint_2026-11-05.RData"
)

# Record the existing diagnostic figure's path alongside the other figure paths.
dr_sa_plot_file <- file.path(dr_within_dir, "DR_GSE87433_mean_variance_diagnostic.png")

# Save full and selected working data, aligned annotation, model, results, figures and ranking.
# Keep the original precheck environment and unrelated study objects out of this checkpoint.
save(
  dr_within_dir, dr_precheck_file, dr_final_checkpoint, dr_within_session,
  dr_within_expr, dr_within_metadata, dr_within_annotation,
  dr_primary_expr, dr_primary_metadata, dr_design, dr_fit, dr_sa_plot_file,
  dr_results, dr_degs, dr_summary, dr_result_files,
  dr_md_data, dr_md_hits, dr_md_plot, dr_md_file,
  dr_heat_genes, dr_heat_expr, dr_heat_samples, dr_heat_z, dr_heat_limit,
  dr_heat_annotations, dr_heat_row_annotations, dr_heat_annotation_colours,
  dr_heatmap, dr_heatmap_files,
  dr_dot_genes, dr_dot_expr, dr_dot_data, dr_dot_gene_labels, dr_dotplot,
  dr_dotplot_file, dr_dotplot_data_file,
  dr_ranked, dr_ranked_file, dr_rank_check,
  file = dr_final_checkpoint
)

# Save a readable record of the R and package versions beside the checkpoint.
capture.output(
  dr_within_session,
  file = file.path(dr_within_dir, "DR_GSE87433_sessionInfo_final_2026-11-05.txt")
)

# Confirm that the checkpoint exists and report its size in bytes.
print(file.exists(dr_final_checkpoint))
print(file.info(dr_final_checkpoint)$size)
