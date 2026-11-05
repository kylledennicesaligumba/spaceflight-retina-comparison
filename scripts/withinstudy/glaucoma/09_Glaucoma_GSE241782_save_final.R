# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Confirm that the glaucoma output folder and completed ranked tables are available.
glaucoma_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "glaucoma")
stopifnot(
  dir.exists(glaucoma_within_dir),
  exists("glaucoma_3D_ranked"),
  exists("glaucoma_2W_ranked"),
  exists("glaucoma_6W_ranked")
)

# Record the R and package versions at completion of the gene-level analysis.
glaucoma_within_session <- sessionInfo()

# Give the completed gene-level checkpoint its own filename in the glaucoma folder.
glaucoma_final_checkpoint <- file.path(
  glaucoma_within_dir,
  "Glaucoma_GSE241782_withinstudy_final_checkpoint_2026-11-04.RData"
)

# Save prepared data, the model, all stage results, plots and full ranked tables.
save(
  glaucoma_within_dir, glaucoma_final_checkpoint, glaucoma_within_session,
  glaucoma_within_dge, glaucoma_within_metadata,
  glaucoma_design, glaucoma_fit,
  glaucoma_3D_test, glaucoma_2W_test, glaucoma_6W_test,
  glaucoma_3D_results, glaucoma_2W_results, glaucoma_6W_results,
  glaucoma_3D_degs, glaucoma_2W_degs, glaucoma_6W_degs,
  glaucoma_3D_summary, glaucoma_stage_summary,
  glaucoma_deg_ids, glaucoma_deg_counts,
  glaucoma_cpm_within, glaucoma_deg_cpm,
  glaucoma_plot_colours, glaucoma_plot_labels,
  glaucoma_dot_data, glaucoma_dotplot,
  glaucoma_heat_log_cpm, glaucoma_heat_z,
  glaucoma_heat_samples, glaucoma_heat_symbols,
  glaucoma_heat_annotations, glaucoma_heat_annotation_colours,
  glaucoma_heatmap,
  glaucoma_md_data, glaucoma_md_hits, glaucoma_md_plot,
  glaucoma_3D_ranked, glaucoma_2W_ranked, glaucoma_6W_ranked,
  glaucoma_rank_tolerance, glaucoma_ranked_files,
  file = glaucoma_final_checkpoint
)

# Save a readable record of the R and package versions alongside the checkpoint.
capture.output(
  glaucoma_within_session,
  file = file.path(
    glaucoma_within_dir,
    "Glaucoma_GSE241782_sessionInfo_final_2026-11-04.txt"
  )
)

# Confirm that the final checkpoint exists and report its size in bytes.
print(file.exists(glaucoma_final_checkpoint))
print(file.info(glaucoma_final_checkpoint)$size)
