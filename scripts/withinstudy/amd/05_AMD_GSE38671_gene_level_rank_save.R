# Repository paths: begin in the repository root, or set RETINA_PROJECT_DIR.
# Existing model calculations are preserved; downloaded checkpoints stay local.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Save full signed rankings and the completed AMD gene-level checkpoint.
# Like DR, limma's moderated t statistic already supplies direction and magnitude.
stopifnot(exists("amd_gene_results"), exists("amd_figure_files"), all(file.exists(amd_figure_files)))
amd_gene_ranked <- lapply(amd_gene_results, function(x) {
  stopifnot(nrow(x) == nrow(amd_gene_expr), all(is.finite(x$t)), all(is.finite(x$logFC)))
  x$rank_score <- x$t
  x[order(-x$rank_score, x$entrez_id, method = "radix"), , drop = FALSE]
})
for (comparison in names(amd_gene_ranked)) {
  write.csv(amd_gene_ranked[[comparison]], file.path(amd_within_dir,
    paste0("AMD_GSE38671_gene_level_", comparison, "_ranked_all_genes.csv")), row.names = FALSE)
}

# Record the installed R/package versions and the fixed client-date analysis identifier.
amd_analysis_date <- "2026-11-05"
amd_within_session <- sessionInfo()
capture.output(amd_within_session, file = file.path(amd_within_dir, "AMD_GSE38671_gene_level_sessionInfo_final_2026-11-05.txt"))
amd_final_checkpoint <- file.path(amd_within_dir, "AMD_GSE38671_gene_level_withinstudy_final_checkpoint_2026-11-05.RData")

# Preserve data, traceability, models, complete results, selected expression and plot objects.
# Keep the loaded precheck environment and redundant raw long-form text out of the final checkpoint.
save(amd_within_dir, amd_precheck_file, amd_analysis_date, amd_input_checksum, amd_within_session,
  amd_within_expr, amd_within_metadata, amd_mapping_audit, amd_feature_selection,
  amd_gene_representatives, amd_gene_expr, amd_gene_annotation, amd_features_per_gene,
  amd_design, amd_contrasts, amd_gene_fit, amd_gene_contrast_fit,
  amd_gene_results, amd_gene_degs, amd_gene_summary, amd_gene_ranked,
  amd_group_colours, amd_group_labels, amd_comparison_labels,
  amd_sa_plot_file, amd_sa_variance_outliers,
  amd_md_data, amd_md_labels, amd_md_plot, amd_md_file,
  amd_heat_samples, amd_heat_annotations, amd_heatmaps,
  amd_dot_key, amd_dot_genes, amd_dot_expr, amd_dot_data, amd_dot_labels, amd_dotplot, amd_dot_file,
  amd_figure_files, amd_final_checkpoint, file = amd_final_checkpoint)
stopifnot(file.exists(amd_final_checkpoint), file.info(amd_final_checkpoint)$size > 0)
print(amd_gene_summary)
print(c(Checkpoint_bytes = file.info(amd_final_checkpoint)$size))
