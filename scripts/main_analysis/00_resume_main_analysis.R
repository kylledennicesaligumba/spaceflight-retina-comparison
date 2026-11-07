# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the latest saved full-panel review without rerunning the analysis.
main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_final_resume_file <- file.path(main_analysis_dir,
  "Main_analysis_final_review_checkpoint_2026-11-07.rds")
main_identity_resume_file <- file.path(main_analysis_dir,
  "Main_analysis_contributor_identity_review_checkpoint_2026-11-07.rds")
main_all_le_resume_file <- file.path(main_analysis_dir,
  "Main_analysis_leading_edge_full_panel_checkpoint_2026-11-07.rds")
main_panel_resume_file <- file.path(main_analysis_dir,
  "Main_analysis_programs_full_panel_checkpoint_2026-11-06.rds")

# Prefer the completed final sensitivity record, then the earlier reviewed stages.
if (file.exists(main_final_resume_file)) {
  main_final_review <- readRDS(main_final_resume_file)
  stopifnot(length(main_final_review$checks) > 0L, all(main_final_review$checks))
  print(main_final_review$comparison_summary, digits = 4, row.names = FALSE)
  print(main_final_review$checks)
  cat(main_final_review$inference_status, "\n")
  cat("Written record:", main_final_review$report_file, "\n")
} else if (file.exists(main_identity_resume_file)) {
  main_identity_review <- readRDS(main_identity_resume_file)
  stopifnot(length(main_identity_review$checks) > 0L, all(main_identity_review$checks))
  print(main_identity_review$program_summary, row.names = FALSE)
  print(main_identity_review$both_gene_FDR[, c("Study", "Comparison", "Gene_symbol",
    "Gene_logFC_pattern")], row.names = FALSE)
  print(main_identity_review$checks)
} else if (file.exists(main_all_le_resume_file)) {
  main_all_le_review <- readRDS(main_all_le_resume_file)
  stopifnot(length(main_all_le_review$checks) > 0L, all(main_all_le_review$checks))
  print(main_all_le_review$summary, row.names = FALSE)
  print(main_all_le_review$gene_FDR_support[
    main_all_le_review$gene_FDR_support$Pathway_group == "All selected pathways",
    c("Study", "Comparison", "Shared_genes", "Both_gene_FDR005", "Flight_only_gene_FDR005",
      "Terrestrial_only_gene_FDR005", "Neither_gene_FDR005")], row.names = FALSE)
  print(main_all_le_review$checks)
} else {
  stopifnot(file.exists(main_panel_resume_file))
  main_panel_review <- readRDS(main_panel_resume_file)
  stopifnot(length(main_panel_review$checks) > 0L, all(main_panel_review$checks))
  print(main_panel_review$summary, row.names = FALSE)
  print(main_panel_review$checks)
}

# Locate the latest completed figures without recalculating them.
main_panel_review_figures <- if (file.exists(main_final_resume_file)) {
  file.path(main_analysis_dir, c("Main_analysis_OSD758_gene_coverage_sensitivity.png",
    "Main_analysis_OSD758_pathway_rank_stability.png"))
} else {
  file.path(main_analysis_dir, c("Main_analysis_OSD758_positive_programs_full_panel.png",
    "Main_analysis_OSD758_negative_programs_full_panel.png"))
}
print(file.exists(main_panel_review_figures))
