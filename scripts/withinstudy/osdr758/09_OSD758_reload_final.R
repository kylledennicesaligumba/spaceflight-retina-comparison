# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Restore the final OSD-758 analysis after its checkpoint has been saved.
# This reloads the existing results without repeating any statistical calculations.

# Locate and restore the completed checkpoint.
osd758_final_checkpoint <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758", "OSD758_withinstudy_final_checkpoint_2026-11-06.RData")
stopifnot(file.exists(osd758_final_checkpoint))
load(osd758_final_checkpoint)
osd758_within_dir <- file.path(portfolio_repo_dir, "results", "withinstudy", "osdr758")

# Display the restored result summary and full ranking count.
print(osd758_summary)
print(nrow(osd758_ranked))

# Check restored sample alignment, full ranking membership and score ordering.
osd758_reload_checks <- c(
  sample_order = identical(
    osd758_within_metadata$sample, colnames(osd758_within_dge$counts)
  ),
  design_order = identical(
    rownames(osd758_design), colnames(osd758_within_dge$counts)
  ),
  all_tested_genes = setequal(osd758_ranked$gene_id, osd758_results$gene_id),
  finite_scores = all(is.finite(osd758_ranked$rank_score)),
  decreasing_scores = !is.unsorted(-osd758_ranked$rank_score)
)
print(osd758_reload_checks)
stopifnot(all(osd758_reload_checks))
