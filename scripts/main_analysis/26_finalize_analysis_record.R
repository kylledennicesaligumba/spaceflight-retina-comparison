# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Consolidate completed sensitivity results and their limits into a restartable analysis record.
# Read saved results only; do not refit models, repeat enrichment or change the primary analysis.

# Format full-panel results for the written record without altering their stored precision.
main_final_markdown_table <- function(x) {
  headers <- names(x)
  rows <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(headers)), collapse = " | "), " |"), rows)
}

# Restore the completed identity view and both final sensitivity checkpoints.
if (!exists("main_analysis_dir")) main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_final_inputs <- file.path(main_analysis_dir, c(
  "Main_analysis_gene_coverage_sensitivity_checkpoint_2026-11-07.rds",
  "Main_analysis_pathway_rank_stability_checkpoint_2026-11-07.rds",
  "Main_analysis_contributor_identity_review_checkpoint_2026-11-07.rds"))
stopifnot(all(file.exists(main_final_inputs)))
main_final_cov <- readRDS(main_final_inputs[1])
main_final_stab <- readRDS(main_final_inputs[2])
main_final_identity <- readRDS(main_final_inputs[3])
stopifnot(all(main_final_cov$checks), all(main_final_stab$checks), all(main_final_identity$checks),
  identical(tools::md5sum(main_final_cov$input_files), main_final_cov$input_md5),
  identical(tools::md5sum(names(main_final_cov$source_input_md5)), main_final_cov$source_input_md5),
  identical(tools::md5sum(main_final_stab$input_files), main_final_stab$input_md5),
  identical(tools::md5sum(main_final_identity$input_file), main_final_identity$input_md5))
main_final_genes <- main_final_cov$results
main_final_pathways <- main_final_stab$summary
main_final_keys <- paste(main_final_pathways$Study, main_final_pathways$Comparison, sep = "__")
g <- match(main_final_keys, paste(main_final_genes$Study, main_final_genes$Comparison, sep = "__"))
stopifnot(!anyNA(g))
main_final_combined <- data.frame(Study = main_final_pathways$Study, Comparison = main_final_pathways$Comparison,
  Gene_rho_primary = main_final_genes$Primary_rho[g], Gene_rho_pairwise = main_final_genes$Pairwise_rho[g],
  Pathway_rho_primary = main_final_pathways$Baseline_rho, Pathway_rank_primary = main_final_pathways$Baseline_rank,
  Pathway_rank_best = main_final_pathways$Best_rank, Pathway_rank_worst = main_final_pathways$Worst_rank,
  Top_deletion_scenarios = main_final_pathways$Top_rank_scenarios)
main_final_combined <- main_final_combined[order(main_final_combined$Pathway_rank_primary), ]
rownames(main_final_combined) <- NULL

# Show the actual effect of deleting the two previously discussed program neighborhoods.
main_final_example_names <- c("Extracellular matrix organization", "Neuronal System")
main_final_examples <- list()
for (name in main_final_example_names) {
  ids <- which(vapply(strsplit(main_final_stab$scenarios$Anchor_pathway_names, ";", fixed = TRUE),
    function(x) name %in% x, logical(1)))
  stopifnot(length(ids) == 1L)
  scenario <- main_final_stab$scenarios[ids, ]
  rows <- main_final_stab$profiles[main_final_stab$profiles$Scenario == scenario$Scenario, ]
  rows$Deleted_neighborhood <- name
  rows$Removed_pathways <- scenario$Removed_pathways
  main_final_examples[[name]] <- rows[, c("Deleted_neighborhood", "Removed_pathways", "Retained_pathways",
    "Study", "Comparison", "Spearman_rho", "Observed_rank")]
}
main_final_examples <- do.call(rbind, main_final_examples)
rownames(main_final_examples) <- NULL
main_final_nscenarios <- length(main_final_stab$principal_scenarios)
main_final_top_counts <- main_final_pathways[main_final_pathways$Top_rank_scenarios > 0,
  c("Study", "Comparison", "Top_rank_scenarios", "Best_rank", "Worst_rank"), drop = FALSE]
main_final_rank_changes <- sum(main_final_genes$Primary_rank != main_final_genes$Pairwise_rank)
main_final_glaucoma <- main_final_genes[main_final_genes$Study == "Glaucoma" &
  main_final_genes$Comparison == "2W_vs_Naive", ]
main_final_glaucoma_path <- main_final_pathways[main_final_pathways$Study == "Glaucoma" &
  main_final_pathways$Comparison == "2W_vs_Naive", ]
main_final_table <- data.frame(Contrast = unname(main_final_cov$labels[paste(main_final_combined$Study,
    main_final_combined$Comparison, sep = "__")]),
  `Gene rho: primary` = sprintf("%.3f", main_final_combined$Gene_rho_primary),
  `Gene rho: pairwise` = sprintf("%.3f", main_final_combined$Gene_rho_pairwise),
  `Pathway rho: primary` = sprintf("%.3f", main_final_combined$Pathway_rho_primary),
  `Pathway rank: primary` = main_final_combined$Pathway_rank_primary,
  `Rank range: group deletions` = paste0(main_final_combined$Pathway_rank_best, " to ", main_final_combined$Pathway_rank_worst),
  check.names = FALSE)
main_final_files <- file.path(main_analysis_dir, c("Main_analysis_results_and_limits_2026-11-07.md",
  "Main_analysis_final_comparison_summary.csv", "Main_analysis_pathway_stability_ECM_neuronal_examples.csv",
  "Main_analysis_final_review_checkpoint_2026-11-07.rds", "Main_analysis_final_review_sessionInfo_2026-11-07.txt"))
main_final_report <- c(
  "# Retinal flight-versus-terrestrial comparison: results and limits",
  "", "Date: 7 November 2026 (Asia/Manila).",
  "", "Status: coverage sensitivity and the conditional ranking-stability follow-up have been executed and verified. Formal biological ranking confidence intervals and a test separating the leading models have NOT been estimated. Conclusions therefore remain descriptive.",
  "", "## Research question and preserved primary analysis",
  "", "Compare the measured retinal responses of seven terrestrial disease/injury model studies (12 distinct contrasts) with the fixed OSD-758 mouse retinal flight0G-versus-HGC reference. SANS motivates the research; these are not retinal transcriptomes from astronauts diagnosed with SANS.",
  "", sprintf("The primary comparison retains the same %s harmonized mouse Entrez genes across all13 contrasts and %s matched Reactome pathways across all13 columns. Every disease contrast remains separate. Original within-study scores, gene FDRs, mapping exclusions, pathway definitions and enrichment FDRs are unchanged.",
    format(length(main_final_cov$common_ids), big.mark = ","), format(length(main_final_stab$pathway_ids), big.mark = ",")),
  "", "## Complete-panel numerical overview", "", main_final_markdown_table(main_final_table),
  "", "The deletion ranges in this table are sensitivity ranges, not confidence intervals. They concern pathway-profile correlations; no averaging of gene and pathway correlations into a disease score was performed.",
  "", "## Coverage check: what was tested and what happened",
  "", "Purpose: determine whether restricting every model to the same all-study common genes materially drives the observed gene-level resemblance. Each flight-model pair was rechecked on all eligible mapped genes available to that pair, without gene-FDR filtering, score changes or refitting. The primary common-gene analysis remains the main fair comparison; different pairwise universes are a sensitivity view.",
  "", sprintf("Pairwise coverage spans %s to %s genes. Two-week glaucoma remains first: rho %.6f on %s primary genes and %.6f on %s pairwise genes. The top six gene-level positions are unchanged. Only %s near-zero lower-ranking positions exchange order (young AMD and combined LHON); their sign changes are close to zero and do not alter the leading pattern.",
    format(min(main_final_genes$Pairwise_genes), big.mark = ","), format(max(main_final_genes$Pairwise_genes), big.mark = ","),
    main_final_glaucoma$Primary_rho, format(main_final_glaucoma$Primary_genes, big.mark = ","),
    main_final_glaucoma$Pairwise_rho, format(main_final_glaucoma$Pairwise_genes, big.mark = ","), main_final_rank_changes),
  "", "## Pathway ranking-stability check: purpose and rules",
  "", "Purpose: assess whether repeatedly represented, overlapping biological programs drive the observed pathway ranking. For each of the799 retained pathways, form a direct neighborhood containing that anchor and every pathway sharing at least half the genes of the smaller matched set. Delete the entire neighborhood together, then recompute all12 reference correlations using the SAME remaining pathway rows. No enrichment is rerun and no pathway is selected by FDR for this check. Identical deletion sets are evaluated once.",
  "", "The overlap coefficient is established for representing gene-set relationships (Merico et al.,2010). Using its neighborhoods as deletion perturbations is a project-specific exploratory sensitivity adaptation, not a published confidence-interval estimator. Its rules were recorded before examining these new outcomes, after the primary ranking had already been reviewed.",
  "", sprintf("There are %s distinct deletion scenarios. They remove %s to %s pathways, leaving %s to %s; all retain at least half of the original799 and are estimable. Two-week glaucoma ranks first in %s scenarios and second in the remaining%s. Uveitis is first in4 and RP in3. Glaucoma's rho ranges from%.6f to%.6f. These scenario counts are NOT probabilities, bootstrap support or a confidence level: the neighborhoods overlap and are deterministic.",
    main_final_nscenarios, min(main_final_stab$scenarios$Removed_pathways), max(main_final_stab$scenarios$Removed_pathways),
    min(main_final_stab$scenarios$Retained_pathways), max(main_final_stab$scenarios$Retained_pathways),
    main_final_glaucoma_path$Top_rank_scenarios, main_final_nscenarios - main_final_glaucoma_path$Top_rank_scenarios,
    main_final_glaucoma_path$Minimum_rho, main_final_glaucoma_path$Maximum_rho),
  "", "Examples: deleting the ECM-organization neighborhood removes15 related pathways and leaves two-week glaucoma first. Deleting the Neuronal System neighborhood removes55 pathways and puts RP first. The latter is consistent with the previously reviewed neuronal opposition in RP: removing that differentiating program can change the relative resemblance. These examples illustrate feature dependence, not an independent test of disease mechanisms.",
  "", "## Biological interpretation already established",
  "", "Two-week glaucoma has the highest observed gene-score and pathway-profile correlations. Matrix/collagen-positive and neuronal/synaptic-negative enrichment is shared with other models, particularly uveitis; it is not glaucoma-specific. RP shares matrix-related programs but opposes the neural program, while light damage has substantial opposing pathway enrichment.",
  "", "The reviewed contributor analysis found174 unique shared leading-edge genes for flight-glaucoma2W; none meets the individual-gene FDR cutoff in BOTH studies (eight flight-only,166 neither). They are shared contributors, not174 replicated DEGs. Across selected contributors in the complete panel,10 both-gene-FDR contrast-records represent9 distinct genes. Pdia3 and Hspa8 oppose flight in their respective terrestrial comparisons. Col4a1 recurs as a higher-in-both DEG for RP and uveitis.",
  "", "In the two illustrative programs,13 ECM and33 neuronal gene identities recur across the flight-glaucoma2W, flight-RP and flight-uveitis overlaps within the SAME pathway. These sets are not gene-FDR filtered. The reviewed alphabetical previews show matching positive ECM estimates and opposing RP neuronal estimates; do not claim all individual-gene signs/FDRs for unreviewed rows from a preview alone.",
  "", "## What uncertainty remains unquantified",
  "", "No biological-sample confidence interval, model-selection probability or formal test of the difference between disease correlations has been produced. Genes and overlapping pathways are not independent biological replicates. Treating their large row counts as the biological sample size would not answer the question of biological ranking precision.",
  "", "The glaucoma model uses eight pooled libraries across four groups, one female and one male pooled library per group, with a common-sex-effect model (~sex+group) and nominal residual df3. An ordinary bootstrap within group-by-sex strata has only one library per stratum and cannot generate informative within-stratum resampling variation. Pooling also prevents recovering individual-animal variation. This does not make every possible model-based uncertainty method impossible; a separate design-aware resampling/modeling plan would be required and has not been implemented here.",
  "", "The present scope is completed as an exploratory/descriptive comparison with explicit coverage and pathway-group robustness checks. A claim that glaucoma is statistically proven closer than uveitis/RP, a validated clinical SANS equivalent, or a causal gene mechanism is unsupported. If formal biological rank uncertainty remains a mandatory project objective, that specific objective remains unmet rather than being replaced by these sensitivity ranges.",
  "", "## Defensible conclusion",
  "", "Among the tested terrestrial retinal contrasts, two-week experimental glaucoma showed the highest observed transcriptomic resemblance to the OSD-758 flight-associated retinal response. The leading result persisted with expanded pairwise gene coverage and in most predefined pathway-neighborhood deletions, although uveitis or RP led after some program deletions. Shared matrix-related and neuronal expression programs explain components of the resemblance, with limited individual-gene confirmation in the glaucoma contributor set. These findings identify exploratory molecular parallels in mouse retina relevant to SANS research, without establishing a statistically separated disease winner or clinical equivalence to SANS.",
  "", "## Saved outputs and reproducibility",
  "", "Scripts24 and25 reproduce the new checks; script26 rebuilds this record from saved checkpoints. The short 00_resume_main_analysis.R restores the latest completed stage. Do not source the growing main script just to reload.",
  "", "- Main_analysis_OSD758_gene_coverage_sensitivity.csv and corresponding PNG/PDF.",
  "- Main_analysis_pairwise_gene_coverage_row_key.csv preserves the original row links.",
  "- Main_analysis_pathway_rank_stability_summary.csv, scenarios.csv and profiles.csv, plus PNG/PDF and saved neighborhood memberships.",
  "- Main_analysis_pathway_stability_ECM_neuronal_examples.csv contains exact example scores for all12 contrasts.",
  "- Main_analysis_final_sensitivity_run_log_2026-11-07.txt records successful execution and checks.",
  "- Each checkpoint includes method settings and input hashes; full numerical precision is retained in CSV/RDS outputs.",
  "", "## Method sources and what they support",
  "", "- [R cor documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html): Spearman calculation and tie ranking.",
  "- [Merico et al.,2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2981572/): gene-set overlap coefficient/network representation; not this project's deletion protocol or sampling uncertainty.",
  "- [Wu and Smyth,2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/): importance of inter-gene dependence in gene-set analysis; CAMERA itself was not substituted for the agreed fgsea analysis.",
  "- [Goeman and Buhlmann,2007](https://pubmed.ncbi.nlm.nih.gov/17303618/): gene-versus-subject sampling and gene-set testing assumptions; not validation of the complete project workflow.",
  "", "Preservation: source tables and all primary results were verified unchanged after execution. No synced project source material was modified."
)
writeLines(main_final_report, main_final_files[1])
utils::write.csv(main_final_combined, main_final_files[2], row.names = FALSE)
utils::write.csv(main_final_examples, main_final_files[3], row.names = FALSE)
main_final_checks <- c(Coverage_passed = all(main_final_cov$checks), Stability_passed = all(main_final_stab$checks),
  Primary_sources_unchanged = identical(tools::md5sum(names(main_final_cov$source_input_md5)), main_final_cov$source_input_md5),
  All_contrasts_retained = nrow(main_final_combined) == 12L,
  Identity_review_preserved = all(main_final_identity$checks),
  Explicit_uncertainty_limits = any(grepl("Formal biological ranking confidence intervals", main_final_report, fixed = TRUE)))
stopifnot(all(main_final_checks))
saveRDS(list(comparison_summary = main_final_combined, coverage_summary = main_final_genes,
  pathway_stability = main_final_pathways, key_program_examples = main_final_examples,
  top_scenario_counts = main_final_top_counts, checks = main_final_checks,
  inference_status = "Descriptive analysis and sensitivity follow-up complete; formal biological ranking uncertainty NOT estimated",
  report_file = main_final_files[1], input_files = main_final_inputs, input_md5 = tools::md5sum(main_final_inputs)), main_final_files[4])
utils::capture.output(utils::sessionInfo(), file = main_final_files[5])
print(main_final_checks)
print(file.exists(main_final_files))
