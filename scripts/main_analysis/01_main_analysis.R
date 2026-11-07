# Resolve repository paths before running the recorded analysis blocks.
if (!exists("portfolio_repo_dir", inherits = TRUE)) {
  portfolio_repo_dir <- normalizePath(Sys.getenv("RETINA_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
}

# Compare the saved OSD-758 retinal response with the seven terrestrial studies.
# Build the input manifest first; review the user's output before importing tables.

# Set the new output folder and the location of the completed within-study results.
main_analysis_dir <- Sys.getenv("MAIN_ANALYSIS_DIR", unset = file.path(portfolio_repo_dir, "results", "main_analysis"))
main_within_dir <- file.path(portfolio_repo_dir, "tables", "withinstudy")
dir.create(main_analysis_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(main_analysis_dir), dir.exists(main_within_dir))

# Register the flight reference and retain each terrestrial stage/genotype/age separately.
main_input_files <- data.frame(
  Study = c(
    "OSD758", rep("Glaucoma", 3), "RP", rep("LHON", 3),
    "DR", rep("AMD", 2), "Autoimmune_uveitis", "Light_damage"
  ),
  Comparison = c(
    "FLT_0G_vs_HGC",
    "3D_vs_Naive", "2W_vs_Naive", "6W_vs_Naive",
    "P23H_vs_WT",
    "YARS2_vs_WT", "COI_vs_WT", "YARS2_COI_vs_WT",
    "Diabetic_vs_Normal",
    "Young_CfhKO_vs_WT", "Aged_CfhKO_vs_WT",
    "EAU_vs_CON", "LD_vs_CTRL"
  ),
  Role = c("Reference", rep("Terrestrial", 12)),
  File = file.path(
    main_within_dir,
    c(
      "osdr758/OSD758_FLT_0G_vs_HGC_ranked_all_genes.csv",
      "glaucoma/Glaucoma_GSE241782_3D_vs_Naive_ranked_all_genes.csv",
      "glaucoma/Glaucoma_GSE241782_2W_vs_Naive_ranked_all_genes.csv",
      "glaucoma/Glaucoma_GSE241782_6W_vs_Naive_ranked_all_genes.csv",
      "retinitis_pigmentosa/RP_GSE156533_P23H_vs_WT_ranked_all_genes.csv",
      "lhon/LHON_GSE304707_YARS2_vs_WT_ranked_all_genes.csv",
      "lhon/LHON_GSE304707_COI_vs_WT_ranked_all_genes.csv",
      "lhon/LHON_GSE304707_YARS2_COI_vs_WT_ranked_all_genes.csv",
      "diabetic_retinopathy/DR_GSE87433_Diabetic_vs_Normal_ranked_all_genes.csv",
      "amd/AMD_GSE38671_gene_level_Young_CfhKO_vs_WT_ranked_all_genes.csv",
      "amd/AMD_GSE38671_gene_level_Aged_CfhKO_vs_WT_ranked_all_genes.csv",
      "autoimmune_uveitis/Autoimmune_Uveitis_GSE284505_EAU_vs_CON_ranked_all_genes.csv",
      "light_damage/Light_Damage_GSE314043_LD_vs_CTRL_ranked_all_genes.csv"
    )
  ),
  stringsAsFactors = FALSE
)

# Confirm that every complete ranked-result file can be found from this R session.
main_input_files$Exists <- file.exists(main_input_files$File)
print(main_input_files[, c("Study", "Comparison", "Exists")], row.names = FALSE)
stopifnot(all(main_input_files$Exists))

# Save an explicit source-file record in the new main-analysis folder.
main_manifest_file <- file.path(main_analysis_dir, "Main_analysis_input_manifest.csv")
write.csv(main_input_files, file = main_manifest_file, row.names = FALSE)
print(file.exists(main_manifest_file))

# User confirmed all 13 source files were found and the manifest was saved.

# Import all full ranked tables and review columns plus source-identifier examples.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "02_import_ranked_tables.R"))

# User reviewed the imported dimensions and identifier previews; inventory save is TRUE.
# Use mouse Entrez Gene IDs as the common target and audit mappings before matching genes.

# Install/load the mouse annotation database and record its versions and source metadata.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "03_prepare_gene_mapping.R"))

# User confirmed org.Mm.eg.db 3.22.0, all four key types and both setup-file checks.

# Audit all candidate Entrez mappings for the fixed OSD-758 reference before filtering.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "04_map_OSD758_reference.R"))

# User confirmed 16,615 one-match rows, 104 multiple-match rows and 4,616 no-match rows.
# Among the one-match rows, 29 share a target; both mapping CSV files were saved.

# Review which mapping groups contain significant flight-response genes before eligibility decisions.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "05_review_OSD758_mapping.R"))

# User confirmed 319 DEGs in one-match/no-repeat rows, 51 unmatched and three with
# multiple candidates. Shared-target rows contain no DEGs; both review files were saved.
# Keep unresolved cases documented while auditing the remaining source datasets.

# Apply the same Ensembl-to-Entrez audit to glaucoma stages, RP and light damage.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "06_map_terrestrial_Ensembl.R"))

# User confirmed all five audits and six file saves. No source IDs were missing.
# Glaucoma per stage: 17,653 one-match / 109 multiple / 5,539 no-match / 33 shared-target rows.
# RP: 15,208 / 81 / 2,062 / 6; light damage: 15,021 / 81 / 1,750 / 6.

# Audit the supplied DR/AMD Entrez IDs while keeping AMD's representative features fixed.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "07_audit_existing_Entrez.R"))

# User confirmed DR: 24,685 IDs in the snapshot and 288 absent; AMD per age: 19,743 / 82.
# No malformed, missing or repeated native IDs were found; all four audit files were saved.

# Audit the original LHON/uveitis labels, recording exact-symbol and alias lookup routes.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "08_map_gene_name_labels.R"))

# User confirmed all six saves. LHON per contrast: 17,113 one-match / 9 multiple /
# 738 no-match / 29 shared-target rows; routes: 16,683 official-symbol and 439 alias matches.
# Uveitis: 16,579 / 5 / 1,603 / 24; routes: 16,301 official-symbol, 239 alias and 44 Ensembl.
# No source IDs were missing. Unresolved cases remain documented for follow-up.

# Assess provisional shared-gene coverage with explicit unique-ID eligibility rules.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "09_build_common_gene_coverage.R"))

# User confirmed 12,127 candidate common genes and four saved outputs, then clarified
# the mapping exclusions and pairwise-versus-all-13 distinction before requesting continuation.

# Align the accepted common genes across all response tables and verify their source-row links.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "10_align_common_responses.R"))

# User confirmed all three matrices are 12,127 x 13, all eight alignment/value checks
# are TRUE and all five response-matrix/checkpoint outputs were saved.

# Calculate descriptive common-gene Spearman scores against the fixed OSD-758 reference.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "11_gene_level_Spearman.R"))

# User confirmed all 12 scores on 12,127 genes and three saved outputs.
# Highest observed rho: 2-week glaucoma 0.2429348; most negative: light damage -0.1967997.
# These are descriptive ordering agreements; magnitude and robustness need interpretation.

# Draw and save all observed scores on the full correlation scale with numeric labels.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "12_plot_gene_Spearman.R"))

# User reviewed the Spearman figure and requests gene annotation to start the pathway layer.
# The planned pairwise-coverage sensitivity remains pending for the gene-level results.

# Retrieve gene symbols and full names for review, keeping the accepted common IDs fixed.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "13_annotate_common_genes.R"))

# User confirmed 12,127 common genes, complete symbols/full names, no multiple rows,
# and a saved raw lookup. The select 1:1 message is informational.

# Align the reviewed labels to the common-response matrices and save the annotation key.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "14_save_common_annotation.R"))

# User confirmed 12,127 x 3, all four annotation checks TRUE and both files saved.

# Prepare and record the mouse Reactome pathway collection before checking common-gene coverage.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "15_prepare_Reactome_resource.R"))

# User confirmed reactome.db 1.95.0, 1,822 mouse pathways, 1,815 with Entrez
# members, 8,780 distinct resource genes and all four saved outputs.

# Match pathways to the common genes and retain sets with 15-500 matched members.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "16_match_Reactome_common_genes.R"))

# User confirmed 6,130 common genes have any Reactome membership, 5,997 have none,
# and 5,956 appear in eligible pathways. All 799 eligible pathways, five checks,
# and four saved files were reviewed; all 12,127 genes remain in the ranking universe.

# Run the same preranked enrichment in every contrast and save progress after each one.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "17_run_Reactome_fgsea.R"))

# User confirmed 799 estimated pathways per contrast, no warnings and seven saved files.
# OSD-758 has 44 enriched pathways at within-contrast FDR < 0.05 (19 positive NES,
# 25 negative NES). The all-pathway NES profile comparison remains a separate next block.

# Align all pathway profiles, calculate descriptive reference correlations and save paired evidence flags.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "18_compare_Reactome_profiles.R"))

# User confirmed 799 x 13 matrices, six TRUE checks and seven saved outputs.
# Glaucoma 2W has the highest observed pathway rho (0.302077), followed by
# uveitis (0.273064) and RP (0.260153). These are descriptive profile agreements.
# Both-enriched same/opposite counts: glaucoma 2W 18/0; uveitis 11/1; RP 21/3.
# Review actual matched pathway identities and leading-edge genes before mechanism claims.

# Display the 2-week glaucoma shared pathway names and save name lists for all contrasts.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "19_inspect_shared_pathway_names.R"))

# User reviewed all 18 shared names: nine positive ECM/collagen/cell-matrix sets
# and nine negative sets including neuronal/synaptic, cell-cycle and protein-processing themes.
# Four name-list files were saved. Annotated leading-edge gene comparison remains next.

# Compare and annotate leading-edge genes for all 18 reviewed shared glaucoma pathways.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "20_compare_shared_leading_edge_genes.R"))

# User confirmed 174 unique shared leading-edge genes: 41 higher and 133 lower
# in both contrasts, none with opposite estimated logFC, five TRUE checks and four saves.
# These are contributors, not automatically DEGs; all 20 alphabetical preview genes
# have gene FDR > 0.05 in both studies. Full individual-gene support remains to review.

# Compare every flight-enriched pathway across all terrestrial contrasts using saved NES/FDR values.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "21_compare_programs_full_panel.R"))

# Review the positive and negative pathway heatmaps and the full-panel summary together.
# Original all-pathway correlations remain unchanged; this display retains unstarred NES cells.

# User confirmed all five checks and seven saved files, then requested stopping before interpretation.
# Resume with 00_resume_main_analysis.R and review both full-panel PNGs together next time.

# User reviewed both full-panel figures on 7 November: shared and opposing named programs differ across models.
# Extend the same within-pathway contributor comparison to all 12 contrasts and report original gene FDR support.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "22_compare_leading_edges_full_panel.R"))

# Review the contributor and gene-FDR summaries before interpreting new shared-gene identities.

# User confirmed all eight checks/saves. Glaucoma2W's174 contributors have gene-FDR support:
# zero both, eight flight-only, zero terrestrial-only,166 neither. They are not174 shared DEGs.
# Both-study gene-FDR contributor counts: RP2, uveitis5, lightdamage3; inspect names/logFC next.

# Inspect saved gene names, original fold changes/FDRs and same-pathway overlap across models.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "23_inspect_contributor_gene_identities.R"))

# Begin with both-gene-FDR contributors across all12 contrasts and the reviewed ECM/neuronal programs.
# Common-program examples retain non-DEGs; their actual identities remain pending user output.

# User confirmed five checks/six saves and reviewed10 both-FDR contrast records (9 unique genes).
# Same-pathway intersections: ECM13 (pairwise41/28/29), neuronal33 (58/94/70).
# Only10 alphabetical previews per program shown; retain original gene-sign/FDR distinctions.

# On 7 November the user explicitly authorized the assistant to execute the remaining analyses.
# Check expanded pairwise gene coverage while preserving original mappings and signed scores.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "24_gene_coverage_sensitivity.R"))

# Assess conditional pathway-ranking stability by deleting overlapping program neighborhoods together.
# These ranges are sensitivity results, not biological confidence intervals or winner probabilities.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "25_pathway_rank_stability.R"))

# Consolidate verified results and explicit uncertainty limits without changing primary outcomes.
source(file.path(portfolio_repo_dir, "scripts", "main_analysis", "26_finalize_analysis_record.R"))
