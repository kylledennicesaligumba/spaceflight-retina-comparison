# Retinal flight-versus-terrestrial comparison: results and limits

Date: 7 November 2026 (Asia/Manila).

Status: coverage sensitivity and the conditional ranking-stability follow-up have been executed and verified. Formal biological ranking confidence intervals and a test separating the leading models have NOT been estimated. Conclusions therefore remain descriptive.

## Research question and preserved primary analysis

Compare the measured retinal responses of seven terrestrial disease/injury model studies (12 distinct contrasts) with the fixed OSD-758 mouse retinal flight0G-versus-HGC reference. SANS motivates the research; these are not retinal transcriptomes from astronauts diagnosed with SANS.

The primary comparison retains the same 12,127 harmonized mouse Entrez genes across all13 contrasts and 799 matched Reactome pathways across all13 columns. Every disease contrast remains separate. Original within-study scores, gene FDRs, mapping exclusions, pathway definitions and enrichment FDRs are unchanged.

## Complete-panel numerical overview

| Contrast | Gene rho: primary | Gene rho: pairwise | Pathway rho: primary | Pathway rank: primary | Rank range: group deletions |
| --- | --- | --- | --- | --- | --- |
| Glaucoma: 2 weeks vs naive | 0.243 | 0.233 | 0.302 |  1 | 1 to 2 |
| Autoimmune uveitis: EAU vs control | 0.145 | 0.136 | 0.273 |  2 | 1 to 4 |
| RP model: P23H vs WT | -0.031 | -0.025 | 0.260 |  3 | 1 to 6 |
| Glaucoma: 6 weeks vs naive | 0.129 | 0.123 | 0.179 |  4 | 2 to 5 |
| AMD-related: aged CfhKO vs WT | -0.039 | -0.045 | 0.174 |  5 | 3 to 7 |
| Glaucoma: 3 days vs naive | 0.044 | 0.039 | 0.120 |  6 | 4 to 6 |
| Diabetic retinal injury: diabetic vs normal | 0.115 | 0.096 | 0.050 |  7 | 5 to 7 |
| LHON-linked: COI vs WT | 0.022 | 0.026 | -0.053 |  8 | 8 to 9 |
| Light damage: LD vs control | -0.197 | -0.192 | -0.078 |  9 | 8 to 12 |
| LHON-linked: YARS2 vs WT | -0.019 | -0.023 | -0.106 | 10 | 9 to 12 |
| LHON-linked: YARS2 + COI vs WT | 0.002 | -0.004 | -0.119 | 11 | 10 to 12 |
| AMD-related: young CfhKO vs WT | 0.008 | -0.017 | -0.161 | 12 | 9 to 12 |

The deletion ranges in this table are sensitivity ranges, not confidence intervals. They concern pathway-profile correlations; no averaging of gene and pathway correlations into a disease score was performed.

## Coverage check: what was tested and what happened

Purpose: determine whether restricting every model to the same all-study common genes materially drives the observed gene-level resemblance. Each flight-model pair was rechecked on all eligible mapped genes available to that pair, without gene-FDR filtering, score changes or refitting. The primary common-gene analysis remains the main fair comparison; different pairwise universes are a sensitivity view.

Pairwise coverage spans 13,934 to 15,841 genes. Two-week glaucoma remains first: rho 0.242935 on 12,127 primary genes and 0.233039 on 15,841 pairwise genes. The top six gene-level positions are unchanged. Only 2 near-zero lower-ranking positions exchange order (young AMD and combined LHON); their sign changes are close to zero and do not alter the leading pattern.

## Pathway ranking-stability check: purpose and rules

Purpose: assess whether repeatedly represented, overlapping biological programs drive the observed pathway ranking. For each of the799 retained pathways, form a direct neighborhood containing that anchor and every pathway sharing at least half the genes of the smaller matched set. Delete the entire neighborhood together, then recompute all12 reference correlations using the SAME remaining pathway rows. No enrichment is rerun and no pathway is selected by FDR for this check. Identical deletion sets are evaluated once.

The overlap coefficient is established for representing gene-set relationships (Merico et al.,2010). Using its neighborhoods as deletion perturbations is a project-specific exploratory sensitivity adaptation, not a published confidence-interval estimator. Its rules were recorded before examining these new outcomes, after the primary ranking had already been reviewed.

There are 552 distinct deletion scenarios. They remove 1 to 193 pathways, leaving 606 to 798; all retain at least half of the original799 and are estimable. Two-week glaucoma ranks first in 545 scenarios and second in the remaining7. Uveitis is first in4 and RP in3. Glaucoma's rho ranges from0.261682 to0.429895. These scenario counts are NOT probabilities, bootstrap support or a confidence level: the neighborhoods overlap and are deterministic.

Examples: deleting the ECM-organization neighborhood removes15 related pathways and leaves two-week glaucoma first. Deleting the Neuronal System neighborhood removes55 pathways and puts RP first. The latter is consistent with the previously reviewed neuronal opposition in RP: removing that differentiating program can change the relative resemblance. These examples illustrate feature dependence, not an independent test of disease mechanisms.

## Biological interpretation already established

Two-week glaucoma has the highest observed gene-score and pathway-profile correlations. Matrix/collagen-positive and neuronal/synaptic-negative enrichment is shared with other models, particularly uveitis; it is not glaucoma-specific. RP shares matrix-related programs but opposes the neural program, while light damage has substantial opposing pathway enrichment.

The reviewed contributor analysis found174 unique shared leading-edge genes for flight-glaucoma2W; none meets the individual-gene FDR cutoff in BOTH studies (eight flight-only,166 neither). They are shared contributors, not174 replicated DEGs. Across selected contributors in the complete panel,10 both-gene-FDR contrast-records represent9 distinct genes. Pdia3 and Hspa8 oppose flight in their respective terrestrial comparisons. Col4a1 recurs as a higher-in-both DEG for RP and uveitis.

In the two illustrative programs,13 ECM and33 neuronal gene identities recur across the flight-glaucoma2W, flight-RP and flight-uveitis overlaps within the SAME pathway. These sets are not gene-FDR filtered. The reviewed alphabetical previews show matching positive ECM estimates and opposing RP neuronal estimates; do not claim all individual-gene signs/FDRs for unreviewed rows from a preview alone.

## What uncertainty remains unquantified

No biological-sample confidence interval, model-selection probability or formal test of the difference between disease correlations has been produced. Genes and overlapping pathways are not independent biological replicates. Treating their large row counts as the biological sample size would not answer the question of biological ranking precision.

The glaucoma model uses eight pooled libraries across four groups, one female and one male pooled library per group, with a common-sex-effect model (~sex+group) and nominal residual df3. An ordinary bootstrap within group-by-sex strata has only one library per stratum and cannot generate informative within-stratum resampling variation. Pooling also prevents recovering individual-animal variation. This does not make every possible model-based uncertainty method impossible; a separate design-aware resampling/modeling plan would be required and has not been implemented here.

The present scope is completed as an exploratory/descriptive comparison with explicit coverage and pathway-group robustness checks. A claim that glaucoma is statistically proven closer than uveitis/RP, a validated clinical SANS equivalent, or a causal gene mechanism is unsupported. If formal biological rank uncertainty remains a mandatory project objective, that specific objective remains unmet rather than being replaced by these sensitivity ranges.

## Defensible conclusion

Among the tested terrestrial retinal contrasts, two-week experimental glaucoma showed the highest observed transcriptomic resemblance to the OSD-758 flight-associated retinal response. The leading result persisted with expanded pairwise gene coverage and in most predefined pathway-neighborhood deletions, although uveitis or RP led after some program deletions. Shared matrix-related and neuronal expression programs explain components of the resemblance, with limited individual-gene confirmation in the glaucoma contributor set. These findings identify exploratory molecular parallels in mouse retina relevant to SANS research, without establishing a statistically separated disease winner or clinical equivalence to SANS.

## Saved outputs and reproducibility

Scripts24 and25 reproduce the new checks; script26 rebuilds this record from saved checkpoints. The short 00_resume_main_analysis.R restores the latest completed stage. Do not source the growing main script just to reload.

- Main_analysis_OSD758_gene_coverage_sensitivity.csv and corresponding PNG/PDF.
- Main_analysis_pairwise_gene_coverage_row_key.csv preserves the original row links.
- Main_analysis_pathway_rank_stability_summary.csv, scenarios.csv and profiles.csv, plus PNG/PDF and saved neighborhood memberships.
- Main_analysis_pathway_stability_ECM_neuronal_examples.csv contains exact example scores for all12 contrasts.
- Main_analysis_final_sensitivity_run_log_2026-11-07.txt records successful execution and checks.
- Each checkpoint includes method settings and input hashes; full numerical precision is retained in CSV/RDS outputs.

## Method sources and what they support

- [R cor documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html): Spearman calculation and tie ranking.
- [Merico et al.,2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2981572/): gene-set overlap coefficient/network representation; not this project's deletion protocol or sampling uncertainty.
- [Wu and Smyth,2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/): importance of inter-gene dependence in gene-set analysis; CAMERA itself was not substituted for the agreed fgsea analysis.
- [Goeman and Buhlmann,2007](https://pubmed.ncbi.nlm.nih.gov/17303618/): gene-versus-subject sampling and gene-set testing assumptions; not validation of the complete project workflow.

Preservation: source tables and all primary results were verified unchanged after execution. No synced project source material was modified.
