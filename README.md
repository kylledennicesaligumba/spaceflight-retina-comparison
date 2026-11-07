# Comparative transcriptomics of mouse retinal responses to spaceflight, disease, and injury

Maintained by [Kylle Dennice Saligumba](https://github.com/kylledennicesaligumba) · Completed 7 November 2026

This exploratory study compares the mouse retinal response to spaceflight with **12 disease and injury contrasts from seven terrestrial studies** to identify which models show similar gene-expression changes and where their responses differ. The comparison uses 12,127 shared mouse genes and 799 matched Reactome pathways.

## Main findings

**Two-week experimental glaucoma showed the closest observed resemblance to the flight response**, with the highest rank correlation across both gene scores (Spearman ρ = 0.243) and pathway profiles (ρ = 0.302). [Full comparison table](tables/main_analysis/Main_analysis_final_comparison_summary.csv)

### Gene-level resemblance across all 12 contrasts

![Rank correlations between the flight gene scores and each of the 12 terrestrial contrasts, with two-week glaucoma highest and light damage lowest](figures/main_analysis/Main_analysis_OSD758_gene_Spearman_common.png)

Each point compares a terrestrial contrast with the same flight reference using the same 12,127 genes. Positive values indicate similar gene-score rankings and negative values indicate opposing rankings.

### Shared and opposing pathway responses

Flight, two-week glaucoma and autoimmune uveitis shared positive enrichment in extracellular matrix and collagen programs and negative enrichment in neuronal and synaptic programs. The P23H retinitis pigmentosa model shared the matrix pattern and showed positive enrichment in the neuronal and synaptic pathways that were negatively enriched in flight.

![Positive flight pathway heatmap showing shared matrix and collagen enrichment across several models](figures/main_analysis/Main_analysis_OSD758_positive_programs_full_panel.png)

![Negative flight pathway heatmap showing matching neuronal enrichment directions in glaucoma and uveitis and opposing directions in retinitis pigmentosa and light damage](figures/main_analysis/Main_analysis_OSD758_negative_programs_full_panel.png)

The first column in both heatmaps is the flight reference. Orange indicates positive enrichment and blue indicates negative enrichment. Stars mark pathways with a false discovery rate below 5% in that contrast.

### How stable was the ranking?

Two-week glaucoma also ranked first in the gene comparison when each model was compared on all eligible genes it shared with flight. It ranked first in 545 of 552 pathway sensitivity checks that removed groups of overlapping pathways, and second in the other seven. These checks measure sensitivity to which genes and pathways are included. Biological sampling uncertainty in the ranking remains unquantified. [Gene coverage check](figures/main_analysis/Main_analysis_OSD758_gene_coverage_sensitivity.png) · [Pathway rank stability](figures/main_analysis/Main_analysis_OSD758_pathway_rank_stability.png) · [Full results and limitations](docs/main_analysis/records/Main_analysis_results_and_limits_2026-11-07.md)

## Study design and workflow

Spaceflight-associated neuro-ocular syndrome (SANS) motivates the research. The fixed reference compares mice flown at 0 g with habitat ground controls. The researchers collected the retinal tissue after the mice returned to Earth, so the measured response includes recovery exposures. [NASA OSD-758](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758) · [Kremsky et al. (2024)](https://doi.org/10.3390/ijms252212041)

| Stage | Start here |
| --- | --- |
| Dataset selection and quality checks | [Dataset panel and measurement differences](docs/datasets.md) · [Terrestrial prechecks](#completed-prechecks) · [OSD-758 precheck](docs/osdr758/README.md) |
| Differential expression within each study | [Results for all 12 terrestrial contrasts](docs/withinstudy/README.md#recorded-comparison-results) · [Flight reference results](docs/withinstudy/osdr758.md) · [Scripts](scripts/withinstudy/) |
| Flight-to-terrestrial comparison | [Final analysis and reproduction instructions](docs/main_analysis/README.md) · [Scripts](scripts/main_analysis/) · [Tables](tables/main_analysis/) · [Figures](figures/main_analysis/) |

## Within-study differential expression

<details>
<summary>Study models, gene rankings and execution records</summary>

The [within-study results](docs/withinstudy/README.md) now include study-specific models, complete tested-gene tables, FDR-qualified gene lists, full signed rankings and saved result figures. Glaucoma and RP account for sex; LHON retains three genotype-versus-WT contrasts; DR uses the untreated six-versus-six comparison; AMD has two age-matched gene-level contrasts; uveitis and light damage retain their selected disease/injury comparisons.

RNA-seq models reuse the filtered counts and TMM factors with robust edgeR quasi-likelihood testing. The microarrays use robust limma-trend models. AMD's original array features are mapped conservatively to 19,825 recorded Entrez IDs using fixed maximum-mean representatives, with saved mapping/selection audits. The earlier feature-level results remain identifiable as history. Gene-level FDR calls use BH < 0.05, without an extra fold-change threshold; all tested rows remain in the rankings.

The [overview of 12 contrasts](docs/withinstudy/README.md#recorded-comparison-results) links each dataset's results and figures. The final comparison uses harmonized ranks and pathway profiles. The [OSD-758 result](docs/withinstudy/osdr758.md) supplies the flight reference. [Execution instructions](docs/withinstudy/running-the-analysis.md) describe the required local checkpoints; RP's final checkpoint/session export was not present, while its completed tables, rankings and figures are included.

</details>

## Completed prechecks

<details>
<summary>Quality-check decisions, dataset distinctions and figures</summary>

All seven processed-data prechecks are completed. Each includes input integrity, sample/feature annotation alignment, expression-distribution checks and reviewed MDS/PCA figures. RNA-seq workflows use count filtering and TMM. The microarrays follow distinct branches: **DR uses supplied processed values without additional log2, CPM or TMM; AMD uses a single log2 QC copy of supplied intensities, with no CPM, TMM, extra normalization or expression-row filtering.**

| Dataset | Completed | Samples used for QC | Feature rows used for QC | Record |
| --- | --- | --- | --- | --- |
| Experimental glaucoma, GSE241782 | 28 October 2026 | All eight selected retinal libraries | 23,301 of 55,487 | [Glaucoma report](docs/glaucoma-precheck.md) |
| Retinitis pigmentosa, GSE156533 | 29 October 2026 | All seven libraries: four P23H and three WT | 17,351 of 54,532 | [RP report](docs/retinitis-pigmentosa-precheck.md) |
| LHON-linked model, GSE304707 | 30 October 2026 | All 12 libraries: three per genotype | 17,860 of 49,714 | [LHON report](docs/lhon-precheck.md) |
| Diabetic retinopathy, GSE87433 | 31 October 2026 | All 17 arrays: six normal, six diabetic, five treated | All 24,973 supplied features | [DR report](docs/diabetic-retinopathy-precheck.md) |
| Experimental autoimmune uveitis, GSE284505 | 2 November 2026 | Six selected samples: three CON and three EAU; sleep-deprivation groups excluded by design | 18,187 of 56,859 | [Uveitis report](docs/autoimmune-uveitis-precheck.md) |
| AMD-related CFH-deficiency model, GSE38671 | 3 November 2026 | 16 neuroretina arrays: four per young/aged WT/Cfh-/- group | All 35,556 array features | [AMD report](docs/amd-precheck.md) |
| Light-induced retinal damage, GSE314043 | 3 November 2026 | All 10 libraries: five CTRL-vehicle and five LD-vehicle | 16,852 of 19,999 | [Light-damage report](docs/light-damage-precheck.md) |

The glaucoma decision retains R-08 with its distinct profile documented and omits the optional correlation check. Glaucoma stages and sex remain recorded; each condition/time pair consists of one male pool and one female pool. The RP decision retains all seven libraries and its later genotype model accounts for sex. The RP precheck report includes the recorded checkpoint audit and stored TMM-factor table. LHON retains all 12 libraries across WT, YARS2, COI and YARS2/COI, with its precheck metadata, scaling factors, session record and restore verification documented.

DR retains all 17 arrays for exploratory QC, with the six-versus-six untreated disease contrast and five treated arrays kept conceptually separate. Its native Entrez IDs were aligned with supplied platform gene symbols/names. The boxplot and density plots were reviewed as correct; the author verified all 17 density curves against their corresponding sample columns.

Autoimmune uveitis retains all six selected CON/EAU samples. Original labels mix gene names and Ensembl IDs. After BioMart connections failed, Ensembl REST returned 1,440 names, 24 records without names and 84 unresolved IDs from the 1,548-ID subset. Original labels and lookup gaps are preserved in a separate annotation table; the 55,311 other labels were not queried. Both sleep-deprivation groups remain outside this comparison.

AMD preserves array feature/transcript-cluster IDs and their multiple supplied gene assignments; its expanded annotation is separate from expression measurements. All 16 arrays remain included, with GSM947360 flagged for later inclusion review. Light damage preserves Ensembl IDs, uses a fixed release-116 GTF after REST failures, and retains all 10 libraries, including the control-like LD-vehicle-5. See the [dataset distinctions](docs/datasets.md#measurement-and-identifier-distinctions).

All seven datasets now have group-coloured boxplots. The [shared restyle script and record](docs/boxplot-styling.md) cover the six revised plots; light damage already used this style. Sample keys preserve original IDs and plot order.

These prechecks use deposited expression measurements and metadata for the seven terrestrial datasets. The separate flight precheck records 17 selected libraries. [OSD-758 precheck](docs/osdr758/README.md)

### AMD and light damage prechecks

The [AMD report](docs/amd-precheck.md) documents neuroretina-only selection, the log2 microarray branch, probe-group IDs, ambiguous/truncated annotations and the recorded sample flag. The [light-damage report](docs/light-damage-precheck.md) records count normalization, fixed GTF annotation, all 10 retained libraries and existing independent verification.

Light-damage projections show partial condition separation and variation among LD libraries. LD-vehicle-5 is near controls and remains retained. [Light-damage PCA](figures/light_damage/Light_Damage_GSE314043_TMM_PCA_top500.png)

</details>

## Explore the repository

<details>
<summary>Script index, metadata and dataset-specific guides</summary>

- [Earlier analysis scripts and dated backups](scripts/archive/).
- [Within-study results and 12-comparison summary](docs/withinstudy/README.md), [analysis scripts](scripts/withinstudy/), [result tables](tables/withinstudy/) and [figures](figures/withinstudy/).
- [Glaucoma script](scripts/01_Glaucoma_GSE241782_precheck.R) and [execution guide](docs/running-the-precheck.md).
- [RP script](scripts/02_Retinitis_Pigmentosa_GSE156533_precheck.R), [recovery helper](scripts/03_Retinitis_Pigmentosa_GSE156533_reload.R) and [execution guide](docs/running-the-rp-precheck.md).
- [Shared projection-styling script](scripts/04_Restyle_Retinal_MDS_PCA.R) and [style/reuse instructions](docs/projection-styling.md).
- [Glaucoma metadata](metadata/GSE241782_samples_all.csv) and [RP metadata](metadata/GSE156533_samples_matched.csv), both derived from public source records.
- [LHON teaching script](scripts/05_LHON_GSE304707_precheck.R), [reload helper](scripts/06_LHON_GSE304707_reload.R), [style helper](scripts/07_LHON_projection_style.R), [save/export helper](scripts/08_LHON_GSE304707_save_final.R) and [execution guide](docs/running-the-lhon-precheck.md).
- [LHON matched metadata](metadata/GSE304707_samples_matched.csv), [saved TMM factors](tables/lhon/LHON_GSE304707_TMM_summary_2026-10-30.csv) and [session record](docs/lhon-sessionInfo_2026-10-30.txt).
- [DR teaching script](scripts/09_Diabetic_Retinopathy_precheck.R), [resume helper](scripts/10_Diabetic_Retinopathy_GSE87433_resume.R), [style helper](scripts/11_DR_projection_style.R) and [execution guide](docs/running-the-dr-precheck.md).
- [Uveitis teaching script](scripts/12_Autoimmune_Uveitis_GSE284505_precheck.R), [reload helper](scripts/13_Autoimmune_Uveitis_GSE284505_reload.R), [style helper](scripts/14_Uveitis_projection_style.R), [save/export helper](scripts/15_Autoimmune_Uveitis_GSE284505_save.R) and [execution guide](docs/running-the-uveitis-precheck.md).
- [Uveitis selected metadata](metadata/GSE284505_samples_selected.csv), [gene annotation](tables/autoimmune_uveitis/Autoimmune_Uveitis_GSE284505_gene_annotation_2026-11-02.csv), [TMM summary](tables/autoimmune_uveitis/Autoimmune_Uveitis_GSE284505_TMM_summary_2026-11-02.csv) and [session record](docs/uveitis-sessionInfo_2026-11-02.txt).
- [AMD teaching script](scripts/16_AMD_GSE38671_precheck.R), [reload](scripts/17_AMD_GSE38671_reload.R), [annotation formatter](scripts/18_AMD_annotation_format.R), [projection helper](scripts/19_AMD_projection_style.R), [save helper](scripts/24_AMD_save_session.R) and [execution guide](docs/running-the-amd-precheck.md).
- [Light-damage script](scripts/20_Light_Damage_GSE314043_precheck.R), [reload](scripts/21_Light_Damage_GSE314043_reload.R), [verification](scripts/22_Light_Damage_GSE314043_verify.R), [saved result tables](tables/light_damage/) and [execution guide](docs/running-the-light-damage-precheck.md).
- [Group-coloured boxplots and sample keys](docs/boxplot-styling.md).
- [Dataset panel, measurement/identifier distinctions and input sources](docs/datasets.md).
- [GitHub first steps](docs/github-first-steps.md).

</details>

## Reproducibility and scope

The teaching scripts are designed for interactive execution one block at a time; separate complete workflows and restore helpers are documented where provided. Recorded observations and figures come from the completed project R runs. Repository synchronization reviewed those records, checked script syntax/portable paths and saved-table consistency, and preserved the original figure files without rerunning expression analysis. See [environment details](docs/environment.md).

Original expression inputs remain at their public repositories and on local storage; [data instructions](data/README.md) explain where to obtain them. Downloaded expression inputs, R workspace checkpoints and local files under `results/` are excluded from Git. Selected reviewed QC figures are versioned under `figures/`. The repository includes portable copies of the current scripts and an archive of earlier scripts.
