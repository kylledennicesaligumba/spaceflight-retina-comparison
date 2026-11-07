# Retinal flight and terrestrial disease comparison

The completed analysis compares 12 contrasts from seven terrestrial studies with the OSD-758 mouse retinal flight 0 g versus habitat ground control response. The shared analysis uses 12,127 harmonized mouse Entrez genes and 799 matched Reactome pathways. [Final results and limits](records/Main_analysis_results_and_limits_2026-11-07.md)

Two-week experimental glaucoma had the highest observed gene-score correlation (Spearman ρ = 0.243) and pathway-profile correlation (ρ = 0.302). Its gene correlation ranked first with the larger pairwise gene set as well (ρ = 0.233). It ranked first in 545 of 552 pathway neighborhood deletions and second in the other seven. [Complete comparison table](../../tables/main_analysis/Main_analysis_final_comparison_summary.csv)

The deletion analysis measures sensitivity to overlapping pathway groups. Biological sampling uncertainty and the difference between the leading models remain unquantified. The recorded shared signals include matrix and collagen programs and lower neuronal and synaptic enrichment. [Coverage sensitivity](../../figures/main_analysis/Main_analysis_OSD758_gene_coverage_sensitivity.png) · [Pathway rank stability](../../figures/main_analysis/Main_analysis_OSD758_pathway_rank_stability.png)

| Stage | Materials |
| --- | --- |
| Flight reference | [OSD-758 precheck](../osdr758/README.md) and [within-study results](../withinstudy/osdr758.md) |
| Terrestrial models | [Within-study results for all 12 contrasts](../withinstudy/README.md) |
| Identifier mapping and coverage | [Gene coverage](../../tables/main_analysis/Main_analysis_gene_coverage.csv) and [mapping and alignment tables](../../tables/main_analysis/) |
| Gene-level comparison | [Spearman results](../../tables/main_analysis/Main_analysis_OSD758_gene_Spearman_common.csv) and [gene similarity figure](../../figures/main_analysis/Main_analysis_OSD758_gene_Spearman_common.png) |
| Reactome comparison | [Enrichment results](../../tables/main_analysis/Main_analysis_Reactome_fgsea_all_contrasts.csv) and [pathway-profile correlations](../../tables/main_analysis/Main_analysis_OSD758_pathway_Spearman_common.csv) |
| Shared and opposing programs | [Positive programs](../../figures/main_analysis/Main_analysis_OSD758_positive_programs_full_panel.png), [negative programs](../../figures/main_analysis/Main_analysis_OSD758_negative_programs_full_panel.png) and [leading-edge summary](../../tables/main_analysis/Main_analysis_OSD758_leading_edge_summary_full_panel.csv) |
| Sensitivity | [Pairwise gene coverage](../../tables/main_analysis/Main_analysis_OSD758_gene_coverage_sensitivity.csv) and [pathway rank ranges](../../tables/main_analysis/Main_analysis_pathway_rank_stability_summary.csv) |

## Running the recorded workflow

Start from the repository root. The main entry point registers the 13 ranked tables under `tables/withinstudy/` and calls modules 02 through 26 in order. It writes to `results/main_analysis/` by default. Use `MAIN_ANALYSIS_DIR` to select an absolute output directory. Modules 24 and 25 implement the final sensitivity analyses, and module 26 consolidates their saved results. [All 27 main-analysis scripts](../../scripts/main_analysis/)

```r
# Run the full recorded comparison using the ranked tables in the repository.
source(file.path("scripts", "main_analysis", "01_main_analysis.R"))
```

The source scripts contain package setup and use the versions recorded with the analysis. The recorded Reactome run used R 4.5.3, Bioconductor 3.22, AnnotationDbi 1.72.0, fgsea 1.36.2 and reactome.db 1.95.0. [Version records](../../tables/main_analysis/Main_analysis_Reactome_versions.csv) · [Session and method records](records/)

For a saved local run, set `MAIN_ANALYSIS_DIR` to the absolute directory containing the checkpoints and use `00_resume_main_analysis.R`. The resume script reads the latest available checkpoint. Individual sensitivity modules also check stored file hashes, so checkpoints from the original run expect their recorded input locations. Raw inputs and R workspace files are kept in local storage. [Resume script](../../scripts/main_analysis/00_resume_main_analysis.R)
