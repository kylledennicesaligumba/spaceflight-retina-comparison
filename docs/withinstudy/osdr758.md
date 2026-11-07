# OSD-758 flight versus habitat ground control

The reference models five flight 0 g and 12 habitat ground control retinal libraries with a group-only robust edgeR quasi-likelihood model. The comparison uses the 21,335 rows and TMM factors from the reviewed precheck. The study collected retinal tissue after flight return. [Precheck record](../osdr758/README.md) · [Study methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11593819/)

At BH FDR below 0.05, 373 genes met the threshold, with 286 higher and 87 lower in flight. The signed ranking retains all 21,335 tested rows and supplies the fixed reference for the cross-study analysis. [Complete results](../../tables/withinstudy/osdr758/OSD758_FLT_0G_vs_HGC_all_genes.csv) · [DEG table](../../tables/withinstudy/osdr758/OSD758_FLT_0G_vs_HGC_DEGs_FDR005.csv) · [Full ranking](../../tables/withinstudy/osdr758/OSD758_FLT_0G_vs_HGC_ranked_all_genes.csv)

The figures include the saved dispersion diagnostic, genome-wide MD plot and selected-gene displays. [Dispersion](../../figures/withinstudy/osdr758/OSD758_QL_dispersion.png) · [MD](../../figures/withinstudy/osdr758/OSD758_MD_plot.png) · [Heatmap](../../figures/withinstudy/osdr758/OSD758_DEG_expression_heatmap_top30.png) · [Dot plots](../../figures/withinstudy/osdr758/OSD758_DEG_expression_dotplots_top10.png)

Start R in the repository root and follow the recorded blocks. The entry point reads the saved precheck under `results/osdr758/` or the absolute checkpoint path set by `OSD758_PRECHECK_FILE`. Outputs go under `results/withinstudy/osdr758/`. The nine scripts include separate reload helpers for the precheck input, the fitted model and the final result. [Scripts](../../scripts/withinstudy/osdr758/) · [Recorded session](records/osdr758/OSD758_sessionInfo_final_2026-11-06.txt)
