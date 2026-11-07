# OSD-758 flight reference precheck

The reviewed input had 56,840 rows and 35 libraries. The comparison selected five flight 0 g libraries and 12 habitat ground controls. The precheck separated 92 ERCC controls, then retained 21,335 of 56,748 mouse Ensembl rows through expression filtering. All 17 selected libraries carried forward, including HGC05 with its recorded projection flag. [OSD-758 precheck scripts](../../scripts/prechecks/osdr758/)

The plots show all selected libraries. MDS axes explain 20% and 15% and the top-500 PCA axes explain 25% and 19%, as displayed in the saved figures. [Boxplot](../../figures/osdr758/OSD758_logCPM_boxplot.png) · [Density](../../figures/osdr758/OSD758_logCPM_density.png) · [MDS](../../figures/osdr758/OSD758_TMM_MDS.png) · [PCA](../../figures/osdr758/OSD758_TMM_PCA_top500.png)

Start R in the repository root. The precheck reads `GLDS-664_rna-seq_RSEM_Unnormalized_Counts.csv`, `GLDS-664_rna-seq_SampleTable.csv` and `GLDS-664_rna-seq_MHU-8_RTN_metadata.csv` from `data/raw/osdr758/` or the absolute directory set by `OSD758_DATA_DIR`. Follow the recorded blocks in order. Outputs and the saved precheck checkpoint go under `results/osdr758/`. [Precheck entry point](../../scripts/prechecks/osdr758/01_OSD758_precheck.R) · [NASA dataset record](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758)

The saved model uses the precheck's TMM factors and retained rows. [Flight within-study analysis](../withinstudy/osdr758.md)
