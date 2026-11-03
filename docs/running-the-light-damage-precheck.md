# Running or restoring the light-damage precheck

The [completed report](light-damage-precheck.md) records all 10 retained libraries, the 16,852-row QC copy, fixed GTF annotation and saved verification. This completed workflow was run in a separate R process; uploading it did not rerun the analysis.

Start from the repository root, using a separate project R session/process. Required packages are `xml2`, `edgeR`, `limma`, `ggplot2` and `ggrepel`; see the [actual session record](light-damage-sessionInfo_2026-11-03.txt). The default input folder is `data/raw/Light_Damage_GSE314043/`, containing `GSE314043_genecount.txt` and the locally named `GSE314043_family` archive. No extraction of the already-uncompressed count file is needed.

```r
# Use existing downloaded counts and the GEO family archive.
Sys.setenv(LIGHT_DAMAGE_DATA_DIR = "/path/to/Light_Damage_GSE314043")
```

[20_Light_Damage_GSE314043_precheck.R](../scripts/20_Light_Damage_GSE314043_precheck.R) preserves the source workflow with portable paths and writes under `results/light_damage/` without changing R's working directory. It imports the deposited counts, matches all 10 samples, creates distributions, adds annotation, filters expression, estimates TMM and saves projections/tables/checkpoint.

For a fresh full reconstruction, run this deliberately from the repository root:

```r
# Reconstruct the complete light-damage precheck and replace its named local outputs.
source(file.path("scripts", "20_Light_Damage_GSE314043_precheck.R"))
```

The script downloads the official release-116 GRCm39 GTF if it is absent from the local annotation output folder. That download, reading the GTF, and all expression calculations occur only when the reconstruction is executed. The successful final route does not depend on Ensembl REST availability. Original IDs, missing annotations and all raw rows are preserved; the expression filter creates a separate QC copy.

To reopen an existing completed checkpoint instead, make it available under `results/light_damage/` and use:

```r
# Restore saved results and confirm sample/gene order without refitting QC.
source(file.path("scripts", "21_Light_Damage_GSE314043_reload.R"))
```

The helper expects the dated 3 November checkpoint. Saved path objects retain their original values; a MaxineDrive checkpoint may point subsequent file access/writes back to the original output tree. Keep that tree available, or deliberately relocate files and reset paths before using a moved checkpoint. Saved redraw functions are available but are not called automatically.

[22_Light_Damage_GSE314043_verify.R](../scripts/22_Light_Damage_GSE314043_verify.R) sources the portable reload helper, checks the original files against the locally saved manifest, compares count/QC exports, independently calculates PCA and rewrites its verification record/top-500 list. It needs the full local result tree, including raw-count CSV and manifest, which are not uploaded. The existing verification record is included for provenance; no verification rerun was performed during this upload.

The four PNGs under `figures/light_damage/` and selected readable tables under `tables/light_damage/` are curated copies. The full GTF, original/downloaded expression measurements, normalized-expression RDS and workspace remain excluded from Git.
