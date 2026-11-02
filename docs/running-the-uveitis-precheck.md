# Running or restoring the autoimmune uveitis precheck

The GSE284505 precheck was completed on 2 November 2026. It retains three CON and three EAU samples, excludes the two sleep-deprivation groups from the selected comparison, and preserves the supplied mixed-format gene labels. See the [completed report](autoimmune-uveitis-precheck.md).

## Inputs and packages

Start R in the repository root. The main/reload scripts capture that folder as `precheck_project_dir` before changing the working directory. New outputs go to `results/autoimmune_uveitis/`; helper scripts are loaded from the captured repository's `scripts/` folder.

The default input folder is `data/raw/Autoimmune_Uveitis_GSE284505/`. To use existing downloads, set an absolute folder before starting:

```r
# Use the existing count and metadata downloads without copying them.
Sys.setenv(UVEITIS_DATA_DIR = "/path/to/Autoimmune_Uveitis_GSE284505")
```

The script expects `GSE284505_R.HTSeq.counts.txt` and the locally named `GSE284505_family` archive containing `GSE284505_family.xml`. Preserve filenames, adjusting the archive assignment if your download uses a different name. The script extracts the XML under the local output folder.

Packages used are `xml2`, `edgeR`, `limma`, `ggplot2`, `ggrepel`, `httr` and `jsonlite`. See the [saved session information](uveitis-sessionInfo_2026-11-02.txt) for the actual completed-run versions. `biomaRt` appears in that session after failed connection attempts, but the successful annotation code uses `httr`/`jsonlite` with Ensembl REST. No packages are installed automatically.

## Interactive reconstruction

Open [12_Autoimmune_Uveitis_GSE284505_precheck.R](../scripts/12_Autoimmune_Uveitis_GSE284505_precheck.R) and run one explained block at a time, reviewing outputs before continuing. Historical pending comments refer to the stage at which they were written; use the completed report for the final state. Sourcing the full script repeats import, live annotation queries, calculations and output writes; use the reload helper to reopen existing results.

The main script checks count integrity, matches sample metadata, selects CON/EAU columns, reviews initial all-row log-CPM distributions, preserves mixed labels, retrieves the Ensembl subset's annotation, filters expression, estimates TMM and calculates MDS/PCA. Missing annotation names are kept; they do not define the expression filter.

Ensembl lookup results can change with database releases and service availability. The [saved 2 November annotation table](../tables/autoimmune_uveitis/Autoimmune_Uveitis_GSE284505_gene_annotation_2026-11-02.csv) preserves release 116/GRCm39 results, including all unresolved original IDs. Source processing notes report GRCm38. The lookup is added annotation, not a change to the count reference or a finished cross-study mapping.

The main script loads [14_Uveitis_projection_style.R](../scripts/14_Uveitis_projection_style.R), calculates both projections and explicitly displays/saves their PNGs. The helper defines formatting only. Initial distribution exports precede filtering/TMM; projection exports use the 18,187-row TMM log-CPM matrix.

## Save the completed state

The final main-script block loads [15_Autoimmune_Uveitis_GSE284505_save.R](../scripts/15_Autoimmune_Uveitis_GSE284505_save.R). To save an already completed active session deliberately:

```r
# Export the existing completed objects, metadata, annotation, factors and session record.
source(file.path(precheck_project_dir, "scripts", "15_Autoimmune_Uveitis_GSE284505_save.R"))
```

The helper checks required objects, dimensions and alignments, then saves the dated checkpoint and readable exports. It excludes live XML reader objects and the HTTP response; parsed metadata and annotation results are retained. Rerunning replaces the same 2 November files. It does not recalculate QC or repeat API lookups.

## Restore without repeating analysis

Begin in the repository root and make the completed checkpoint available at `results/autoimmune_uveitis/Autoimmune_Uveitis_GSE284505_precheck_checkpoint_2026-11-02.RData`. Checkpoints are excluded from GitHub.

```r
# Restore saved results, reopen XML readers and check dimensions and alignment.
source(file.path("scripts", "13_Autoimmune_Uveitis_GSE284505_reload.R"))
```

Saved path objects keep their original values. A checkpoint made on MaxineDrive still expects that output folder and XML file there; ensure those locations are available or deliberately update the saved paths before using relocated files or writing new results. The helper's captured repository path resolves the renamed style script independently of saved output paths. It does not refit projections or repeat annotation requests. This upload confirmed checkpoint presence without running a restore test.

The four reviewed PNGs under `figures/autoimmune_uveitis/` are versioned separately from new outputs in `results/`. Updating a local figure or table does not automatically replace its curated GitHub copy. No additional precheck step is required by this completed upload.
