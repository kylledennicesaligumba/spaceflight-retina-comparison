# Running the glaucoma precheck

This script is intended for learning and review **one small block at a time**. Read each block and its output before continuing. The glaucoma processed-data precheck was completed on 28 October 2026; all eight libraries were retained. Recorded results in comments are reference observations from the author's run, rather than output produced by opening the file.

The file retains the interactive learning log. Earlier pending/error comments are superseded by the later result blocks and final closure decision. The optional correlation check was omitted at the agreed endpoint; it is not a pending requirement.

## Requirements

- R; the observed local environment is recorded in [environment.md](environment.md).
- `edgeR` for log-CPM, the QC data container, expression filtering and TMM composition scaling.
- `limma` for the exported `plotMDS()` generic used for MDS and its documented common-gene PCA mode.
- The two manually downloaded, extracted [GSE241782 files](datasets.md#glaucoma-inputs).

If packages are unavailable, use the installation instructions on the official Bioconductor pages for [edgeR](https://bioconductor.org/packages/edgeR/) and [limma](https://bioconductor.org/packages/limma/). Package installation is a separate setup step; the precheck does not install packages or download datasets automatically.

## Set the working directory

Open R or RStudio and set the working directory to your local `spaceflight-retina-comparison` repository. In R, a working directory is the starting folder for relative paths such as `data/raw` and `results/glaucoma`.

For example, replace the placeholder below with your own local repository location:

```r
setwd("/path/to/spaceflight-retina-comparison")
getwd()
```

Open `scripts/01_Glaucoma_GSE241782_precheck.R`. Execute its blocks interactively, starting at Step 1, rather than sourcing the entire file at once.

## Locate the data

Choose either of these approaches:

1. Put the extracted input files in `data/raw/Glaucoma_GSE241782/` under the repository; or
2. Keep the data in an existing folder and set the environment variable in R before Step 1:

```r
Sys.setenv(GLAUCOMA_DATA_DIR = "/path/to/Glaucoma_GSE241782")
```

The folder must contain these exact filenames:

```text
GSE241782_gene_count_matrix.csv
GSE241782_family.soft
```

Step 1 checks that both files exist. A successful file check establishes location, not biological quality. The repository does not include the downloaded expression matrix.

## Checkpoints and generated files

Step 6 creates `results/glaucoma/` if needed and saves metadata plus an R workspace checkpoint there. This checkpoint records the state reached at Step 6; it does not automatically refresh after later steps. A pause block saves a separately dated 27 October checkpoint, and the final block saves the completed 28 October checkpoint. Run each save block only when you intend to preserve that stage. Rerunning a save block replaces its same-named file.

`save.image()` saves the entire active workspace. Use a clean, project-specific R session to avoid mixing unrelated objects. Workspace files are local conveniences; they are excluded from Git. The script and documented input sources are the reproducibility record.

The repository version uses relative output paths. The original working script and checkpoint on the external drive remain separate. To restore a repository-generated checkpoint in a later session, start in the repository root and use the commented `load()` example in Step 6.

## Figure exports and completed report

- Steps 9C–9D display and save the initial log-CPM boxplot.
- Steps 9E–9F display and save log-CPM density curves.
- Step 11C displays the pairwise MDS; Step 11D includes a complete calculation-and-save block.
- Step 11E calculates and saves the complementary common-gene PCA.

All PNG export blocks write to local `results/glaucoma/`. The MDS and PCA blocks use the same filtered/TMM QC object; they do not refilter genes or reestimate normalization factors. Opening a PNG does not refresh a workspace checkpoint.

The reviewed original figures are committed under `figures/glaucoma/` and embedded in the [completed report](glaucoma-precheck.md). Regenerating local files under `results/` does not automatically replace those curated copies. No new figure generation is needed merely to review the completed precheck.
