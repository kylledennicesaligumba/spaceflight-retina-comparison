# Running the glaucoma precheck

This script is intended for learning and review **one small block at a time**. Read each block and its output before continuing. Recorded results in comments are reference observations from the author's run, rather than output produced by opening the file.

## Requirements

- R; the observed local environment is recorded in [environment.md](environment.md).
- The `edgeR` package for the initial log-CPM plot in Step 9C. Earlier blocks use base R.
- The two manually downloaded, extracted [GSE241782 files](datasets.md#glaucoma-inputs).

If `edgeR` is unavailable, use the installation instructions on its [official Bioconductor page](https://bioconductor.org/packages/edgeR/). Package installation is a separate setup step; the precheck does not install packages or download datasets automatically.

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

Step 6 creates `results/glaucoma/` if needed and saves metadata plus an R workspace checkpoint there. This checkpoint records the state reached at Step 6; it does not automatically refresh after later steps. Rerunning the save block replaces the same output files.

`save.image()` saves the entire active workspace. Use a clean, project-specific R session to avoid mixing unrelated objects. Workspace files are local conveniences; they are excluded from Git. The script and documented input sources are the reproducibility record.

The repository version uses relative output paths. The original working script and checkpoint on the external drive remain separate. To restore a repository-generated checkpoint in a later session, start in the repository root and use the commented `load()` example in Step 6.

Step 9C displays a plot in R's graphics window. It does not save a plot file. Review that plot before extending the precheck.
