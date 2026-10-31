# Running or restoring the diabetic retinopathy precheck

The GSE87433 precheck uses **deposited processed microarray expression values**. **No additional log2 transformation, CPM calculation or TMM normalization is performed.** The code also does not add RNA-seq expression filtering or repeat source preprocessing. See the [completed report](diabetic-retinopathy-precheck.md).

## Inputs and packages

Start R in the repository root. The portable scripts capture that folder as `dr_project_dir` before the original workflow changes the working directory to its output folder. This keeps input paths and helper-script paths valid after `setwd()`.

The default input directory is `data/raw/Diabetic_Retinopathy_GSE87433/`. To use existing storage, set an absolute input path first:

```r
# Setup: use the existing DR downloads without copying or altering them.
Sys.setenv(DR_DATA_DIR = "/path/to/Diabetic_Retinopathy_GSE87433")
```

The actual local filenames are `GSE87433_series_matrix` and `GSE87433_family`. The family archive contains `GSE87433_family.xml` and the supplied `GPL20710-tbl-1.txt` annotation table. Preserve downloaded names, adjusting filename assignments if your files use other names. Extracted metadata/annotation and new outputs go under `results/diabetic_retinopathy/`.

The workflow uses base R for import, integrity checks and distributions, `xml2` for sample/platform metadata, `limma` for projections, and `ggplot2`/`ggrepel` for the approved appearance. Count-processing packages are not required for DR. See [environment.md](environment.md).

## Interactive reconstruction

Open [09_Diabetic_Retinopathy_precheck.R](../scripts/09_Diabetic_Retinopathy_precheck.R) and review one explained block and its output at a time. Historical pending notes in the learning script refer to their original stage; use the completed report for the final state. Do not source the full history merely to restore a completed session.

The sequence covers expression import and integrity, sample metadata and processing provenance, exact sample matching, distributions and density verification, supplied feature annotation and alignment, and then MDS/PCA with all 17 arrays. The annotation table has no header row; the script obtains its 15 field names from the reviewed XML definitions before matching IDs.

The main script sources [11_DR_projection_style.R](../scripts/11_DR_projection_style.R) through the captured project folder. This helper defines formatting/export functions only. The main code calculates both projections from `dr_expr_all`, displays them, and explicitly saves the PNGs. No extra transformation, TMM or feature filtering occurs between the supplied matrix and those exploratory projections.

## Checkpoints and restoration

The initial save block targets the 30 October checkpoint; the final block targets `DR_GSE87433_precheck_checkpoint_2026-10-31.RData`. Run a save only when you intend to preserve that state. Rerunning a block replaces its same-named file. The new file's presence was confirmed during synchronization; its saved contents were not loaded for this upload.

To restore a repository-generated final checkpoint, first start in the repository root and make the local checkpoint available under `results/diabetic_retinopathy/`. Then use:

```r
# Restore saved results and rebuild XML/plotting helpers; do not rerun QC.
source(file.path("scripts", "10_Diabetic_Retinopathy_GSE87433_resume.R"))
```

The helper captures the project directory, loads the final checkpoint, restores its working directory, recreates XML pointers and loads the projection helper. It prints expression, metadata and annotation dimensions. Saved paths retain their original values: a checkpoint from MaxineDrive still expects its original output/XML/annotation locations unless you deliberately reset them. This repository's renamed script files correspond to the local main, resume and style scripts; their calculations are preserved.

## Versioned figures

Local results and RData files remain excluded from Git. The five reviewed PNGs are copied to `figures/diabetic_retinopathy/` and embedded in the report. Regenerating a local image does not automatically replace the curated GitHub copy. The large platform table remains an input reference extracted from the source archive, rather than a new portfolio result.
