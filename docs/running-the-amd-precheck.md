# Running or restoring the AMD microarray precheck

The [completed report](amd-precheck.md) distinguishes the 16 selected neuroretina arrays, one log2 QC transformation, array feature IDs and multipart annotations. This workflow uses no CPM, TMM or RNA-seq feature filter.

Start R in the repository root. Required packages are `xml2`, `limma`, `ggplot2` and `ggrepel`. The default input folder is `data/raw/AMD_GSE38671/`, containing the locally named `GSE38671_series_matrix` and `GSE38671_family`. The family archive supplies `GSE38671_family.xml` and `GPL6246-tbl-1.txt`. Adjust filename assignments if your downloads have different extensions.

```r
# Use existing AMD downloads instead of copying the input files.
Sys.setenv(AMD_DATA_DIR = "/path/to/AMD_GSE38671")
```

Open [16_AMD_GSE38671_precheck.R](../scripts/16_AMD_GSE38671_precheck.R) and run one explained block at a time. It captures `precheck_project_dir` before changing R's working directory to `results/amd/`. It loads the renamed annotation/projection helpers from the captured repository's `scripts/` folder. Historical pending/error comments describe the learning sequence; the completed report records final observations.

The current formatter already handles the known source truncation marker. Earlier comments describing the original formatting failure remain as history; they do not mean the current formatter must fail at that point. Unexplained malformed records still stop for review. Full formatting creates a large table because original strings are repeated for each assignment; allow space for the existing approximately 1.82 GB output. Assignment rows must remain separate from expression rows.

The main script generates the original initial boxplot. Run the [shared restyler](boxplot-styling.md) separately to reproduce the versioned group-coloured replacement. Its default is the repository's `results/` folder, or the output tree selected by `PREQUALITYCHECKS_DIR`.

The main script's last block loads [24_AMD_save_session.R](../scripts/24_AMD_save_session.R). It creates timestamped teaching-script/checkpoint files, saves package versions, attempts a separate console-history save, and updates `AMD_GSE38671_latest_checkpoint.txt` after checkpoint creation. XML readers are excluded. The portable helper uses the captured repository folder for its teaching-script source and `results/amd/` output. When running the save helper independently, start at the repository root; after the main script, use its existing `precheck_project_dir`. No new timestamped checkpoint was confirmed during this upload; the existing 2 November checkpoint predates the later plots.

```r
# Restore the existing checkpoint in a project-specific R session.
# Check the saved stage against the completion report before continuing.
source(file.path("scripts", "17_AMD_GSE38671_reload.R"))
```

The portable reload helper reads the latest-checkpoint pointer under `results/amd/` when available and falls back to the 2 November checkpoint otherwise. It restores objects and reopens XML readers. A copied pointer contains the path written when saved; update it deliberately if the checkpoint has moved. Saved path objects retain their original values: a MaxineDrive checkpoint requires the existing original output/XML locations unless you deliberately relocate and reset them. No completed-state restore test was run during this upload.

Downloaded expression, original/expanded annotation and binary checkpoints remain excluded from Git. The four reviewed PNGs and selected metadata are separate curated copies; regenerating local outputs does not automatically change GitHub.
