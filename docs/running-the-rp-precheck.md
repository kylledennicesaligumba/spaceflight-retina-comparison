# Running the retinitis pigmentosa precheck

The GSE156533 processed-data precheck was completed on 29 October 2026. All seven libraries were retained; the filtered QC copy contains 17,351 gene rows. The [completion report](retinitis-pigmentosa-precheck.md) records the observations and limitations.

The main script retains the interactive learning log. Work through its blocks in order and review the output. Earlier pending notes are superseded by later results and the closure decision. Merely opening or parsing a script does not run it.

## Inputs and packages

Start R in the repository root. Use R with `readxl`, `xml2`, `edgeR`, `limma`, `ggplot2` and `ggrepel`; observed versions are recorded in [environment.md](environment.md). Use the packages' official installation instructions if they are unavailable. The scripts do not install packages or download data automatically.

The [input-source guide](datasets.md#retinitis-pigmentosa-inputs) identifies the count workbook and GEO XML archive. By default the portable scripts look in `data/raw/Retinitis_Pigmentosa_GSE156533/`. To keep the files elsewhere, set this in R before starting:

```r
Sys.setenv(RP_DATA_DIR = "/path/to/Retinitis_Pigmentosa_GSE156533")
```

The actual downloaded filenames used in this workflow are:

```text
GSE156533_readcounts_organized_by_largest_to_smallest_count.xlsx
GSE156533_family
```

The second file is a tar archive containing `GSE156533_family.xml`, despite its extensionless local name. Preserve the downloaded file; if your archive has a different name, adjust the `rp_metadata_archive` filename assignment. Step 6A extracts only the XML member into local `results/retinitis_pigmentosa/metadata/`.

## Main workflow

Open [02_Retinitis_Pigmentosa_GSE156533_precheck.R](../scripts/02_Retinitis_Pigmentosa_GSE156533_precheck.R). Execute the blocks interactively. The script creates `results/retinitis_pigmentosa/` for local outputs.

- Steps 0–8 inspect files, import counts, extract XML metadata and align the seven sample records.
- Steps 9A–9D cover library totals, detected genes and initial boxplot/density exports.
- Steps 10–12 check supplied gene annotations, construct the QC data object and filter low-expression rows.
- Steps 13–15 cover TMM-adjusted log-CPM, pairwise MDS and common-top-500 PCA.
- The shared styling step updates RP MDS and glaucoma MDS/PCA using existing scores. See [projection-styling.md](projection-styling.md) for its session/checkpoint prerequisites.

The main RP script draws the final PCA appearance directly. Its initial MDS block uses the earlier appearance; the shared styling script supplies the revised MDS. Both approaches preserve the stored projection coordinates.

## Recovery and checkpoints

[03_Retinitis_Pigmentosa_GSE156533_reload.R](../scripts/03_Retinitis_Pigmentosa_GSE156533_reload.R) reconstructs the initial-boxplot stage by rereading the original workbook and already extracted XML. It is not a loader for the completed checkpoint. Do not source it at the end of a completed session merely to obtain the latest state.

The main script's final save block initializes its base checkpoint object list if needed, then extends it with completed QC and plotting objects. The base list is the same list used in the recovery helper. Complete the shared styling step first so `rp_mds_plot` and its drawing helper exist. This checkpoint bookkeeping was made explicit in the portable copy; expression calculations were preserved.

The final checkpoint filename is `RP_GSE156533_precheck_checkpoint_2026-10-29.RData`. The original local file was confirmed present during synchronization; its contents were not loaded. The saved object selection excludes XML external pointers and the isolated glaucoma environment. Recreate XML documents from the source file if XML inspection is needed after a restart.

For a checkpoint generated under this repository's relative paths, restore it in a project-specific R session from the repository root:

```r
load(file.path("results", "retinitis_pigmentosa",
               "RP_GSE156533_precheck_checkpoint_2026-10-29.RData"))
```

An original checkpoint from another location may contain its original file paths; reset path objects deliberately before writing further outputs. Do not load an unrelated checkpoint into a working analysis session.

## Versioned artifacts

Local outputs and checkpoints under `results/` are ignored by Git. Reviewed figure copies are versioned in `figures/retinitis_pigmentosa/`; the matched sample table is under `metadata/`. Regenerating local results does not automatically replace those curated copies. No new calculations are needed simply to review the completed report.
