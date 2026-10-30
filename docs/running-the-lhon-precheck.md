# Running or restoring the LHON precheck

The GSE304707 processed-data precheck was completed on 30 October 2026 with all 12 libraries retained. The [completed report](lhon-precheck.md) records the 17,860-gene QC copy, final figures, exported tables and limitations.

## Fresh block-by-block reconstruction

Start R in the repository root. Required packages are `xml2`, `edgeR`, `limma`, `ggplot2` and `ggrepel`; the [saved session information](lhon-sessionInfo_2026-10-30.txt) records the actual completed-session versions. No package installation or data download happens automatically.

The input directory defaults to `data/raw/LHON_GSE304707/`. To use existing storage, set this before beginning:

```r
# Step 0: Point to the downloaded LHON inputs without copying them.
Sys.setenv(LHON_DATA_DIR = "/path/to/LHON_GSE304707")
```

The main script expects `GSE304707_gene_count_matrix.csv` and the locally named archive `GSE304707_family`, containing `GSE304707_family.xml`. Preserve original filenames; adjust the archive filename assignment if your download uses another name.

Open [05_LHON_GSE304707_precheck.R](../scripts/05_LHON_GSE304707_precheck.R) and execute one explained block at a time. It creates local `results/lhon/`, extracts the XML into its `metadata/` subfolder, and preserves the original counts while constructing separate QC objects. It does not choose a final disease comparison or run differential expression.

Step 12 sources [07_LHON_projection_style.R](../scripts/07_LHON_projection_style.R) from the repository's `scripts/` directory. This helper defines formatting using existing scores; the main script calculates the MDS/PCA, displays both and explicitly saves the two PNGs. Initial boxplot/density exports are included in Step 7. Reviewed copies under `figures/lhon/` are separate from regenerated local outputs.

## Save the completed stage

After all completed-stage objects and plots exist in the active R session, use:

```r
# Final save: export existing LHON objects, matched metadata, factors and session info.
# This block checks the completed state and does not repeat expression calculations.
source(file.path("scripts", "08_LHON_GSE304707_save_final.R"))
```

The helper requires the full and filtered data, aligned annotations, normalized log-CPM and projections/plots. It saves an explicit object list excluding XML document/node pointers, then exports CSV tables and `sessionInfo()`. Four file-existence checks confirm the writes. Rerunning replaces the same dated local files and preserves the earlier checkpoint.

## Restore the final checkpoint

Use [06_LHON_GSE304707_reload.R](../scripts/06_LHON_GSE304707_reload.R) rather than sourcing the whole teaching script. The portable copy targets local `results/lhon/LHON_GSE304707_precheck_checkpoint_2026-10-30.RData`; make your completed checkpoint available there before running it. GitHub does not contain the binary checkpoint.

```r
# Restore: load saved results, reopen their XML metadata and check dimensions/modes.
# Use a project-specific R session; this is a restore operation, not a QC rerun.
source(file.path("scripts", "06_LHON_GSE304707_reload.R"))
```

Saved path objects retain the values used when the checkpoint was created. An original MaxineDrive checkpoint therefore expects its existing XML path on that drive. Ensure that file is available, and deliberately reset paths before subsequent writes if moving the project. The reload helper reopens XML pointers, checks original/filtered dimensions and projection modes, and prints alignment checks; it does not refit the projections.

The original 29 October checkpoint is partial and is not the completed restart source. No new precheck steps are required simply to review or restore this closed work.
