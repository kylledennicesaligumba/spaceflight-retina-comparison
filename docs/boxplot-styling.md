# Group-coloured initial boxplots

Updated 3 November 2026. All seven terrestrial datasets now have boxplots coloured by their recorded groups. [23_Restyle_Study_Boxplots.R](../scripts/23_Restyle_Study_Boxplots.R) is the portable copy of the shared `Restyle_Study_Boxplots.R` from the larger prequalitychecks folder. It updates six existing studies; light damage already had the chosen style and is excluded from this script.

| Dataset | Display groups | Expression used |
| --- | --- | --- |
| Glaucoma | Naive baseline, glaucoma 3D, 2W, 6W; sex and library ID in labels | Original all-row log-CPM, before filtering/TMM |
| RP | WT, P23H; sample ID and sex in labels | Original all-row log-CPM, before filtering/TMM |
| LHON-linked model | WT, YARS2, COI, YARS2/COI | Original all-row log-CPM, before filtering/TMM |
| DR | Normal, diabetic, diabetic + transplantation | Supplied processed microarray values; no added log2/CPM/TMM |
| Autoimmune uveitis | CON, EAU | Original all-row log-CPM, before filtering/TMM |
| AMD-related model | Young WT, young Cfh-/-, aged WT, aged Cfh-/- | Single log2 QC copy of neuroretina microarray intensities |
| Light damage | CTRL-vehicle, LD-vehicle | Original all-row log-CPM, before filtering/TMM; already styled |

The palette uses blue `#0072B2`, orange `#D55E00`, green `#009E73` and purple `#CC79A7` in each study's documented group order. Colours label groups within a study; a shared colour does not mean two studies measure the same condition. Control groups appear first, retaining original order within each group. Individual extreme-value dots are hidden, as before, without deleting feature rows.

The [six saved sample keys](../metadata/boxplot_sample_keys/) retain displayed positions, original matrix-column positions, original IDs, labels, groups, colours and recorded metadata. They provide the complete text when long labels are difficult to read in a PNG. The [restyle record](../tables/boxplots/Boxplot_restyle_record_2026-11-03.csv) records `box_statistics_unchanged = TRUE` for each of the six studies. Only its machine-specific paths were adapted to repository-relative `results/` paths for version control.

The source restyle script loads each checkpoint into a private environment and checks sample alignment, dimensions and unchanged boxplot statistics after reordering. AMD's older checkpoint predates its plot, so that branch reconstructs only the reviewed neuroretina subset and applies log2 once to the original supplied positive intensities. It does not apply log2 to a previously logged matrix. No filtering, TMM, differential expression or projection fitting is added by this restyler.

## Reproduce deliberately

Start R in the repository root and make the six original dated checkpoints available in their corresponding folders under `results/`. Alternatively, point to an existing complete output tree:

```r
# Use an existing prequalitychecks tree containing the six required checkpoints.
Sys.setenv(PREQUALITYCHECKS_DIR = "/path/to/prequalitychecks")
```

```r
# Replace the six named boxplots and write sample keys/export records.
# Existing PNGs are backed up once in each study's dated archive folder.
source(file.path("scripts", "23_Restyle_Study_Boxplots.R"))
```

This is an output-writing operation. The script writes a temporary PNG before replacing each usual filename, closes each device, and preserves the first existing PNG in a dated archive. It does not change the original analysis scripts or the user's active R objects. Its default output root is `results/`; setting the environment variable directs writes to that selected tree.

The initial teaching scripts retain their original plotting blocks. Reconstructing an earlier workflow produces its earlier boxplot appearance until this restyler is run. Versioned figures under `figures/` are separate curated copies, and are not replaced automatically by rerunning a local script.

Repository synchronization copied the existing revised images byte-for-byte and checked saved keys/records. It did not load the checkpoints, rerun the restyle or generate new scientific plots. Earlier figures remain available in Git history and in the original local archives.
