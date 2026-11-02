# Shared MDS/PCA appearance

The approved appearance uses actual sample points, readable F/M labels, connector lines and a right-side legend outside the plotting panel. `ggplot2` draws the existing scores and `ggrepel` separates text labels. See the [official ggrepel examples](https://ggrepel.slowkow.com/articles/examples.html) for the label controls.

The final revised exports were made on 29 October 2026. The repository stores the original PNGs from that session, without regenerating them during synchronization.

## What each script draws

The RP main script draws its PCA in the approved style. [04_Restyle_Retinal_MDS_PCA.R](../scripts/04_Restyle_Retinal_MDS_PCA.R) applies the same style to RP MDS and glaucoma MDS/PCA. It aligns annotations using the stored projection's sample IDs, retains the WT/P23H and treatment/time colours, and reads axis percentages from the existing projection objects.

The styling script does not repeat filtering, normalization or projection fitting. Text labels can move for readability; the points remain at their saved coordinates.

## Prerequisites for reproducing the exports

Run from the repository root in an R session with `ggplot2` and `ggrepel` available. The current RP session must contain `rp_mds` and `rp_dge_qc` from the completed RP workflow.

The script uses `glaucoma_mds`, `glaucoma_pca` and `glaucoma_dge_qc` from the same session when present. Otherwise it loads the local `results/glaucoma/Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData` into an isolated environment. It does not load that large workspace into the RP global environment. Checkpoints are not included on GitHub.

For a fresh reconstruction of both datasets, complete the glaucoma precheck and RP precheck through PCA in the same project-specific session before running the styling script. For an existing saved session, deliberately restore the required objects or make the local checkpoint available at the documented path. Uploading this repository did not copy or load either checkpoint.

Both local output directories must already exist, as created by the precheck scripts. Then run:

```r
source(file.path("scripts", "04_Restyle_Retinal_MDS_PCA.R"))
```

This displays and replaces the three named local PNGs under `results/`, at 2400 by 1800 pixels and 300 dpi. The save helper closes each device. Three file-existence checks should return `TRUE`. This export does not save a new R checkpoint or automatically replace the curated image copies under `figures/`.

## LHON specialization

[07_LHON_projection_style.R](../scripts/07_LHON_projection_style.R) uses the same point/repelled-label/outside-legend pattern for the four LHON genotype colours. It defines a formatting function only. The LHON main script calculates both projections and calls the helper before its explicit display/save blocks; this is separate from the RP/glaucoma restyling script. See [the LHON execution guide](running-the-lhon-precheck.md).

## DR microarray specialization

[11_DR_projection_style.R](../scripts/11_DR_projection_style.R) defines matching point/repelled-label/outside-legend formatting for normal, diabetic and treated arrays. It is loaded by the portable DR scripts using the captured project folder. DR projections use the supplied processed microarray values without additional log2, CPM or TMM. MDS axes use the generic names MDS1/MDS2, preserving the documented scale limitation. See the [DR execution guide](running-the-dr-precheck.md).

## Autoimmune uveitis specialization

[14_Uveitis_projection_style.R](../scripts/14_Uveitis_projection_style.R) defines the matching appearance for CON/EAU, with sample labels, unchanged calculated coordinates and an outside condition legend. The main script fits pairwise MDS and common-top-500 PCA on filtered TMM log-CPM, then explicitly displays and saves both plots. The reload helper loads the formatting definition without refitting projections. See the [uveitis execution guide](running-the-uveitis-precheck.md).
