# Spaceflight retina and terrestrial disease responses

**Status: glaucoma and retinitis pigmentosa prechecks completed; remaining dataset prechecks in progress.**

Maintained by [Kylle Dennice Saligumba](https://github.com/kylledennicesaligumba).

Last documented update: 29 October 2026.

This bioinformatics project asks: **which selected terrestrial eye-disease models most closely resemble the measured mouse retinal response to spaceflight, and which biological programs are shared or change in opposite directions?**

The fixed spaceflight reference is [NASA OSD-758](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758), associated with [Kremsky et al. (2024)](https://doi.org/10.3390/ijms252212041). The planned primary comparison is flight 0 g versus habitat ground controls. The study collected tissue after the mice returned to Earth, so this contrast represents a flight-associated response with recovery exposures included.

Spaceflight-associated neuro-ocular syndrome (SANS) motivates the research. The measurements are from mouse retina; conclusions will concern resemblance among the tested models. A closest model has not yet been established in this workflow.

## Current work

Two processed-data prechecks are completed. Each includes input integrity, sample matching, gene-label checks, initial expression distributions, expression filtering, TMM scaling, and reviewed MDS/PCA figures.

| Dataset | Completed | Libraries retained | Gene rows in the QC copy | Record |
| --- | --- | --- | --- | --- |
| Experimental glaucoma, GSE241782 | 28 October 2026 | All eight selected retinal libraries | 23,301 of 55,487 | [Glaucoma report](docs/glaucoma-precheck.md) |
| Retinitis pigmentosa, GSE156533 | 29 October 2026 | All seven libraries: four P23H and three WT | 17,351 of 54,532 | [RP report](docs/retinitis-pigmentosa-precheck.md) |

The glaucoma decision retains R-08 with its distinct profile documented and omits the optional correlation check. Glaucoma stages and sex remain recorded; each condition/time pair consists of one male pool and one female pool. The RP decision retains all seven libraries and preserves sex for later genotype modelling. The RP report includes the recorded checkpoint audit and stored TMM-factor table.

These prechecks concern deposited counts and metadata. Raw sequence QC was not repeated, and disease–spaceflight similarity has not yet been calculated. The seven-dataset panel remains the agreed project scope.

## Latest completed precheck: retinitis pigmentosa

The saved RP PCA uses one common set of 500 variable genes after TMM. The genotype and sex-associated patterns are descriptive within-study observations. Separate studies' PCA axes are not a common disease–spaceflight similarity scale.

![RP retinal-library PCA after TMM using a common set of 500 variable genes](figures/retinitis_pigmentosa/RP_GSE156533_TMM_PCA_top500.png)

The [RP report](docs/retinitis-pigmentosa-precheck.md) includes all four figures. The [glaucoma report](docs/glaucoma-precheck.md) now shows the revised MDS/PCA appearance: sample points, readable F/M labels, connector lines and outside legends. Stored coordinates and retention decisions were preserved during that formatting update.

## Explore the repository

- [Glaucoma script](scripts/01_Glaucoma_GSE241782_precheck.R) and [execution guide](docs/running-the-precheck.md).
- [RP script](scripts/02_Retinitis_Pigmentosa_GSE156533_precheck.R), [recovery helper](scripts/03_Retinitis_Pigmentosa_GSE156533_reload.R) and [execution guide](docs/running-the-rp-precheck.md).
- [Shared projection-styling script](scripts/04_Restyle_Retinal_MDS_PCA.R) and [style/reuse instructions](docs/projection-styling.md).
- [Glaucoma metadata](metadata/GSE241782_samples_all.csv) and [RP metadata](metadata/GSE156533_samples_matched.csv), both derived from public source records.
- [Dataset panel and input sources](docs/datasets.md).
- [GitHub first steps](docs/github-first-steps.md).

## Planned comparison

Each study will be assessed against its own controls using a method appropriate to its measurement format. The resulting responses will then be compared against the same spaceflight reference. The objective includes relative resemblance, shared programs and opposing responses; no condition is assumed to be the closest in advance.

The panel contains RNA-seq counts and microarray measurements. Exact statistical models and similarity metrics remain to be finalized after the individual prechecks. Separate studies' raw measurements will not be pooled into one disease/control analysis.

## Reproducibility and scope

The precheck scripts are designed for interactive execution one block at a time. Recorded numerical observations and figures come from the author's interactive R runs. Repository synchronization reviewed those records, checked R syntax and preserved the original figure files without rerunning expression analysis. See [environment details](docs/environment.md).

Original expression inputs remain at their public repositories and on local storage; [data instructions](data/README.md) explain where to obtain them. Downloaded counts, R workspace checkpoints and local files under `results/` are excluded from Git. Selected reviewed QC figures are versioned under `figures/`. This repository contains a portable copy of the working script; changes to a separate working copy need to be brought into this repository deliberately.
