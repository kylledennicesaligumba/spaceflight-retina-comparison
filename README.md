# Spaceflight retina and terrestrial disease responses

**Status: dataset prechecks in progress; glaucoma precheck completed.**

Maintained by [Kylle Dennice Saligumba](https://github.com/kylledennicesaligumba).

Last documented update: 28 October 2026.

This bioinformatics project asks: **which selected terrestrial eye-disease models most closely resemble the measured mouse retinal response to spaceflight, and which biological programs are shared or change in opposite directions?**

The fixed spaceflight reference is [NASA OSD-758](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758), associated with [Kremsky et al. (2024)](https://doi.org/10.3390/ijms252212041). The planned primary comparison is flight 0 g versus habitat ground controls. The study collected tissue after the mice returned to Earth, so this contrast represents a flight-associated response with recovery exposures included.

Spaceflight-associated neuro-ocular syndrome (SANS) motivates the research. The measurements are from mouse retina; conclusions will concern resemblance among the tested models. A closest model has not yet been established in this workflow.

## Current work

The processed-data precheck for **GSE241782, experimental glaucoma**, was completed on **28 October 2026**. It covers input integrity, sample matching and selection, gene-label structure, expression distributions, expression filtering, TMM composition normalization, and exploratory MDS/PCA plots.

| Component | Current status |
| --- | --- |
| Dataset panel | Seven terrestrial condition/model datasets selected for individual prechecks |
| Glaucoma counts and metadata | Imported and matched in the interactive workflow |
| Retinal glaucoma/control selection | 55,487 gene rows and eight pooled libraries; sample order checked |
| Library totals, detection, boxplots and density plots | Reviewed and documented |
| Gene-label structure and row alignment | Checked; cross-study reference harmonization remains a later step |
| Expression filtering and TMM | 23,301 genes retained in the QC copy; all eight libraries retained |
| MDS and PCA | Completed, reviewed and saved; four QC figures included |
| Glaucoma precheck decision | Retain all eight libraries, including R-08; optional sample correlation omitted |
| Other dataset prechecks, disease signatures and cross-study comparisons | In progress / planned |

The glaucoma libraries comprise two naïve controls and two libraries at each of three glaucoma stages. Each pair contains one male pool and one female pool. Retaining stages and sex is important: six glaucoma libraries across three stages are not six replicates at one stage. See the [precheck record](docs/glaucoma-precheck.md) and the [original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531004/).

R-08's distinct profile remains documented without classifying it as a failed library. The completed precheck concerns deposited processed counts and metadata; raw sequencing QC was not repeated.

## Glaucoma QC preview

The saved PCA shows the eight retained libraries after TMM, using one common set of 500 variable genes. The displayed axes describe these expression profiles; they are not disease–spaceflight similarity scores.

![Glaucoma retinal-library PCA after TMM using a common set of 500 variable genes](figures/glaucoma/Glaucoma_GSE241782_TMM_PCA_top500.png)

See the [completed report and all four figures](docs/glaucoma-precheck.md) for the observations, sample-retention decision and design limitations.

## Explore the repository

- [Glaucoma precheck script](scripts/01_Glaucoma_GSE241782_precheck.R): the current educational script, adapted to portable paths.
- [Sample metadata](metadata/GSE241782_samples_all.csv): 36 library records parsed from public GEO metadata, before retinal selection.
- [Dataset panel and sources](docs/datasets.md): accessions, measurement types and intended comparisons.
- [Completed glaucoma precheck](docs/glaucoma-precheck.md): checks, observed summaries, four figures and the retained-library decision.
- [Running the precheck](docs/running-the-precheck.md): inputs, R requirements and block-by-block instructions.
- [GitHub first steps](docs/github-first-steps.md): what commits and pushes mean, and how to save future updates.

## Planned comparison

Each study will be assessed against its own controls using a method appropriate to its measurement format. The resulting responses will then be compared against the same spaceflight reference. The objective includes relative resemblance, shared programs and opposing responses; no condition is assumed to be the closest in advance.

The panel contains RNA-seq counts and microarray measurements. Exact statistical models and similarity metrics remain to be finalized after the individual prechecks. Separate studies' raw measurements will not be pooled into one disease/control analysis.

## Reproducibility and scope

The script is designed for interactive execution one block at a time. Recorded numerical observations and figures come from the author's interactive R runs. Repository synchronization reviewed those records, checked R syntax and preserved the original figure files without rerunning expression analysis. See [environment details](docs/environment.md).

Original expression inputs remain at their public repositories and on local storage; [data instructions](data/README.md) explain where to obtain them. Downloaded counts, R workspace checkpoints and local files under `results/` are excluded from Git. Selected reviewed QC figures are versioned under `figures/`. This repository contains a portable copy of the working script; changes to a separate working copy need to be brought into this repository deliberately.
