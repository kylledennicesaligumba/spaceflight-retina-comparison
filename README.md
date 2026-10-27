# Spaceflight retina and terrestrial disease responses

**Status: preprocessing and quality checks in progress.**

Maintained by [Kylle Dennice Saligumba](https://github.com/kylledennicesaligumba).

Last documented update: 27 October 2026.

This bioinformatics project asks: **which selected terrestrial eye-disease models most closely resemble the measured mouse retinal response to spaceflight, and which biological programs are shared or change in opposite directions?**

The fixed spaceflight reference is [NASA OSD-758](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758), associated with [Kremsky et al. (2024)](https://doi.org/10.3390/ijms252212041). The planned primary comparison is flight 0 g versus habitat ground controls. The study collected tissue after the mice returned to Earth, so this contrast represents a flight-associated response with recovery exposures included.

Spaceflight-associated neuro-ocular syndrome (SANS) motivates the research. The measurements are from mouse retina; conclusions will concern resemblance among the tested models. A closest model has not yet been established in this workflow.

## Current work

The first documented dataset is **GSE241782, experimental glaucoma**. Its precheck script covers importing gene counts, matching sample descriptions, checking count integrity, selecting retinal libraries, and inspecting initial library summaries.

| Component | Current status |
| --- | --- |
| Dataset panel | Seven terrestrial condition/model datasets selected for individual prechecks |
| Glaucoma counts and metadata | Imported and matched in the interactive workflow |
| Retinal glaucoma/control selection | 55,487 gene rows and eight pooled libraries; sample order checked |
| Library totals and detected genes | Interactive results recorded |
| Initial expression distributions | Plotting code included; plot review pending in this snapshot |
| Remaining QC, disease signatures and cross-study comparisons | Planned |

The glaucoma libraries comprise two naïve controls and two libraries at each of three glaucoma stages. Each pair contains one male pool and one female pool. Retaining stages and sex is important: six glaucoma libraries across three stages are not six replicates at one stage. See the [precheck record](docs/glaucoma-precheck.md) and the [original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531004/).

## Explore the repository

- [Glaucoma precheck script](scripts/01_Glaucoma_GSE241782_precheck.R): the current educational script, adapted to portable paths.
- [Sample metadata](metadata/GSE241782_samples_all.csv): 36 library records parsed from public GEO metadata, before retinal selection.
- [Dataset panel and sources](docs/datasets.md): accessions, measurement types and intended comparisons.
- [Precheck observations](docs/glaucoma-precheck.md): completed checks, observed summaries and outstanding steps.
- [Running the precheck](docs/running-the-precheck.md): inputs, R requirements and block-by-block instructions.
- [GitHub first steps](docs/github-first-steps.md): what commits and pushes mean, and how to save future updates.

## Planned comparison

Each study will be assessed against its own controls using a method appropriate to its measurement format. The resulting responses will then be compared against the same spaceflight reference. The objective includes relative resemblance, shared programs and opposing responses; no condition is assumed to be the closest in advance.

The panel contains RNA-seq counts and microarray measurements. Exact statistical models and similarity metrics remain to be finalized after the individual prechecks. Separate studies' raw measurements will not be pooled into one disease/control analysis.

## Reproducibility and scope

The script is designed for interactive execution one block at a time. Recorded numerical observations come from the author's interactive R runs. Repository setup checked R syntax and file organization without rerunning expression analysis. See [environment details](docs/environment.md).

Original expression inputs remain at their public repositories and on local storage; [data instructions](data/README.md) explain where to obtain them. R workspace checkpoints and generated results are excluded from Git. This repository contains a portable copy of the working script; changes to a separate working copy need to be brought into this repository deliberately.
