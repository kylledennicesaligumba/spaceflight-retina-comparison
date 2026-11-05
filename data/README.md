# Local data

Download the original inputs manually using the [dataset source guide](../docs/datasets.md).

Use these defaults, or set the corresponding R environment variable to an existing absolute input directory:

| Dataset | Default folder | Override |
| --- | --- | --- |
| Glaucoma | `data/raw/Glaucoma_GSE241782/` | `GLAUCOMA_DATA_DIR` |
| RP | `data/raw/Retinitis_Pigmentosa_GSE156533/` | `RP_DATA_DIR` |
| LHON-linked model | `data/raw/LHON_GSE304707/` | `LHON_DATA_DIR` |
| DR | `data/raw/Diabetic_Retinopathy_GSE87433/` | `DR_DATA_DIR` |
| Autoimmune uveitis | `data/raw/Autoimmune_Uveitis_GSE284505/` | `UVEITIS_DATA_DIR` |
| AMD-related CFH-deficiency model | `data/raw/AMD_GSE38671/` | `AMD_DATA_DIR` |
| Light damage | `data/raw/Light_Damage_GSE314043/` | `LIGHT_DAMAGE_DATA_DIR` |

See the [glaucoma guide](../docs/running-the-precheck.md), [RP guide](../docs/running-the-rp-precheck.md), [LHON guide](../docs/running-the-lhon-precheck.md), [DR guide](../docs/running-the-dr-precheck.md) and [uveitis guide](../docs/running-the-uveitis-precheck.md).

The final two guides cover [AMD](../docs/running-the-amd-precheck.md) and [light damage](../docs/running-the-light-damage-precheck.md). All seven terrestrial prechecks are documented. See the [measurement/identifier distinctions](../docs/datasets.md#measurement-and-identifier-distinctions).

Downloaded files under `data/` are excluded from Git; this README is the exception. Derived sample descriptions suitable for version control are stored under `metadata/`. Generated outputs and R checkpoints go under `results/`, which is also excluded. Selected reviewed QC PNGs are copied to `figures/glaucoma/`, `figures/retinitis_pigmentosa/`, `figures/lhon/` and `figures/diabetic_retinopathy/` for version control and the completed reports.

The completed uveitis exports are under `figures/autoimmune_uveitis/`, `metadata/GSE284505_samples_selected.csv` and `tables/autoimmune_uveitis/`. The saved annotation preserves original mixed labels and missing lookup results; it contains no expression measurements. Original counts, XML inputs and RData stay local.

This repository documents secondary analysis. Experimental data credit belongs to the original study authors, linked through the public accession records.

The saved LHON normalization summary is versioned under `tables/lhon/`; its original counts and binary workspace remain local.

DR uses a processed microarray matrix and supplied Entrez-based platform annotation. Its precheck adds no log2 transformation, CPM or TMM; original input/annotation files and R checkpoints remain local.

AMD uses continuous microarray intensities and array feature IDs. A separate QC copy is transformed once with log2; no CPM, TMM or extra normalization is applied. The 88.8 MB platform table, 1.82 GB expanded assignment export and 26.45 MB viewing copy stay local. Selected metadata and four reviewed figures are versioned; the formatter preserves missing/multiple/truncated annotations without expanding expression rows.

Light-damage results include selected sample metadata, readable annotation and QC/projection summaries under `tables/light_damage/`, with four figures under `figures/light_damage/`. Original gene-count CSV/RDS, normalized-expression RDS, the full fixed GTF and workspace are local outputs/references and are not uploaded.

The shared [boxplot restyler](../docs/boxplot-styling.md) reads six local dated checkpoints under `results/` by default, or an existing output tree selected by `PREQUALITYCHECKS_DIR`. Its saved sample keys and export record are versioned separately from those workspaces.

Within-study models start from the completed precheck checkpoints; they do not redownload measurements or repeat precheck filtering/normalization. [Within-study execution instructions](../docs/withinstudy/running-the-analysis.md) list the seven checkpoint override variables. Local model outputs/checkpoints go to `results/withinstudy/<dataset>/`, which remains excluded from Git. Selected result/DEG/ranking CSVs and plot values are curated under `tables/withinstudy/`, original figure copies under `figures/withinstudy/`, and actual session/verification exports under `docs/withinstudy/records/`. Small DEG/heatmap expression extracts describe the saved result displays; original complete expression inputs and binary models remain local.
