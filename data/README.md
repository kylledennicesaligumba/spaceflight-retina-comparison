# Local data

Download the original inputs manually using the [dataset source guide](../docs/datasets.md).

Use these defaults, or set the corresponding R environment variable to an existing absolute input directory:

| Dataset | Default folder | Override |
| --- | --- | --- |
| Glaucoma | `data/raw/Glaucoma_GSE241782/` | `GLAUCOMA_DATA_DIR` |
| RP | `data/raw/Retinitis_Pigmentosa_GSE156533/` | `RP_DATA_DIR` |
| LHON-linked model | `data/raw/LHON_GSE304707/` | `LHON_DATA_DIR` |
| DR | `data/raw/Diabetic_Retinopathy_GSE87433/` | `DR_DATA_DIR` |

See the [glaucoma guide](../docs/running-the-precheck.md), [RP guide](../docs/running-the-rp-precheck.md), [LHON guide](../docs/running-the-lhon-precheck.md) and [DR guide](../docs/running-the-dr-precheck.md).

Downloaded files under `data/` are excluded from Git; this README is the exception. Derived sample descriptions suitable for version control are stored under `metadata/`. Generated outputs and R checkpoints go under `results/`, which is also excluded. Selected reviewed QC PNGs are copied to `figures/glaucoma/`, `figures/retinitis_pigmentosa/`, `figures/lhon/` and `figures/diabetic_retinopathy/` for version control and the completed reports.

This repository documents secondary analysis. Experimental data credit belongs to the original study authors, linked through the public accession records.

The saved LHON normalization summary is versioned under `tables/lhon/`; its original counts and binary workspace remain local.

DR uses a processed microarray matrix and supplied Entrez-based platform annotation. Its precheck adds no log2 transformation, CPM or TMM; original input/annotation files and R checkpoints remain local.
