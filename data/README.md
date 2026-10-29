# Local data

Download the original inputs manually using the [dataset source guide](../docs/datasets.md).

The default input folders are `data/raw/Glaucoma_GSE241782/` and `data/raw/Retinitis_Pigmentosa_GSE156533/`. Alternatively, set `GLAUCOMA_DATA_DIR` or `RP_DATA_DIR` in R to use existing directories without copying the data. See the [glaucoma guide](../docs/running-the-precheck.md) and [RP guide](../docs/running-the-rp-precheck.md).

Downloaded files under `data/` are excluded from Git; this README is the exception. Derived sample descriptions suitable for version control are stored under `metadata/`. Generated outputs and R checkpoints go under `results/`, which is also excluded. Selected reviewed QC PNGs are copied to `figures/glaucoma/` and `figures/retinitis_pigmentosa/` for version control and the completed reports.

This repository documents secondary analysis. Experimental data credit belongs to the original study authors, linked through the public accession records.
