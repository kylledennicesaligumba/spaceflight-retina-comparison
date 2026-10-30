# Local data

Download the original inputs manually using the [dataset source guide](../docs/datasets.md).

The default input folders are `data/raw/Glaucoma_GSE241782/`, `data/raw/Retinitis_Pigmentosa_GSE156533/` and `data/raw/LHON_GSE304707/`. Alternatively, set `GLAUCOMA_DATA_DIR`, `RP_DATA_DIR` or `LHON_DATA_DIR` in R to use existing directories without copying the data. See the [glaucoma guide](../docs/running-the-precheck.md), [RP guide](../docs/running-the-rp-precheck.md) and [LHON guide](../docs/running-the-lhon-precheck.md).

Downloaded files under `data/` are excluded from Git; this README is the exception. Derived sample descriptions suitable for version control are stored under `metadata/`. Generated outputs and R checkpoints go under `results/`, which is also excluded. Selected reviewed QC PNGs are copied to `figures/glaucoma/`, `figures/retinitis_pigmentosa/` and `figures/lhon/` for version control and the completed reports.

This repository documents secondary analysis. Experimental data credit belongs to the original study authors, linked through the public accession records.

The saved LHON normalization summary is versioned under `tables/lhon/`; its original counts and binary workspace remain local.
