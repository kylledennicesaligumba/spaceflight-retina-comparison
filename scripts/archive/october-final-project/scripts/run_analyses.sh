#!/bin/zsh
set -euo pipefail
TASK_ROOT="$(cd -- "$(dirname -- "$0")/.." && pwd)"
cd -- "$TASK_ROOT"
/usr/local/bin/Rscript scripts/analyze.R
/usr/local/bin/Rscript scripts/separate_glaucoma_crush.R
/usr/local/bin/Rscript scripts/sensitivity.R
/usr/local/bin/Rscript scripts/alternative_glaucoma.R
/usr/local/bin/Rscript scripts/alternative_glaucoma_date_sensitivity.R
/usr/local/bin/Rscript scripts/alternative_uveitis.R
/usr/local/bin/Rscript scripts/shared_universe.R
/usr/local/bin/Rscript scripts/pca_overview.R
