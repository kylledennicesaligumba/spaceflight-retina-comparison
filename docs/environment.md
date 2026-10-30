# Observed R environment

Checked on 29 October 2026 without loading expression inputs or saved R workspaces.

| Component | Observed value |
| --- | --- |
| R | 4.5.3 (2026-03-11) |
| Platform | aarch64-apple-darwin20 |
| edgeR | 4.8.2 |
| limma | 3.66.0 |
| readxl | 1.5.0 |
| xml2 | 1.5.2 |
| ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 |
| BiocManager | 1.30.27 |

Base R handles import, metadata matching, initial count summaries, gene-label checks and initial distribution plots. `edgeR` supplies `cpm()`, `DGEList()`, `filterByExpr()` and `normLibSizes()`. The exported `limma::plotMDS()` generic dispatches to the registered DGEList method for the pairwise MDS and common-gene PCA views. BiocManager is an installation utility, not a calculation dependency of the script.

RP uses `readxl` to read the deposited count workbook and `xml2` to parse GEO sample metadata. Its MDS/PCA calls take the already prepared TMM log-CPM matrix. The revised projection appearance uses `ggplot2` and `ggrepel` with stored scores.

This is a record of the available environment, not a dependency lockfile or a claim that the full workflow has been rerun on other systems. Record `sessionInfo()` alongside future analysis outputs when those steps are executed.

## Repository preparation checks

- Parsed the R script without executing its expressions.
- Compared parsed calculation expressions with the original working scripts, allowing portable file-path changes, output-directory creation and explicit RP checkpoint-list initialization.
- Checked the included metadata table's structure, unique identifiers and selected sample labels.
- Checked local Markdown file links and the list of files selected for Git.
- Verified that raw data, R workspace files and generated results are ignored by Git.
- Inspected the revised glaucoma projections and all four RP figures; checked that the uploaded image copies match the original files.

The 28 October glaucoma and 29 October RP checkpoint files were confirmed on local storage; their binary contents were not loaded during synchronization. The small RP sample table was extracted from the already available public XML; the expression workbook was not read during this update.

No expression analysis was run to prepare this repository. Recorded numerical observations retain their provenance in the [glaucoma record](glaucoma-precheck.md) and [RP record](retinitis-pigmentosa-precheck.md).

## LHON completed-session record — 30 October 2026

The [exported sessionInfo](lhon-sessionInfo_2026-10-30.txt) records the actual R session used for the completed LHON precheck. It includes R 4.5.3, edgeR 4.8.2, limma 3.66.0, xml2 1.5.2, ggplot2 4.0.3 and ggrepel 0.9.8. System-library paths are retained; trailing whitespace was normalized for version control.

The LHON update copies the final metadata and TMM exports, original figures and portable scripts. Syntax and source-expression comparison are performed without running the workflow. The existing completion record reports a restore-verified 47-object checkpoint; the GitHub update confirms its file presence but does not load it. See [the completed LHON report](lhon-precheck.md).
