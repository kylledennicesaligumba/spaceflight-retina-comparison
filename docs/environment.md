# Observed R environment

Checked on 27 October 2026 without loading expression inputs or saved R workspaces.

| Component | Observed value |
| --- | --- |
| R | 4.5.3 (2026-03-11) |
| Platform | aarch64-apple-darwin20 |
| edgeR | 4.8.2 |
| BiocManager | 1.30.27 |

Base R handles the import, metadata matching and initial count summaries. `edgeR` supplies `cpm()` for Step 9C. BiocManager is an installation utility, not a calculation dependency of the script.

This is a record of the available environment, not a dependency lockfile or a claim that the full workflow has been rerun on other systems. Record `sessionInfo()` alongside future analysis outputs when those steps are executed.

## Repository preparation checks

- Parsed the R script without executing its expressions.
- Checked the included metadata table's structure, unique identifiers and selected sample labels.
- Checked local Markdown file links and the list of files selected for Git.
- Verified that raw data, R workspace files and generated results are ignored by Git.

No expression analysis was run to prepare this repository. Recorded numerical observations retain their provenance in [the precheck record](glaucoma-precheck.md).
