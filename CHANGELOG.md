# Project progress

## 3 November 2026 — All seven terrestrial prechecks complete

- Added AMD and light-damage scripts, reports, execution guides, selected metadata and eight original saved figures.
- Documented AMD's neuroretina-only selection, age/genotype groups, one log2 QC copy, array feature/transcript-cluster IDs, multiple assignments and source truncation. All 16 arrays remain retained; GSM947360 keeps its recorded review flag. No CPM, TMM, extra normalization or expression-row filtering was added.
- Documented light damage's five-versus-five vehicle comparison, Ensembl IDs, fixed release-116 GRCm39 GTF annotation after REST failures, 16,852 retained QC genes and all 10 libraries retained, including control-like LD-vehicle-5. Added existing annotation/QC/projection exports, PCA feature list, session information and completion verification record.
- Replaced the five previously versioned boxplots with their existing group-coloured exports; AMD's newly added boxplot is also restyled and light damage already used group colours. Added the shared restyle script, six sample keys and the saved unchanged-statistics record with portable paths.
- Updated every completed report's boxplot provenance, the README to seven completed prechecks, and the dataset guide with measurement, identifier, transformation and sample-selection distinctions.
- Preserved original files. Confirmed both dated checkpoints exist without loading them; AMD's older checkpoint is not represented as a restore-verified completed state. Large AMD annotation, original expression inputs, the full GTF and binary results remain local.

Synchronization checked script syntax/source consistency, saved tables/keys, links and original PNG identity. It did not import expression inputs, inspect unreviewed source metadata, load RData, rerun analysis or annotation, refit projections or regenerate scientific figures.

## 2 November 2026 — Experimental autoimmune uveitis precheck complete

- Added four portable teaching/reload/style/save scripts, a completion report and an execution guide, with four unchanged saved QC figures.
- Recorded the selected three-CON versus three-EAU comparison; both sleep-deprivation groups are excluded by design, while original full-count objects are preserved.
- Documented the mixed gene-name/Ensembl labels and failed BioMart connections followed by successful Ensembl REST annotation: 1,440 named records, 24 records without names and 84 unresolved IDs; 55,311 other labels were not queried.
- Added the dated full annotation, selected sample metadata, TMM-factor table and actual session record. Preserved original IDs, missing values and release 116/GRCm39 provenance separately from the source GRCm38 processing notes.
- Recorded 18,187 retained QC rows, all six selected samples retained, and completed MDS/PCA review.
- Confirmed the dated final checkpoint is present locally without loading it or testing restoration. Updated the overview to five completed prechecks.

Synchronization preserved original external-drive files and did not import expression inputs, inspect unreviewed metadata, load RData, rerun analysis, make annotation queries or regenerate figures.

## 31 October 2026 — DR processed-microarray precheck complete

- Added portable DR teaching/resume/style scripts, completion report and execution guide, with five unchanged saved QC PNGs.
- Recorded all 17 arrays and 24,973 supplied features used for exploratory QC; the six-versus-six untreated disease comparison remains separate from five treated arrays.
- Made the platform distinction explicit: no additional log2 transformation, CPM calculation, TMM normalization or RNA-seq feature filtering was applied.
- Documented Entrez row IDs and corresponding gene symbols/names from supplied platform annotation.
- Preserved the reviewed boxplot and overlapping density views, including the author's successful check of all 17 curves against their own expression columns.
- Confirmed the new 31 October checkpoint file exists locally without loading it; originals and large inputs/annotation remain outside Git.
- Updated the overview to four completed prechecks.

Synchronization did not import expression inputs, inspect unreviewed metadata, load RData, rerun analysis or regenerate figures.

## 30 October 2026 — LHON-linked model precheck complete

- Added the main LHON teaching script, completed-state reload helper, projection-formatting helper and final save/export helper with portable paths.
- Added the completion report, four original saved QC figures, final matched metadata, TMM-factor export and actual session information.
- Recorded all 12 libraries retained, with 17,860 of 49,714 rows in the QC copy; principal contrast and cross-study gene mapping remain later steps.
- Preserved the original labels, documented timing-information scope, and recorded the existing 47-object checkpoint restore verification. Workspace binaries and raw expression inputs remain local.
- Updated the project overview and dataset/execution documentation to show three completed prechecks.

Synchronization preserved the original MaxineDrive files and did not import expression inputs, load RData, recompute QC or regenerate plots.

## 29 October 2026 — RP precheck complete and projection appearance updated

- Added the completed RP GSE156533 precheck and recovery scripts, report, matched seven-library metadata and four saved QC figures.
- Recorded the retained seven-library comparison and 17,351-gene QC copy; preserved sex and the documented source limitations. Included the existing checkpoint audit and stored TMM-factor table.
- Replaced glaucoma MDS/PCA PNGs with the revised original exports and added the shared styling script. Stored projection scores, percentages and sample-retention decisions were preserved.
- Added portable input/output paths, execution guides and the observed readxl/xml2/ggplot2/ggrepel versions.
- Made the RP final checkpoint base-list initialization explicit in the portable script, using the list from the existing recovery helper; expression calculations were preserved.
- Confirmed the dated RP checkpoint file exists locally and kept both datasets' workspace binaries outside Git.

Original MaxineDrive files were preserved. Synchronization did not load workspaces, import expression inputs, rerun QC or generate figures.

## 28 October 2026 — Completed glaucoma processed-data precheck

- Synced the completed glaucoma script while preserving portable input and output paths.
- Recorded gene-label structural checks, the 23,301-gene filtered QC copy, TMM factors and completed MDS/PCA review.
- Added four unchanged original QC exports: boxplot, density, MDS and common-gene PCA.
- Updated the report and README to record the agreed endpoint: retain all eight retinal libraries, including R-08, and omit the optional sample-correlation check.
- Confirmed the completed 28 October checkpoint file exists locally; kept workspace binaries and downloaded counts outside Git.
- Updated execution instructions and package records to include limma.
- Preserved the existing GitHub edit describing this as a bioinformatics project.

The original external-drive files were preserved. This synchronization did not rerun expression analysis, generate new figures or begin differential-expression/disease-similarity analysis.

## 27 October 2026 — Initial private portfolio snapshot

- Documented the research objective and agreed seven-dataset terrestrial panel against OSD-758.
- Added the current glaucoma precheck script with its educational comments and recorded interactive observations.
- Replaced computer-specific input/output paths in the repository copy with a configurable input folder and relative output folder.
- Added creation of the local checkpoint output directory; retained the original calculations and sample-selection steps.
- Included the saved 36-library metadata table derived from public GEO records.
- Recorded completed library-total and detected-gene checks; left expression-distribution review and subsequent analysis pending.
- Added setup instructions, observed R/package versions and a beginner guide to saving future GitHub updates.

The original working script and saved workspace on the external drive were preserved. No expression analysis was rerun during repository preparation.
