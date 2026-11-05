# Running or restoring within-study analyses

The [overview](README.md) describes the saved results. This upload does not require rerunning analysis: figures and complete/DEG/ranked tables are already versioned. The scripts reconstruct the recorded work when deliberately executed with its local inputs.

## Repository layout and inputs

Start a project-specific R session in the repository root. The portable scripts capture `portfolio_repo_dir` once, before working-directory changes, or use the absolute root set by `RETINA_PROJECT_DIR`. Scripts live under `scripts/withinstudy/<dataset>/`; generated local results go under `results/withinstudy/<dataset>/`. Versioned tables/figures are curated copies under `tables/withinstudy/` and `figures/withinstudy/` and are not replaced automatically.

Each model starts from its completed precheck checkpoint. Those binaries contain the required existing filtered/TMM objects or processed microarray values/annotations and are not included on GitHub. Default locations are under `results/<dataset>/`; use an environment variable to point directly to an existing checkpoint:

| Dataset folder | Input override | Recorded precheck checkpoint |
| --- | --- | --- |
| `glaucoma` | `GLAUCOMA_PRECHECK_FILE` | `Glaucoma_GSE241782_precheck_checkpoint_2026-10-28.RData` |
| `retinitis_pigmentosa` | `RP_PRECHECK_FILE` | `RP_GSE156533_precheck_checkpoint_2026-10-29.RData` |
| `lhon` | `LHON_PRECHECK_FILE` | `LHON_GSE304707_precheck_checkpoint_2026-10-30.RData` |
| `diabetic_retinopathy` | `DR_PRECHECK_FILE` | `DR_GSE87433_precheck_checkpoint_2026-10-31.RData` |
| `amd` | `AMD_PRECHECK_FILE` | `AMD_GSE38671_precheck_checkpoint_2026-11-03_130914.RData` |
| `autoimmune_uveitis` | `UVEITIS_PRECHECK_FILE` | `Autoimmune_Uveitis_GSE284505_precheck_checkpoint_2026-11-02.RData` |
| `light_damage` | `LIGHT_DAMAGE_PRECHECK_FILE` | `Light_Damage_GSE314043_precheck_checkpoint_2026-11-03.RData` |

```r
# Point to an existing completed precheck when it is stored outside the repository.
Sys.setenv(RP_PRECHECK_FILE = "/path/to/RP_GSE156533_precheck_checkpoint_2026-10-29.RData")
```

The recorded timestamped AMD precheck is newer than its earlier partial 2 November checkpoint and includes the annotation/log2 state needed for the final mapping. If rebuilding prechecks creates a different timestamped filename, point `AMD_PRECHECK_FILE` to that completed file rather than substituting the earlier partial state.

Use the actual session records for version context: the count models use edgeR 4.8.2, with limma 3.66.0 for microarray modelling/related helpers. Plot modules also use packages such as `ggplot2`, `ggrepel`, `pheatmap` and `statmod`. No installation is performed by these copies. Some scripts explicitly check the recorded edgeR version; running another version is not an audited reproduction of these outputs.

## Main scripts and final reload helpers

| Dataset | Main entry point | Final reload helper |
| --- | --- | --- |
| Glaucoma | [01_Glaucoma_GSE241782_withinstudy.R](../../scripts/withinstudy/glaucoma/01_Glaucoma_GSE241782_withinstudy.R) | [10_Glaucoma_GSE241782_reload_final.R](../../scripts/withinstudy/glaucoma/10_Glaucoma_GSE241782_reload_final.R) |
| RP | [01_RP_GSE156533_withinstudy.R](../../scripts/withinstudy/retinitis_pigmentosa/01_RP_GSE156533_withinstudy.R) | [07_RP_GSE156533_reload_final.R](../../scripts/withinstudy/retinitis_pigmentosa/07_RP_GSE156533_reload_final.R) |
| LHON | [01_LHON_GSE304707_withinstudy.R](../../scripts/withinstudy/lhon/01_LHON_GSE304707_withinstudy.R) | [07_LHON_GSE304707_reload_final.R](../../scripts/withinstudy/lhon/07_LHON_GSE304707_reload_final.R) |
| DR | [01_DR_GSE87433_withinstudy.R](../../scripts/withinstudy/diabetic_retinopathy/01_DR_GSE87433_withinstudy.R) | [08_DR_GSE87433_reload_final.R](../../scripts/withinstudy/diabetic_retinopathy/08_DR_GSE87433_reload_final.R) |
| AMD gene-level | [00_AMD_GSE38671_run_gene_level.R](../../scripts/withinstudy/amd/00_AMD_GSE38671_run_gene_level.R) | [06_AMD_GSE38671_reload_final.R](../../scripts/withinstudy/amd/06_AMD_GSE38671_reload_final.R) |
| Uveitis | [01_Autoimmune_Uveitis_GSE284505_withinstudy.R](../../scripts/withinstudy/autoimmune_uveitis/01_Autoimmune_Uveitis_GSE284505_withinstudy.R) | [02_Autoimmune_Uveitis_GSE284505_reload_final.R](../../scripts/withinstudy/autoimmune_uveitis/02_Autoimmune_Uveitis_GSE284505_reload_final.R) |
| Light damage | [01_Light_Damage_GSE314043_withinstudy.R](../../scripts/withinstudy/light_damage/01_Light_Damage_GSE314043_withinstudy.R) | [02_Light_Damage_GSE314043_reload_final.R](../../scripts/withinstudy/light_damage/02_Light_Damage_GSE314043_reload_final.R) |

Glaucoma, RP, LHON and DR retain the manual teaching sequence: review explained blocks and outputs in order. Historical pending comments describe their original stage; use the saved result reports for the upload's current state. Main scripts eventually source their plot/rank/save helpers. AMD scripts 01–02 preserve the earlier feature-level analysis; its current gene-level entry point is script 00, which sources 03–05. Uveitis/light-damage main scripts are the existing complete reproducible workflows. Running a full entry point repeats fitting and writes its named local outputs.

```r
# Deliberately reconstruct the recorded AMD gene-level workflow from its precheck.
source(file.path("scripts", "withinstudy", "amd", "00_AMD_GSE38671_run_gene_level.R"))
```

## Restoring and verifying saved results

Make the required final checkpoint available in its local `results/withinstudy/<dataset>/` folder, then use the relevant final reload helper. Glaucoma's intermediate fit/results reload scripts are retained for history; they are not the final checkpoint entry point. RP's final checkpoint/session export were not present during this upload, despite the existing completed tables/figures; its save helper must run in a session containing the completed RP objects before that final reload can work.

```r
# Restore an existing final checkpoint rather than repeating model fitting.
source(file.path("scripts", "withinstudy", "glaucoma", "10_Glaucoma_GSE241782_reload_final.R"))
```

Original checkpoints retain the path values used when created. Portable reloads reset the main within-study output directory after global loads, but other saved input/output filename objects can still contain original paths. Keep original inputs/result trees accessible, or deliberately relocate/reset the relevant saved paths before new file access or writes. Verification scripts may rely on the original input checkpoint path saved inside the final workspace.

Existing AMD, uveitis and light-damage verification scripts restore saved objects and compare numerical outputs/plot values; running them requires the complete local checkpoint tree. They were not executed for this upload. Their saved verification records are included as evidence of checks performed when the analyses were completed. No new annotation retrieval is part of the within-study upload.

Portable changes are limited to root capture, local-directory creation, input/output paths, helper-script locations and restoring the main output path after global loads. All 46 scripts were syntax-parsed, and reversing the documented adaptations recovered the exact original source expressions. Models, contrasts, selection rules, tests, ranking definitions and figure calculations were preserved.
