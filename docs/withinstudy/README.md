# Within-study differential-expression analyses

Saved results through 5 November 2026 cover seven terrestrial datasets and 12 predefined comparisons. Each study is analysed against its own controls before comparison with the fixed OSD-758 flight reference. The repository contains full tested-gene tables, FDR-qualified subsets, all-gene rankings, model/plot scripts and unchanged saved figures.

## Recorded comparison results

| Dataset | Comparison | Tested genes | DEGs | Higher | Lower |
| --- | --- | ---: | ---: | ---: | ---: |
| [Glaucoma](glaucoma.md) | 3D versus naive | 23,301 | 6 | 5 | 1 |
| [Glaucoma](glaucoma.md) | 2W versus naive | 23,301 | 2 | 2 | 0 |
| [Glaucoma](glaucoma.md) | 6W versus naive | 23,301 | 3 | 1 | 2 |
| [Retinitis pigmentosa](retinitis_pigmentosa.md) | P23H_vs_WT | 17,351 | 2,576 | 1,361 | 1,215 |
| [LHON-linked model](lhon.md) | YARS2_vs_WT | 17,860 | 3 | 0 | 3 |
| [LHON-linked model](lhon.md) | COI_vs_WT | 17,860 | 1 | 0 | 1 |
| [LHON-linked model](lhon.md) | YARS2_COI_vs_WT | 17,860 | 1,209 | 324 | 885 |
| [Diabetic retinopathy](diabetic_retinopathy.md) | Diabetic_vs_Normal | 24,973 | 3 | 2 | 1 |
| [AMD-related CFH-deficiency model](amd.md) | Young_CfhKO_vs_WT | 19,825 | 3 | 0 | 3 |
| [AMD-related CFH-deficiency model](amd.md) | Aged_CfhKO_vs_WT | 19,825 | 23 | 5 | 18 |
| [Experimental autoimmune uveitis](autoimmune_uveitis.md) | EAU_vs_CON | 18,187 | 5,464 | 2,966 | 2,498 |
| [Light-induced retinal damage](light_damage.md) | LD_vs_CTRL | 16,852 | 2,246 | 911 | 1,335 |

These are BH FDR < 0.05 counts within each contrast, without a separate fold-change threshold. They are copied from existing summaries, rather than recalculated by refitting. Different sample sizes, models, platforms and tested gene universes make counts unsuitable as a disease-severity or spaceflight-resemblance ranking.

## Methods and interpretation

The five RNA-seq analyses reuse their precheck-filtered raw counts and saved TMM factors in the current edgeR 4.8.2 robust quasi-likelihood pipeline. Glaucoma and RP use `~ sex + group`; LHON, uveitis and light damage use their recorded group designs. The direct v4 `glmQLFit` pipeline estimates its required dispersions; an additional legacy dispersion-estimation step is not introduced. See the [official edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf).

DR and AMD use limma linear models with `eBayes(trend = TRUE, robust = TRUE)`. DR uses the supplied processed microarray values and the untreated 6-versus-6 comparison; AMD reuses its precheck log2 intensities after the documented feature-to-gene eligibility/representative selection. Neither uses CPM, TMM or an additional log transform in this stage. The [official limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf) documents these modelling components. AMD's conservative incomplete-annotation rule and representative selection are project-specific choices, documented with saved audits in its report.

Every tested row is retained in a full ranked list. For edgeR, the project-specific score is `sign(logFC) * sqrt(pmax(F, 0))`, with an explicit tolerance guarding negligible negative F values. Original statistics stay intact; the numerical guard does not define DEGs. For limma, the signed moderated `t` statistic is the ranking score. These rankings are inputs for later comparison work, not completed disease-similarity scores or harmonized gene lists.

MD figures show abundance versus estimated change, with the FDR criterion defining colours. Heatmaps use row z-scores; their red/blue colours describe each gene relative to its own mean. RNA-seq dot plots show observed non-log TMM CPM; microarray dots show the respective supplied/log2 expression values. Selecting genes using the same tests means these displays are descriptive, not independent validation. QL/mean–variance diagnostics concern model variability and do not themselves call DEGs or failed samples.

## What is saved

- [Scripts](../../scripts/withinstudy/): 46 source-preserving copies, arranged by dataset, including existing save/reload/verification modules.
- [Tables](../../tables/withinstudy/): full results, DEG subsets, rankings, exact plotted values and available mapping/design/sample audits.
- [Figures](../../figures/withinstudy/): 29 current figures plus AMD's earlier feature-level diagnostic kept as history.
- [Actual session and verification records](records/): existing exports, with trailing whitespace normalized. RP's prepared final checkpoint/session exports were not present; its completed result tables, rankings and figures are included.

AMD's current gene-level results test 19,825 recorded Entrez IDs. Its older feature-level summaries/subsets and scripts remain identifiable as history; two approximately 93 MB complete legacy feature tables stay local. Six final within-study checkpoints were confirmed present without loading them. Existing AMD, uveitis and light-damage verification records describe earlier checks; this upload did not rerun them.

Provenance CSV path fields are adapted to repository-relative local paths; original source checksums remain recorded. AMD's original-run output manifest describes original source filenames/hashes, including local-only files, and is not a checksum manifest for the adapted repository scripts. Original downloaded expression and RData remain outside Git.

See [running or restoring the analyses](running-the-analysis.md). This synchronization reviewed saved outputs and code, checked table consistency, syntax, links and figure identity, and preserved external-drive originals. It did not import expression inputs, load binary models, refit tests, repeat annotation or generate new scientific figures. The completed cross-study comparison uses these ranked tables and the OSD-758 reference. [Final comparison](../main_analysis/README.md) · [Flight reference](osdr758.md)
