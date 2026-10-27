# Glaucoma GSE241782: precheck record

Status recorded on 27 October 2026. This record concerns the current interactive precheck, rather than earlier exploratory analyses in other folders.

## Objective

Identify the retinal glaucoma and control libraries, verify the input structure and sample correspondence, and inspect initial expression summaries before choosing the final analysis design.

The original study is [Keuthan et al. (2023), Regional Gene Expression in the Retina, Optic Nerve Head, and Optic Nerve of Mice with Optic Nerve Crush and Experimental Glaucoma](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531004/).

## Completed checks

- The source table contains 55,487 gene rows and 36 library columns, plus a gene-identifier column.
- Recorded integrity checks found no missing values or duplicate complete gene identifiers. The sample values were numeric, finite, nonnegative and whole-valued.
- All 36 sample identifiers matched the GEO metadata uniquely. Metadata rows were reordered to match count columns, and exact identity/order checks returned `TRUE`.
- The study contains three tissue categories. Retinal selection was made using the explicit `Retina` label.
- The eight selected glaucoma/control libraries retain treatment, collection time and sex. The selected count matrix has 55,487 rows and eight columns, with exact metadata alignment confirmed.

These observations were recorded from interactive R output during 26–27 October. They were not independently recomputed during repository preparation. No gene filtering has been applied at this point.

## Selected libraries and observed summaries

Here a **library** means one sequenced RNA preparation. **Total counts** are the sum of the gene-matrix values within that library. **Detected gene rows** means rows whose count is greater than zero; this is an initial summary, not the final expression-filtering rule.

| Library | Treatment | Time | Sex | Total counts | Detected gene rows |
| --- | --- | --- | --- | ---: | ---: |
| R-03 | Glaucoma | 2W | Male | 40,946,004 | 28,154 |
| R-04 | Glaucoma | 2W | Female | 32,514,596 | 27,497 |
| R-07 | Glaucoma | 3D | Male | 46,940,497 | 28,191 |
| R-08 | Glaucoma | 3D | Female | 19,931,516 | 25,647 |
| R-09 | Glaucoma | 6W | Male | 26,597,848 | 26,320 |
| R-10 | Glaucoma | 6W | Female | 34,211,389 | 28,001 |
| R-11 | Naïve control | 0D | Female | 43,795,226 | 29,611 |
| R-12 | Naïve control | 0D | Male | 37,153,668 | 27,855 |

R-08 has the lowest total counts and detected-row count. These summaries alone do not justify excluding it. All eight libraries are retained for the remaining QC checks.

## Pooling and replication

The [original paper, Methods 4.8–4.9](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531004/), and the GEO protocols describe RNA pooled within sex. Four original RNA samples per tissue/treatment/time group produced two measured libraries: one male pool and one female pool.

The count matrix therefore contains two observed libraries per group. The original samples within a pool cannot be treated as separately measured expression profiles. This small number of libraries limits precision and needs to inform later statistical modeling. The three glaucoma stages share a two-library naïve baseline in the deposited comparison.

## Current stopping point

Step 9C includes an initial boxplot using `edgeR::cpm(..., log = TRUE, prior.count = 2)`. This scales by library totals for a log-CPM view; it does not compute TMM composition-normalization factors. Plot review remains pending in this snapshot.

Next work is to review expression distributions, then continue gene-identifier/annotation and normalized sample-relationship checks. No final sample exclusions, analysis-eligibility decision, differential-expression results or disease-similarity ranking have been established in this precheck.
