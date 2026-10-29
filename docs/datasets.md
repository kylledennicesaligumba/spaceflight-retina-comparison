# Dataset panel and sources

Panel recorded on 26 October 2026; progress updated on 29 October 2026.

The fixed reference is **[OSD-758, mouse retina](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758)**. The planned primary contrast is flight 0 g versus habitat ground controls. [Kremsky et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11593819/) describes the flight experiment.

The seven terrestrial datasets below are selected for retrieval and individual prechecks. Selection does not mean every dataset has passed quality and experimental-design checks.

| Condition/model | Accession | Measurement | Intended comparison to inspect |
| --- | --- | --- | --- |
| Experimental glaucoma | [GSE241782](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241782) | RNA-seq gene counts | Retinal glaucoma versus naïve controls; retain 3-day, 2-week and 6-week stages |
| Retinitis pigmentosa, P23H rhodopsin model | [GSE156533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156533) | RNA-seq read counts | Mutant versus wild-type littermate retina at postnatal day 90 |
| Diabetic retinal injury, STZ diabetes model | [GSE87433](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87433) | Microarray expression | Untreated diabetic versus normoglycemic retina at the 12-week endpoint |
| LHON-linked mitochondrial optic neuropathy model | [GSE304707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304707) | RNA-seq gene counts | YARS2 G186V, mitochondrial COI V421A and double-mutant groups versus wild type; principal contrast pending |
| AMD-related complement-factor-H deficiency model | [GSE38671](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38671) | Microarray expression | Neural-retina Cfh knockout versus age-matched wild type; keep ages separate |
| Experimental autoimmune uveitis | [GSE284505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284505) | RNA-seq HTSeq counts | EAU-only versus healthy control; exclude sleep-deprivation groups |
| Light-induced retinal damage | [GSE314043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314043) | RNA-seq gene counts | Light-damaged vehicle group versus control vehicle group |

All entries concern mouse retinal measurements. LHON and AMD entries are specific experimental models; light damage is an experimental injury. These labels should be retained when presenting future results.

This snapshot contains completed glaucoma and RP processed-data prechecks, each with four saved QC figures. Glaucoma retains all eight selected libraries, including R-08, and omits the optional correlation check. RP retains all seven libraries and records sex for later genotype modelling. Other datasets' prechecks and the formal disease–spaceflight comparison remain separate project stages. Source authors generated the experimental data; this project concerns secondary analysis. Accession links document the source records, not a completed independent validation of every comparison.

## Glaucoma inputs

Download these two files manually from the GEO record:

| File | Purpose |
| --- | --- |
| [GSE241782_gene_count_matrix.csv.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241782/suppl/GSE241782_gene_count_matrix.csv.gz) | Gene-level counts for the study's 36 libraries |
| [GSE241782_family.soft.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241782/soft/GSE241782_family.soft.gz) | Study design, sample descriptions and supplementary filenames |

Extract each compressed file once. The script expects `GSE241782_gene_count_matrix.csv` and `GSE241782_family.soft`. Preserve source filenames. The [execution guide](running-the-precheck.md) describes local input paths.

The included `metadata/GSE241782_samples_all.csv` is the saved, parsed 36-library table from the interactive workflow. It retains original GEO titles and filenames alongside extracted fields. It is a derived metadata table, not a newly generated biological measurement.

## Retinitis pigmentosa inputs

The source is [GSE156533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156533), associated with [Leinonen et al. (2020)](https://elifesciences.org/articles/59422).

| Input | Purpose |
| --- | --- |
| [GSE156533_readcounts_organized_by_largest_to_smallest_count.xlsx](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156533/suppl/GSE156533_readcounts_organized_by_largest_to_smallest_count.xlsx) | Gene counts and supplied annotations; only the seven sample columns enter the count matrix |
| GEO MINiML family archive from the series record, locally saved as `GSE156533_family` | Study/sample XML metadata; this extensionless local file is a tar archive containing `GSE156533_family.xml` |

Preserve downloaded filenames. If the archive downloads under a different name, point `rp_metadata_archive` to that file rather than treating it as SOFT text. The portable default directory is `data/raw/Retinitis_Pigmentosa_GSE156533/`; `RP_DATA_DIR` can point to existing storage.

The included `metadata/GSE156533_samples_matched.csv` was extracted from the existing public XML during repository synchronization. Its seven rows follow the count-column order recorded in the completed precheck, preserving accessions, titles, genotype, sex, age, tissue and strain. The expression workbook was not imported during synchronization. See the [RP execution guide](running-the-rp-precheck.md).
