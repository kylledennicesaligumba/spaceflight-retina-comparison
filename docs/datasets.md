# Dataset panel and sources

Panel recorded on 26 October 2026; progress updated on 3 November 2026.

The fixed reference is **[OSD-758, mouse retina](https://osdr.nasa.gov/bio/repo/data/studies/OSD-758)**. The planned primary contrast is flight 0 g versus habitat ground controls. [Kremsky et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11593819/) describes the flight experiment.

All seven terrestrial datasets below have completed processed-data prechecks. Recorded sample flags, design-information gaps and later modelling decisions remain documented in their reports; completion does not mean the disease–spaceflight comparison has been performed.

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

This snapshot contains all seven completed terrestrial prechecks and 29 saved figures. RNA-seq prechecks each have four figures; AMD has four and DR has five, including its density close-up. Every selected sample is currently retained. AMD flags GSM947360 for later inclusion review; light damage preserves the control-like LD-vehicle-5. Source authors generated the experimental measurements; this project concerns secondary analysis.

## Measurement and identifier distinctions

| Dataset | Deposited measurement and native row labels | Precheck treatment and annotation |
| --- | --- | --- |
| Glaucoma | RNA-seq counts; combined Ensembl/gene-symbol labels | Expression filter, TMM, log-CPM; preserve supplied labels |
| RP | RNA-seq counts; Ensembl gene IDs | Expression filter, TMM, log-CPM; supplied workbook annotation kept separate from sample columns |
| LHON-linked model | RNA-seq counts; symbol-style gene names | Expression filter, TMM, log-CPM; no external identifier lookup in this precheck |
| DR | Processed microarray; Entrez-based IDs from the custom platform | Supplied values used without additional log2, CPM, TMM or count filtering; supplied platform gene names matched by ID |
| Autoimmune uveitis | RNA-seq counts; mixed symbol-style labels and Ensembl IDs | Expression filter, TMM, log-CPM; BioMart failed, REST annotated the Ensembl subset; unresolved IDs preserved |
| AMD-related CFH-deficiency model | Microarray intensities; numeric array feature/transcript-cluster IDs | One log2 QC copy; no CPM/TMM/extra normalization or expression-row filter; multipart platform assignments preserved separately |
| Light damage | RNA-seq counts; Ensembl gene IDs | Expression filter, TMM, log-CPM; failed REST route replaced by fixed Ensembl release-116 GRCm39 GTF; unnamed/unmatched IDs preserved |

Numeric labels do not imply a shared identifier system: DR's Entrez IDs and AMD's array IDs have different meanings. In AMD, one measured array feature can have several annotation assignments; expanding an annotation table does not create extra expression measurements. Probe groups, gene IDs and symbols have not yet been harmonized across studies. The [AMD report](amd-precheck.md) explains this distinction with examples.

Selection also differs by study: glaucoma uses eight retina libraries across baseline and three stages; RP uses seven P23H/WT libraries; LHON uses 12 libraries in four genotypes; DR explores all 17 arrays while keeping the six-versus-six untreated disease contrast separate from five treated arrays; uveitis excludes both sleep-deprivation groups; AMD excludes RPE/choroid and keeps age-matched genotype comparisons separate; light damage uses all 10 supplied vehicle libraries. These choices do not pool heterogeneous study measurements into one analysis.

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

## LHON-linked model inputs

The source is [GSE304707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304707). The manually reviewed inputs are `GSE304707_gene_count_matrix.csv` and the extensionless local family archive `GSE304707_family`, containing `GSE304707_family.xml`.

The portable default directory is `data/raw/LHON_GSE304707/`; `LHON_DATA_DIR` can point to existing storage. Preserve downloaded filenames, adjusting the archive assignment if your download uses a different name.

The repository includes the final exported matched metadata and TMM table from the completed session. They retain the sample IDs WT-1–3, Y-1–3, C-1–3 and YC-1–3 in count-column order. No original expression input or unreviewed metadata was imported during synchronization. The principal genotype contrast and cross-study gene-name harmonization remain later analysis decisions. See the [LHON report](lhon-precheck.md) and [execution guide](running-the-lhon-precheck.md).

## Diabetic retinopathy inputs

The source is [GSE87433](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87433). The reviewed local filenames are `GSE87433_series_matrix` and `GSE87433_family`; the family archive supplies `GSE87433_family.xml` and `GPL20710-tbl-1.txt`. Preserve source filenames and adapt the script's filename assignments if your download uses different extensions.

This is a processed microarray matrix, not RNA-seq counts. The precheck uses the supplied values without additional log2 transformation, CPM calculation, TMM normalization or count filtering. Entrez-based feature IDs are matched to gene symbols/names in the supplied platform annotation. The default input directory is `data/raw/Diabetic_Retinopathy_GSE87433/`, or the absolute directory set by `DR_DATA_DIR`.

All 17 arrays were reviewed: six normal controls, six untreated diabetic samples and five treated samples. The six-versus-six disease comparison remains distinct from the treated condition. See the [DR report](diabetic-retinopathy-precheck.md) and [execution guide](running-the-dr-precheck.md). The original processed matrix and large platform annotation file are not uploaded as portfolio outputs.

## Autoimmune uveitis inputs

The source is [GSE284505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284505). The reviewed local files are `GSE284505_R.HTSeq.counts.txt` and `GSE284505_family`, an archive containing `GSE284505_family.xml`. The default input directory is `data/raw/Autoimmune_Uveitis_GSE284505/`, or the absolute directory set by `UVEITIS_DATA_DIR`.

The original count table contains 12 samples across CON, EAU, CSD-CON and CSD-EAU. The selected disease comparison retains only CON and EAU, three samples each; exclusion of sleep-deprivation conditions is by design. Mixed symbol-style/Ensembl labels are preserved. The saved annotation records REST lookup results after failed BioMart connections, including missing names and unmatched IDs; it is not a completed cross-study harmonization. See the [uveitis report](autoimmune-uveitis-precheck.md) and [execution guide](running-the-uveitis-precheck.md).

## AMD-related CFH-deficiency model inputs

Use [GSE38671](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38671). Reviewed local inputs are `GSE38671_series_matrix` and `GSE38671_family`; the archive supplies `GSE38671_family.xml` and `GPL6246-tbl-1.txt`. The portable default is `data/raw/AMD_GSE38671/`, or the absolute folder set by `AMD_DATA_DIR`.

Select only the 16 neuroretina arrays from 32 original arrays, preserving four young/aged WT/Cfh-/- groups. Platform coverage is complete for the 35,556 expression IDs, with one extra unused platform record. The numeric feature IDs remain distinct from gene names/Entrez annotations and the separate 298,709-row assignment table. Large platform/expanded annotation exports remain local. See the [AMD report](amd-precheck.md) and [execution guide](running-the-amd-precheck.md).

## Light-damage inputs

Use [GSE314043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314043). Reviewed local inputs are the uncompressed `GSE314043_genecount.txt` and extensionless `GSE314043_family` archive containing `GSE314043_family.xml`. The portable default is `data/raw/Light_Damage_GSE314043/`, or the absolute folder set by `LIGHT_DAMAGE_DATA_DIR`.

All 10 supplied CTRL-vehicle/LD-vehicle libraries are used. Final gene annotations come from the fixed official release-116 GRCm39 GTF, rather than the unsuccessful partial REST cache. The GTF itself remains local; its URL and checksum are versioned alongside saved annotations and QC tables. See the [light-damage report](light-damage-precheck.md) and [execution guide](running-the-light-damage-precheck.md).
