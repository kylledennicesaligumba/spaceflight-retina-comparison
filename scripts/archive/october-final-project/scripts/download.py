from pathlib import Path
import urllib.request, json, csv, concurrent.futures, hashlib, datetime

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / 'downloads'; D.mkdir(parents=True, exist_ok=True)
O = ROOT / 'outputs'; O.mkdir(exist_ok=True)
tasks = []
for acc in ['GSE241782','GSE156533','GSE284505','GSE87433','GSE314043','GSE26299','GSE144168']:
    tasks.append((acc+'_family.soft', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc+'&targ=all&view=full&form=text'))
tasks += [
 ('GSE241782_gene_count_matrix.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241782/suppl/GSE241782_gene_count_matrix.csv.gz'),
 ('GSE156533_counts.xlsx','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156533/suppl/GSE156533_readcounts_organized_by_largest_to_smallest_count.xlsx'),
 ('GSE284505_counts.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284505/suppl/GSE284505_R.HTSeq.counts.txt.gz'),
 ('GSE314043_counts.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314043/suppl/GSE314043_genecount.txt.gz'),
 ('GSE87433_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87433/matrix/GSE87433_series_matrix.txt.gz'),
 ('GSE26299_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26299/matrix/GSE26299_series_matrix.txt.gz'),
 ('GSE144168_normalized_CPM.xlsx','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144168/suppl/GSE144168_Normalized_counts_CPM_Lipski_171130.xlsx'),
 ('OSD-758_metadata.json','https://osdr.nasa.gov/osdr/data/osd/meta/758'),
 ('OSD-758_files.json','https://osdr.nasa.gov/osdr/data/osd/files/758'),
 ('OSD-758_assays.json','https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-758/assay/*/sample/*/'),
 ('PMC11593819.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11593819/fullTextXML'),
 ('PMC10531004.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10531004/fullTextXML'),
 ('PMC7529457.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7529457/fullTextXML'),
 ('PMC7076950.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7076950/fullTextXML'),
 ('diabetes_original.pdf','https://link.springer.com/content/pdf/10.1007/s00125-017-4254-y.pdf'),
 ('JCI44646.html','https://www.jci.org/articles/view/44646'),
 ('JCI44646_supplement.pdf','https://www.jci.org/articles/view/44646/sd/pdf/render/1'),
 ('EAU_original_crossref.json','https://api.crossref.org/works/10.1016/j.intimp.2025.114556'),
 ('EAU_original_pubmed.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:40163942&format=json'),
 ('glaucoma_supplement.zip','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10531004/supplementaryFiles'),
]

def fetch(task):
    name,url=task; p=D/name
    record={'file':name,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (scientific dataset audit)'})
        with urllib.request.urlopen(req,timeout=90) as r:
            content=r.read(); record['content_type']=r.headers.get('Content-Type'); record['resolved_url']=r.url
        p.write_bytes(content)
        record.update(status='downloaded',bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
        print(name,len(content),flush=True)
    except Exception as e:
        record.update(status='failed',error=str(e)); print(name,str(e),flush=True)
    return record

with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
    records=list(ex.map(fetch,tasks))
(O/'download_manifest.json').write_text(json.dumps(records,indent=2))

f=json.loads((D/'OSD-758_files.json').read_text())
def walk(obj):
    if isinstance(obj,dict):
        if 'file_name' in obj or 'fileName' in obj: yield obj
        for v in obj.values(): yield from walk(v)
    elif isinstance(obj,list):
        for v in obj: yield from walk(v)
files=list(walk(f))
wanted=['GLDS-664_rna-seq_RSEM_Unnormalized_Counts.csv','GLDS-664_rna-seq_SampleTable.csv',
 'GLDS-664_rna-seq_MHU-8_RTN_metadata.csv','GLDS-664_rna-seq_contrasts.csv',
 'GLDS-664_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv',
 'GLDS-664_rna_seq_SampleTable_GLbulkRNAseq.csv','GLDS-664_rna_seq_contrasts_GLbulkRNAseq.csv',
 'GLDS-664_rna_seq_bulkRNASeq_v2_runsheet.csv','GLDS-664_rna_seq_processing_info_GLbulkRNAseq.zip',
 'GLDS-664_rna_seq_raw_multiqc_GLbulkRNAseq_data.zip','GLDS-664_rna_seq_trimmed_multiqc_GLbulkRNAseq_data.zip',
 'GLDS-664_rna_seq_align_multiqc_GLbulkRNAseq_data.zip','GLDS-664_rna_seq_read_dist_multiqc_GLbulkRNAseq_data.zip']
extra=[]
for obj in files:
    name=obj.get('file_name',obj.get('fileName'))
    if name in wanted:
        url=obj.get('remote_url',obj.get('remoteUrl'))
        if url and url.startswith('/'): url='https://osdr.nasa.gov'+url
        extra.append((name,url))
if not extra:
    print('NASA keys',str(f)[:1500],flush=True)
else:
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: records+=list(ex.map(fetch,extra))
(O/'download_manifest.json').write_text(json.dumps(records,indent=2))

# Persist only metadata fields, with platform annotation in a separate TSV.
for p in sorted(D.glob('GSE*_family.soft')):
    meta={'series':{},'samples':{},'platforms':{}}; section=None; current=None
    platform_lines=[]; in_platform=False
    with p.open(errors='replace') as handle:
        for line in handle:
            line=line.rstrip('\n')
            if line.startswith('^SERIES = '): section='series'; current=meta[section]
            elif line.startswith('^SAMPLE = '): section='samples'; current={};meta[section][line.split(' = ',1)[1]]=current
            elif line.startswith('^PLATFORM = '): section='platforms'; current={};meta[section][line.split(' = ',1)[1]]=current
            elif line=='!platform_table_begin': in_platform=True
            elif line=='!platform_table_end': in_platform=False
            elif in_platform: platform_lines.append(line)
            elif line.startswith('!') and ' = ' in line and current is not None:
                key,val=line[1:].split(' = ',1); current.setdefault(key,[]).append(val)
    (O/(p.stem.replace('_family','')+'_metadata.json')).write_text(json.dumps(meta,indent=2))
    if platform_lines: (O/(p.stem.replace('_family','')+'_platform.tsv')).write_text('\n'.join(platform_lines)+'\n')
print('Download and metadata preparation complete',flush=True)
