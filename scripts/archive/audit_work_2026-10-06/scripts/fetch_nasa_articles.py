import urllib.request,concurrent.futures,pathlib,json,xml.etree.ElementTree as ET
p=pathlib.Path(__file__).resolve().parents[1];dl=p/'downloads';jobs={}
d=json.loads((dl/'OSD758_files_all.json').read_text())['studies']['OSD-758']['study_files']
for x in d:
 n=x['file_name']
 if n in ['GLDS-664_rna-seq_RSEM_Unnormalized_Counts.csv','GLDS-664_rna-seq_SampleTable.csv','GLDS-664_rna-seq_MHU-8_RTN_metadata.csv','GLDS-664_rna-seq_contrasts.csv','GLDS-664_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv','GLDS-664_rna_seq_SampleTable_GLbulkRNAseq.csv','GLDS-664_rna_seq_contrasts_GLbulkRNAseq.csv','GLDS-664_rna_seq_processing_info_GLbulkRNAseq.zip','GLDS-664_rna_seq_bulkRNASeq_v2_runsheet.csv']:
  jobs[n]='https://osdr.nasa.gov'+x['remote_url']
for f in dl.glob('PMC*.xml'):
 root=ET.parse(f).getroot();body=root.find('body');txt=[]
 for el in body.iter():
  if el.tag in ['title','p','caption']:txt.append(' '.join(el.itertext()))
 (p/'outputs'/(f.stem+'_article.txt')).write_text('\n\n'.join(txt))
 print('\nARTICLE',f.stem)
 for s in body.findall('.//sec'):
  ti=s.find('title');title=''.join(ti.itertext()) if ti is not None else ''
  if any(k in title.lower() for k in ['animal','sequencing','transcriptom','sample','ocular hypertension','optic nerve crush','rna-seq']):
   print(title,' '.join(s.itertext())[:11000])
 for x in root.findall('.//supplementary-material'):print('SUPP',' '.join(x.itertext()),[(a.tag,a.attrib) for a in x.iter() if a.attrib])
def get(kv):
 k,u=kv
 try:
  b=urllib.request.urlopen(u,timeout=90).read();(dl/k).write_bytes(b);return k,len(b)
 except Exception as e:return k,str(e)
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
 for x in ex.map(get,jobs.items()):print('DOWNLOAD',x,flush=True)
