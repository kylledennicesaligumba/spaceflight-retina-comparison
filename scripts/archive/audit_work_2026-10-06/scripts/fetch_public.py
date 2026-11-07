import urllib.request,concurrent.futures,pathlib,json
p=pathlib.Path(__file__).resolve().parents[1]/'downloads'
jobs={}
for acc in ['GSE241782','GSE156533','GSE284505','GSE314043']:
 d=json.loads((p.parent/'outputs'/(acc+'_family_metadata.json')).read_text());files=d[acc].get('Series_supplementary_file',[])
 for u in files:
  if any(x in u for x in ['gene_count_matrix','readcounts','All_genes_FPKMS','HTSeq.counts','genecount']):jobs[u.split('/')[-1]]=u.replace('ftp://','https://')
jobs['GSE87433_series_matrix.txt.gz']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87433/matrix/GSE87433_series_matrix.txt.gz'
jobs['OSD758_files_all.json']='https://osdr.nasa.gov/osdr/data/osd/files/758'
jobs['OSD758_biodata.json']='https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-758/assay/*/sample/*/'
for pmid in ['39596110','37762022','32960171','28321468']:
 jobs[pmid+'_epmc.json']='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:'+pmid+'&format=json'
def get(kv):
 k,u=kv
 try:
  b=urllib.request.urlopen(u,timeout=60).read();(p/k).write_bytes(b);return k,len(b)
 except Exception as e:return k,str(e)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
 for x in ex.map(get,jobs.items()):print(x,flush=True)
