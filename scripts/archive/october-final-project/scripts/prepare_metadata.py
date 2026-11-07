from pathlib import Path
import json,csv,xml.etree.ElementTree as ET,re
R=Path(__file__).resolve().parents[1]; D=R/'downloads'; O=R/'outputs'
for p in D.glob('GSE*_family.soft'):
    meta={'series':{},'samples':{},'platforms':{}}; current=None; platform=[];inside=False
    with p.open(errors='replace') as h:
        for l in h:
            l=l.rstrip('\n')
            if l.startswith('^SERIES = '): current=meta['series']
            elif l.startswith('^SAMPLE = '): current={};meta['samples'][l.split(' = ',1)[1]]=current
            elif l.startswith('^PLATFORM = '): current={};meta['platforms'][l.split(' = ',1)[1]]=current
            elif l=='!platform_table_begin': inside=True
            elif l=='!platform_table_end': inside=False
            elif inside: platform.append(l)
            elif l.startswith('!') and ' = ' in l and current is not None:
                k,v=l[1:].split(' = ',1);current.setdefault(k,[]).append(v)
    acc=p.name.split('_')[0]
    (O/(acc+'_metadata.json')).write_text(json.dumps(meta,indent=2))
    if platform: (O/(acc+'_platform.tsv')).write_text('\n'.join(platform)+'\n')
    rows=[]
    for gsm,m in meta['samples'].items():
        row={'accession':gsm}
        for k,v in m.items(): row[k.removeprefix('Sample_')]=' | '.join(v)
        rows.append(row)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (O/(acc+'_sample_metadata.csv')).open('w') as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)
for p in D.glob('PMC*.xml'):
    root=ET.parse(p).getroot()
    text='\n\n'.join(' '.join(e.itertext()) for e in root.findall('.//body'))
    (O/(p.stem+'_article_text.txt')).write_text(text)
    links=[]
    for tag in ['supplementary-material','ext-link','media']:
        for e in root.findall('.//'+tag):
            href=e.get('{http://www.w3.org/1999/xlink}href')
            if href: links.append({'tag':tag,'url':href,'description':' '.join(e.itertext())})
    (O/(p.stem+'_supplement_links.json')).write_text(json.dumps(links,indent=2))
print('Metadata and article text prepared')
