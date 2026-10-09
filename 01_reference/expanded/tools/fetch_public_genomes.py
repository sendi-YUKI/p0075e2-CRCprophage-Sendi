"""Resumable NCBI Datasets v2 metadata only. No FASTA download. Python stdlib."""
import argparse, hashlib, json, time, urllib.request, urllib.parse, urllib.error
from pathlib import Path
from datetime import datetime, timezone

BASE=Path(__file__).resolve().parents[1]
API='https://api.ncbi.nlm.nih.gov/datasets/v2'

def get(url):
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'CRC-PHIRE-reference-inventory/2026-10-08'})
            with urllib.request.urlopen(req,timeout=60) as r: return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code not in (429,500,502,503,504): raise
            if attempt==3: raise
        except (TimeoutError,urllib.error.URLError):
            if attempt==3: raise
        time.sleep(2**attempt)

def main():
    cfg=json.loads((BASE/'selection_scope.json').read_text(encoding='utf-8'))
    raw=BASE/'raw';raw.mkdir(parents=True,exist_ok=True)
    receipts=[]
    for tax in cfg['taxa']:
        key=tax['key']; dest=raw/key; dest.mkdir(exist_ok=True)
        params={'page_size':1000,'filters.is_metagenome_derived':'metagenome_derived_exclude'}
        params.update(tax.get('api_filters',{}))
        total=0; page=1; tokens=set(); expected=None; error=''
        try:
            while True:
                url=API+'/genome/taxon/'+urllib.parse.quote(tax['query'],safe='')+'/dataset_report?'+urllib.parse.urlencode(params)
                f=dest/f'page_{page:04d}.json';meta=dest/f'page_{page:04d}.receipt.json'
                if f.exists() and meta.exists() and json.loads(meta.read_text())['url']==url:
                    j=json.loads(f.read_text(encoding='utf-8'))
                else:
                    j=get(url); b=json.dumps(j,ensure_ascii=False,indent=2).encode('utf-8');f.write_bytes(b)
                    meta.write_text(json.dumps({'url':url,'retrieved_utc':datetime.now(timezone.utc).isoformat(),
                        'sha256':hashlib.sha256(b).hexdigest()},indent=2),encoding='utf-8')
                    time.sleep(.4)
                total+=len(j.get('reports',[]));expected=j.get('total_count',expected)
                token=j.get('next_page_token')
                if not token:break
                if token in tokens:raise RuntimeError('pagination token repeated')
                tokens.add(token);params['page_token']=token;page+=1
                if page>200:raise RuntimeError('safety page limit exceeded; scope review required')
            status='complete' if total==expected else 'count_mismatch'
        except Exception as e:
            status='failed';error=str(e)
        entry={'key':key,'query':tax['query'],'records':total,'reported_total':expected,'pages':page,
               'status':status,'error':error,'filters':{k:v for k,v in params.items() if k!='page_token'}}
        receipts.append(entry)
        (BASE/'query_receipts.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps(entry,ensure_ascii=True),flush=True)
    if any(x['status']!='complete' for x in receipts):raise SystemExit(2)

if __name__=='__main__':main()
