#!/usr/bin/env python3
"""Own-sequence KO evidence; never propagate a representative's annotation."""
import argparse
from pathlib import Path
from collections import defaultdict
from mainline_common import read_table,write_table,save,sha,unique
def summarize(catalog,genes,kofam):
    catalog,genes,out=map(Path,(catalog,genes,kofam))
    gs=unique(read_table(genes/'gene_table.tsv',['gene_id','viral_sequence_id']),'gene_id')
    ss=unique(read_table(out/'gene_ko_status.tsv',['gene_id','status']),'gene_id')
    if set(gs)!=set(ss):raise ValueError('KO gene scope differs from own viral gene catalogue')
    members=read_table(catalog/'vOTU_members.tsv',['viral_sequence_id','votu_id'])
    mapping=defaultdict(set)
    for r in members:
        if r['votu_id']:mapping[r['viral_sequence_id']].add(r['votu_id'])
    counts=defaultdict(set); links=[]
    for h in read_table(out/'gene_ko_hits.tsv'):
        g=gs[h['gene_id']];sid=g['viral_sequence_id']
        for vid in sorted(mapping[sid]) or ['']:
            links.append(dict(h,viral_sequence_id=sid,votu_id=vid))
            if h['passes_kofamscan_rule']=='true':counts[(vid,sid,h['ko'])].add(h['gene_id'])
    import core_kofam
    write_table(out/'viral_member_ko_evidence.tsv',links,core_kofam.HIT+['viral_sequence_id','votu_id'])
    write_table(out/'viral_member_ko_counts.tsv',[dict(votu_id=v,viral_sequence_id=s,ko=k,distinct_gene_count=len(g),gene_ids=';'.join(sorted(g))) for (v,s,k),g in sorted(counts.items())],['votu_id','viral_sequence_id','ko','distinct_gene_count','gene_ids'])
    save(out/'viral_ko_manifest.json',dict(status='completed',genes_sha256=sha(genes/'gene_table.tsv'),membership_sha256=sha(catalog/'vOTU_members.tsv'),semantics='Own viral sequence KO homology; no representative broadcast; not AMG/activity; lack of assignment is not biological absence',outputs={p.name:sha(p) for p in out.glob('*.tsv')}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('catalog','genes','kofam'):p.add_argument('--'+k,required=True)
    a=p.parse_args();summarize(a.catalog,a.genes,a.kofam)
