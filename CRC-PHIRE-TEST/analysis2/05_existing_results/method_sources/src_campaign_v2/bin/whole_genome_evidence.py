#!/usr/bin/env python3
"""Full assembled-genome target search; evidence only, never automatic absence."""
import argparse,hashlib,json,os,sys
from pathlib import Path
from collections import defaultdict
import sys
sys.path.insert(0, str(Path(os.environ.get("CRC_SPARK_REPO","/srv/CRC-PHIRE/projects/crc-pipeline"))/"bin"))
from mainline_common import fasta,write_fasta,read_table,write_table,sha,save,new_output,run,unique
from core_locus import prepared_map,genome
import spark_runtime
FIELDS=['catalog_version','genome_id','votu_id','target_candidate_id','target_sha256','prepared_genome_sha256','assembled_bp','contig_count','non_acgt_bp','reported_hsps','query_covered_bp_union','query_covered_fraction_union','max_reported_identity','edge_hits','search_status','carriage','callability','reason']
def union_size(intervals):
    end=0;total=0
    for a,b in sorted(intervals):
        total+=max(0,b-max(end,a));end=max(end,b)
    return total
def config(path,execute=False):
    rule=json.loads(Path(path).read_text())
    if rule.get('schema_version')!=1 or rule.get('status') not in ('proposed','synthetic','frozen'):raise ValueError('Unknown evidence configuration')
    if rule.get('absence_rule'):
        from genome_callability import validate_policy
        validate_policy(rule['absence_rule'])
    if not execute:return rule
    for k in ('blast_evalue','word_size','edge_bp','rule_version'):
        if rule.get(k) is None:raise ValueError('Explicit search settings required: '+k)
    if not 1e-50<=float(rule['blast_evalue'])<=1 or not 4<=int(rule['word_size'])<=64 or not 0<=int(rule['edge_bp'])<=10000:raise ValueError('Invalid search settings')
    return rule

def execute(a):
    alloc=spark_runtime.allocation(require_lock=True)
    if not 1<=a.cpus<=min(6,alloc['cpus']):raise ValueError('CPU request exceeds allocation')
    rule=config(a.config,True)
    catalog=Path(a.catalog);manifest=json.loads((catalog/'catalog_manifest.json').read_text())
    reps=read_table(catalog/'vOTU_representatives.tsv',['votu_id','candidate_id'])
    unique(reps,'votu_id');seq=fasta(catalog/'candidates.fna')
    queries={'q'+str(i):seq[r['candidate_id']] for i,r in enumerate(reps)}
    qmeta={'q'+str(i):r for i,r in enumerate(reps)}
    preps=prepared_map(a.prepared_dirs)
    dest=Path(a.outdir).resolve()
    for source in [catalog,Path(a.config),*preps.values()]:
        pp=source.resolve()
        if dest==pp or dest.is_relative_to(pp) or pp.is_relative_to(dest):raise ValueError('Output overlaps source')
    out=new_output(dest);write_fasta(out/'targets.fna',queries)
    records=[];sources=[]
    for gi,(gid,p) in enumerate(sorted(preps.items())):
        ss=genome(p,gid);renamed={'s'+str(i):s for i,s in enumerate(ss.values())}
        st={name:key for name,key in zip(renamed,ss)}
        work=out/('genome_'+str(gi));work.mkdir();write_fasta(work/'subject.fna',renamed)
        raw=work/'hits.tsv';hits=defaultdict(list)
        if queries and renamed:
            run(['blastn','-query',out/'targets.fna','-subject',work/'subject.fna','-task','blastn','-word_size',str(rule['word_size']),
                 '-evalue',str(rule['blast_evalue']),'-dust','no','-num_threads',str(a.cpus),
                 '-max_target_seqs',str(max(5,len(renamed))),'-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
                 '-out',raw],work,'whole_genome_blast',timeout=a.timeout)
        else:raw.write_text('')
        normalized=[]
        for line in raw.read_text().splitlines():
            cells=line.split('\t')
            if len(cells)!=12:raise ValueError('Invalid BLAST columns')
            q,t=cells[:2]
            if q not in queries or t not in renamed:raise ValueError('Foreign BLAST identity')
            ident=float(cells[2]);qa,qb=map(int,cells[4:6]);sa,sb=map(int,cells[6:8])
            if not 1<=min(qa,qb)<=max(qa,qb)<=len(queries[q]) or not 1<=min(sa,sb)<=max(sa,sb)<=len(renamed[t]):raise ValueError('BLAST coordinate bounds')
            h=dict(query_id=q,contig_id=st[t],identity=ident,query_start0=min(qa,qb)-1,query_end0=max(qa,qb),
                subject_start0=min(sa,sb)-1,subject_end0=max(sa,sb),strand='+' if sa<=sb else '-',evalue=cells[8],bitscore=cells[9],
                edge_hit=min(sa,sb)-1<rule['edge_bp'] or len(renamed[t])-max(sa,sb)<rule['edge_bp'])
            hits[q].append(h);normalized.append(h)
        write_table(work/'normalized_hits.tsv',normalized,['query_id','contig_id','identity','query_start0','query_end0','subject_start0','subject_end0','strand','evalue','bitscore','edge_hit'])
        sources.append(dict(genome_id=gid,prepared_genome_sha256=sha(p/'genome.fna'),contig_map_sha256=sha(p/'contig_map.tsv'),
            original_source_note='prepared map and parent manifest retain original genome identity; prepared bytes are not original bytes',
            search_directory=work.name,raw_sha256=sha(raw)))
        for q,target in qmeta.items():
            hh=hits[q];covered=union_size([(h['query_start0'],h['query_end0']) for h in hh])
            records.append(dict(catalog_version=manifest['catalog_version'],genome_id=gid,votu_id=target['votu_id'],target_candidate_id=target['candidate_id'],
                target_sha256=hashlib.sha256(queries[q].encode()).hexdigest(),prepared_genome_sha256=sha(p/'genome.fna'),assembled_bp=sum(map(len,ss.values())),
                contig_count=len(ss),non_acgt_bp=sum(sum(c not in 'ACGT' for c in seq) for seq in ss.values()),reported_hsps=len(hh),
                query_covered_bp_union=covered,query_covered_fraction_union=covered/len(queries[q]) if queries[q] else '',
                max_reported_identity=max((h['identity'] for h in hh),default=''),edge_hits=sum(h['edge_hit'] for h in hh),
                search_status='reported_homology' if hh else 'no_reported_homology',carriage='unknown',callability='unassessable',
                reason='assembled_sequence_search_only; scientific_detection_and_genome_absence_rule_pending'))
    write_table(out/'whole_genome_target_evidence.tsv',records,FIELDS)
    from genome_callability import decide
    decisions=decide(records,catalog,preps,rule,out)
    save(out/'absence_review_input.json',dict(schema_version=1,catalog_version=manifest['catalog_version'],rule=rule,claims=[],
        targets=len(reps),sources=sources,automatic_absence=False,
        review_requirements=['original genome hash binding','target detection scope/version','assembly completeness and gaps','homologous fragment/edge evidence','validated detection and negative rule','reviewer/date and scoped claim'],
        limitations=['BLAST search is not proof of complete phage absence','union across loci is search coverage, not one intact insertion','unassembled sequence is not observable']))
    save(out/'whole_genome_evidence_manifest.json',dict(status='completed',producer='assembled_genome_BLASTn_target_search_v1',
        catalog_manifest_sha256=sha(catalog/'catalog_manifest.json'),config_sha256=sha(a.config),config_status=rule['status'],
        biological_absence_producer_ready=False,operational_callability_consumer=True,
        decision_counts={s:sum(r['carriage']==s for r in decisions) for s in ('present','absent_assessable','unknown')},
        scientific_calibration=False,outputs={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('catalog','config','outdir'):p.add_argument('--'+k,required=True)
    p.add_argument('--prepared-dirs',nargs='+',required=True);p.add_argument('--cpus',type=int,default=1);p.add_argument('--timeout',type=int,default=900)
    a=p.parse_args()
    if not 1<=a.timeout<=3600:raise ValueError('Timeout outside campaign bound')
    execute(a)
if __name__=='__main__':main()
