#!/usr/bin/env python3
"""Source-flank extraction and conservative locus evidence, using locked BLASTn."""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from catalog_policy import config,bundle,stable,reviewed
from mainline_common import read_table,write_table,save,sha,unique,fasta,write_fasta,run,new_output,safe_id
from core_votu import parse_blast

FLANK=['locus_id','anchor_candidate_id','source_host_genome','contig_id','side','query_id','start0','end0','length','extraction_status','reason','sequence_sha256']
EVIDENCE=['catalog_version','policy_version','locus_id','genome_id','anchor_candidate_id','anchor_votu_id','run_status','automatic_state','callability','carriage','value','reason','target_contig','orientation','interval_start0','interval_end0','observed_candidate_ids','evidence_id','evidence_sha256','review_status','reviewed_by','reviewed_on','review_source','review_source_sha256']


def prepared_map(directories):
    result={}
    for d in directories:
        p=Path(d); rr=read_table(p/'contig_map.tsv',['genome_id','contig_id','length','sequence_sha256'])
        gids={r['genome_id'] for r in rr}
        if len(gids)!=1:raise ValueError('Prepared mapping must have exactly one genome')
        gid=gids.pop();safe_id(gid)
        if gid in result:raise ValueError('Duplicate prepared genome')
        result[gid]=p
    return result


def genome(p,gid):
    seqs=fasta(p/'genome.fna'); rr=unique(read_table(p/'contig_map.tsv'),'contig_id')
    if set(seqs)!=set(rr):raise ValueError('Prepared sequence/map IDs differ')
    import hashlib
    for cid,s in seqs.items():
        if rr[cid]['genome_id']!=gid or int(rr[cid]['length'])!=len(s) or hashlib.sha256(s.encode()).hexdigest()!=rr[cid]['sequence_sha256']:raise ValueError('Prepared contig identity differs')
    return seqs


def extract(loci,master,preps,cfg):
    queries={}; meta=[]; definitions=unique(loci,'locus_id')
    by_genome=defaultdict(list)
    for l in loci:
        safe_id(l['locus_id']);cid=l['anchor_candidate_id']
        if cid not in master:raise ValueError('Unknown anchor definition')
        by_genome[master[cid]['genome_id']].append(l)
    length=cfg['flanks']['length_bp']
    for gid,ll in sorted(by_genome.items()):
        if gid not in preps:raise ValueError('Locus source genome missing')
        seqs=genome(preps[gid],gid)
        for l in ll:
            r=master[l['anchor_candidate_id']];contig=r['contig_id']
            if contig not in seqs:raise ValueError('Candidate contig foreign key differs')
            start,end=int(r['start0']),int(r['end0']);s=seqs[contig]
            if not 0<=start<end<=len(s) or end-start!=int(r['length']):raise ValueError('Candidate coordinate bounds differ')
            if __import__('hashlib').sha256(s[start:end].encode()).hexdigest()!=r['sequence_sha256']:raise ValueError('Candidate/source contig slice identity differs')
            # A1 slices are in original contig orientation; caller strand is separate.
            for side in ('left','right'):
                q='flank_'+stable([l['locus_id'],side])[:24]
                a,b=(max(0,start-length),start) if length and side=='left' else (end,min(len(s),end+length)) if length else (start,start)
                seq=s[a:b]
                reason='parameters_pending' if length is None else 'contig_edge_truncated' if len(seq)!=length else 'non_acgt_flank' if set(seq)-set('ACGT') else ''
                ok=not reason
                row=dict(locus_id=l['locus_id'],anchor_candidate_id=r['candidate_id'],source_host_genome=gid,contig_id=contig,side=side,query_id=q,start0=a,end0=b,length=len(seq),extraction_status='ready' if ok else 'unassessable',reason=reason,sequence_sha256=__import__('hashlib').sha256(seq.encode()).hexdigest())
                meta.append(row)
                if ok:queries[q]=seq
    return queries,meta


def normalize_hits(path,queries,seqs,cfg):
    if set(queries)&set(seqs):raise ValueError('Flank/source namespace collision')
    grouped=parse_blast(path,{**queries,**seqs}); accepted=defaultdict(dict)
    for (q,s),hits in grouped.items():
        if q not in queries or s not in seqs:raise ValueError('Wrong flank alignment direction')
        for h in hits:
            if h['qstart']>h['qend']:raise ValueError('Unsupported reversed query coordinates')
            covered=h['qend']-h['qstart']+1
            if h['pident']<cfg['identity_percent'] or covered*100/h['qlen']<cfg['coverage_percent'] or h['evalue']>cfg['blast_evalue']:continue
            # No coordinate projection through unaligned inner flank ends.
            hit=dict(contig_id=s,start0=min(h['sstart'],h['send'])-1,end0=max(h['sstart'],h['send']),strand='+' if h['send']>=h['sstart'] else '-',qstart=h['qstart'],qend=h['qend'],qlen=h['qlen'],identity_percent=h['pident'],coverage_percent=covered*100/h['qlen'])
            key=(s,hit['start0'],hit['end0'],hit['strand'],h['qstart'],h['qend'])
            accepted[q][key]=hit
    return {q:list(v.values()) for q,v in accepted.items()}


def assess(locus,flanks,hits,seqs,candidates,complete_ids,maxgap):
    base=dict(automatic_state='unassessable',reason='',target_contig='',orientation='',interval_start0='',interval_end0='',observed_candidate_ids='')
    if any(f['extraction_status']!='ready' for f in flanks):return dict(base,reason='source_flank_unassessable')
    left,right=[next(f for f in flanks if f['side']==side) for side in ('left','right')]
    lh,rh=hits.get(left['query_id'],[]),hits.get(right['query_id'],[])
    if not lh or not rh:return dict(base,reason='one_or_both_flanks_not_detected')
    if len(lh)!=1 or len(rh)!=1:return dict(base,reason='multiple_flank_locations')
    l,r=lh[0],rh[0]
    if l['qend']!=l['qlen'] or r['qstart']!=1:return dict(base,reason='inner_flank_end_not_aligned')
    if l['contig_id']!=r['contig_id']:return dict(base,reason='flanks_on_different_contigs')
    if l['strand']!=r['strand']:return dict(base,reason='inconsistent_flank_orientation')
    start,end=(l['end0'],r['start0']) if l['strand']=='+' else (r['end0'],l['start0'])
    base.update(target_contig=l['contig_id'],orientation=l['strand'],interval_start0=start,interval_end0=end)
    if start>end:return dict(base,reason='overlapping_or_reversed_flanks')
    span=seqs[l['contig_id']][min(l['start0'],r['start0']):max(l['end0'],r['end0'])]
    if set(span)-set('ACGT'):return dict(base,reason='non_acgt_or_gap_between_flanks')
    overlap=[c for c in candidates if c['contig_id']==l['contig_id'] and int(c['start0'])<end and int(c['end0'])>start]
    inside=[c for c in overlap if start<=int(c['start0'])<int(c['end0'])<=end and c['candidate_id'] in complete_ids]
    base['observed_candidate_ids']=';'.join(sorted(c['candidate_id'] for c in overlap))
    if inside and len(inside)==len(overlap):return dict(base,automatic_state='occupied_supported',reason='qualified_candidate_inside_paired_flanks')
    if overlap:return dict(base,reason='unqualified_or_boundary_spanning_candidate_in_interval')
    if end-start<=maxgap:return dict(base,automatic_state='empty_site_supported',reason='unique_paired_continuous_empty_site')
    return dict(base,reason='intervening_sequence_not_explained')


def execute(a):
    out=new_output(a.outdir); cfg=config(a.policy);ev,bh=bundle(a.evidence)
    catalog=Path(a.catalog);manifest=json.loads((catalog/'catalog_manifest.json').read_text())
    if manifest.get('policy')!=cfg or manifest.get('evidence_bundle_sha256')!=bh:raise ValueError('Catalog policy/bundle version mismatch')
    master=unique(read_table(catalog/'prophage_master.tsv'),'candidate_id');members=unique(read_table(catalog/'vOTU_members.tsv'),'candidate_id')
    preps=prepared_map(a.prepared_dirs);loci=ev['locus_definitions.tsv']
    if not {r['genome_id'] for r in master.values()}<=set(preps):raise ValueError('Candidate source genome missing')
    for l in loci:
        m=members.get(l['anchor_candidate_id'],{})
        if not m.get('votu_id') or m.get('complete_element_eligible')!='true':raise ValueError('Locus must reference a qualified, unambiguous anchor')
    queries,flanks=extract(loci,master,preps,cfg)
    write_fasta(out/'flanks.fna',queries);write_table(out/'flank_sources.tsv',flanks,FLANK)
    complete={c for c,m in members.items() if m['complete_element_eligible']=='true'}
    reviews={(r['locus_id'],r['genome_id']):r for r in ev['locus_reviews.tsv']}
    if not set(reviews)<={(l['locus_id'],g) for l in loci for g in preps}:raise ValueError('Locus review foreign key differs')
    results=[];executions=[];all_hits=[]
    ready=all(cfg['flanks'].get(k) is not None for k in ('length_bp','identity_percent','coverage_percent','empty_site_max_gap_bp','blast_evalue'))
    for gid,p in sorted(preps.items()):
        seqs=genome(p,gid);dd=out/'raw'/gid;dd.mkdir(parents=True)
        hits={};status='not_assessed';err='parameters_pending' if not ready else 'no_usable_flanks';raw_hash=''
        if ready and queries:
            try:
                target=dd/'alignments.tsv'
                # Subject mode retains all target contigs; one genome at a time.
                run([a.blastn,'-task','blastn','-query',out/'flanks.fna','-subject',p/'genome.fna','-evalue',str(cfg['flanks']['blast_evalue']),'-max_target_seqs',str(max(1,len(seqs))),'-num_threads',str(a.cpus),'-outfmt','6 std qlen slen qseq sseq','-out',target],dd,'blastn',timeout=a.timeout)
                hits=normalize_hits(target,queries,seqs,cfg['flanks']);raw_hash=sha(target);status='success';err=''
            except (RuntimeError,OSError) as e:status='failed';err=str(e)
        executions.append(dict(genome_id=gid,run_status=status,reason=err,genome_fasta_sha256=sha(p/'genome.fna'),alignment_sha256=raw_hash))
        for q,hs in hits.items():
            for h in hs:all_hits.append(dict(genome_id=gid,query_id=q,**h))
        for l in loci:
            cid=l['anchor_candidate_id'];lf=[f for f in flanks if f['locus_id']==l['locus_id']]
            rr=assess(l,lf,hits,seqs,[r for r in master.values() if r['genome_id']==gid],complete,cfg['flanks']['empty_site_max_gap_bp']) if status=='success' else dict(automatic_state='unassessable',reason=err,target_contig='',orientation='',interval_start0='',interval_end0='',observed_candidate_ids='')
            # Review is bound to the generated evidence hash, not merely locus/genome names.
            payload=dict(catalog_version=manifest['catalog_version'],policy=cfg,bundle_input_tables={k:v for k,v in ev.items() if k!='locus_reviews.tsv'},locus=l,genome_id=gid,source_sha256=executions[-1]['genome_fasta_sha256'],alignment_sha256=raw_hash,result=rr)
            eh=stable(payload);eid='locus_ev_'+eh[:24];review=reviews.get((l['locus_id'],gid),{})
            reviewed_ok=reviewed(review) and review['evidence_id']==eid and review['decision']==rr['automatic_state']
            # Synthetic demonstration can generate binary states; proposed production rules cannot.
            rule_ready=cfg['status'] in ('synthetic','frozen')
            call='assessable' if rule_ready and rr['automatic_state']=='occupied_supported' else 'unassessable'
            carriage='present' if call=='assessable' else 'unknown';value=1 if carriage=='present' else 'NA'
            if rule_ready and rr['automatic_state']=='empty_site_supported' and reviewed_ok:call='assessable';carriage='absent_assessable';value=0
            results.append(dict(catalog_version=manifest['catalog_version'],policy_version=cfg['policy_version'],locus_id=l['locus_id'],genome_id=gid,anchor_candidate_id=cid,anchor_votu_id=members[cid]['votu_id'],run_status=status,**rr,callability=call,carriage=carriage,value=value,evidence_id=eid,evidence_sha256=eh,review_status='reviewed' if reviewed_ok else 'stale_or_conflicting_review' if review else 'not_reviewed',reviewed_by=review.get('reviewed_by',''),reviewed_on=review.get('reviewed_on',''),review_source=review.get('evidence_source',''),review_source_sha256=review.get('evidence_sha256','')))
    write_table(out/'locus_evidence.tsv',results,EVIDENCE)
    write_table(out/'locus_execution.tsv',executions,['genome_id','run_status','reason','genome_fasta_sha256','alignment_sha256'])
    write_table(out/'flank_hits.tsv',all_hits,['genome_id','query_id','contig_id','start0','end0','strand','qstart','qend','qlen','identity_percent','coverage_percent'])
    write_table(out/'locus_review_requests.tsv',[dict(locus_id=r['locus_id'],genome_id=r['genome_id'],evidence_id=r['evidence_id'],decision=r['automatic_state'],review_status='pending',reviewed_by='',reviewed_on='',evidence_source='',evidence_sha256='') for r in results if r['automatic_state']=='empty_site_supported'],['locus_id','genome_id','evidence_id','decision','review_status','reviewed_by','reviewed_on','evidence_source','evidence_sha256'])
    save(out/'locus_manifest.json',dict(status='completed_with_failures' if any(r['run_status']=='failed' for r in executions) else 'completed',catalog_version=manifest['catalog_version'],policy_sha256=sha(a.policy),bundle_sha256=bh,scope='locus occupancy; never whole-genome vOTU absence',coordinates='0-based half-open',review_truthfulness='declared reviewer provenance; not independently verified',scientific_calibration='not_performed',outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))


def apply_reviews(a):
    """Pure metadata review application, avoiding another BLAST invocation."""
    import shutil
    out=new_output(a.outdir);cfg=config(a.policy);ev,bh=bundle(a.evidence)
    src=Path(a.review_existing);old=json.loads((src/'locus_manifest.json').read_text())
    cat=json.loads((Path(a.catalog)/'catalog_manifest.json').read_text())
    if old['policy_sha256']!=sha(a.policy) or old['catalog_version']!=cat['catalog_version']:raise ValueError('Review cannot change policy/catalog identity')
    for name,h in old['outputs'].items():
        p=(src/name).resolve()
        if p.parent!=src.resolve() or sha(p)!=h:raise ValueError('Changed source locus evidence')
        shutil.copyfile(p,out/name)
    rows=read_table(src/'locus_evidence.tsv');reviews={(r['locus_id'],r['genome_id']):r for r in ev['locus_reviews.tsv']}
    if not set(reviews)<={(r['locus_id'],r['genome_id']) for r in rows}:raise ValueError('Review refers to unknown evidence row')
    for r in rows:
        v=reviews.get((r['locus_id'],r['genome_id']),{})
        ok=reviewed(v) and v['evidence_id']==r['evidence_id'] and v['decision']==r['automatic_state']
        if r['automatic_state']=='empty_site_supported':
            allow=ok and cfg['status'] in ('frozen','synthetic') and r['run_status']=='success'
            r.update(callability='assessable' if allow else 'unassessable',carriage='absent_assessable' if allow else 'unknown',value='0' if allow else 'NA')
        r.update(review_status='reviewed' if ok else 'stale_or_conflicting_review' if v else 'not_reviewed',reviewed_by=v.get('reviewed_by',''),reviewed_on=v.get('reviewed_on',''),review_source=v.get('evidence_source',''),review_source_sha256=v.get('evidence_sha256',''))
    write_table(out/'locus_evidence.tsv',rows,EVIDENCE)
    old.update(review_application='metadata_only_no_tool_execution',source_locus_manifest=str((src/'locus_manifest.json').resolve()),source_locus_manifest_sha256=sha(src/'locus_manifest.json'),bundle_sha256=bh,outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    save(out/'locus_manifest.json',old)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('catalog','policy','evidence','outdir'):p.add_argument('--'+n,required=True)
    p.add_argument('--prepared-dirs',nargs='+');p.add_argument('--review-existing');p.add_argument('--blastn',default='blastn');p.add_argument('--cpus',type=int,default=1,choices=range(1,7));p.add_argument('--timeout',type=int)
    a=p.parse_args()
    if a.review_existing:apply_reviews(a)
    elif a.prepared_dirs:execute(a)
    else:p.error('Supply prepared directories or existing evidence for metadata-only review')


if __name__=='__main__':
    try:main()
    except Exception as e:print('ERROR: '+str(e),file=sys.stderr);raise SystemExit(1)
