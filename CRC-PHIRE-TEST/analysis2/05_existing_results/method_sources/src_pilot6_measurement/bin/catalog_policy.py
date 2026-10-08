#!/usr/bin/env python3
"""Opt-in anchor/fragment policy and bounded metadata contracts; no tool calls."""
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from mainline_common import read_table, write_table, save, sha, unique

LIMIT = 1024 * 1024
TABLES = {
    'candidate_reviews.tsv': ['candidate_id','sequence_sha256','boundary_version','viral_identity_support','boundary_support','integration_support','truncated','review_status','reviewed_by','reviewed_on','evidence_id','evidence_source','evidence_sha256','supplemental_approved'],
    'locus_definitions.tsv': ['locus_id','anchor_candidate_id'],
    'locus_reviews.tsv': ['locus_id','genome_id','evidence_id','decision','review_status','reviewed_by','reviewed_on','evidence_source','evidence_sha256'],
    'genome_metadata.tsv': ['genome_id','species','lineage','study_id'],
}
MEMBER = ['catalog_version','votu_id','candidate_id','genome_id','representative_id','membership_status','matching_representative_ids','catalog_tier','member_role','complete_element_eligible','eligibility_reason']


def stable(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def small_json(path):
    p=Path(path)
    if not p.is_file() or p.stat().st_size>LIMIT: raise ValueError('Missing/oversized metadata: '+str(p))
    def pairs(items):
        d={}
        for k,v in items:
            if k in d: raise ValueError('Duplicate JSON key: '+k)
            d[k]=v
        return d
    return json.loads(p.read_text(),object_pairs_hook=pairs,parse_constant=lambda v: (_ for _ in ()).throw(ValueError('Nonfinite number')))


def config(path, frozen=False):
    c=small_json(path)
    if c.get('schema_version')!=1 or c.get('status') not in ('proposed','synthetic','frozen'): raise ValueError('Invalid catalog configuration')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',c.get('policy_version','')): raise ValueError('Missing policy version')
    if c.get('ani_percent')!=95 or c.get('af_shorter_percent')!=85: raise ValueError('This policy retains ANI95/AF85')
    if c.get('require_review_for_absence') is not True: raise ValueError('Reviewed absence evidence is required')
    f=c.get('flanks',{})
    for k in ('length_bp','identity_percent','coverage_percent','empty_site_max_gap_bp','blast_evalue'):
        x=f.get(k)
        if x is None:
            if c['status']!='proposed': raise ValueError('Missing flank setting: '+k)
            continue
        if type(x) not in (int,float) or not math.isfinite(x): raise ValueError('Invalid flank setting: '+k)
        if k in ('length_bp','empty_site_max_gap_bp') and (type(x)!=int or x<(1 if k=='length_bp' else 0)): raise ValueError('Invalid flank integer')
        if k.endswith('percent') and not 0<x<=100: raise ValueError('Flank percent must be 0..100')
        if k=='blast_evalue' and not 0<x<=1: raise ValueError('Invalid E-value')
    if frozen and c['status']!='frozen': raise ValueError('Production catalog policy must be explicitly frozen')
    if c['status']=='frozen' and (not c.get('approved_by') or not c.get('approved_on') or not c.get('calibration_reference')):
        raise ValueError('Frozen policy requires approval and calibration reference')
    return c


def bundle(path):
    root=Path(path).resolve(); m=small_json(root/'manifest.json')
    if m.get('schema_version')!=1 or m.get('purpose')!='catalog_locus_evidence': raise ValueError('Unknown evidence bundle')
    declared={r['path']:r for r in m.get('files',[])}
    if len(declared)!=len(m.get('files',[])) or set(declared)!=set(TABLES): raise ValueError('Evidence bundle must declare four unique tables')
    rows={}
    for name,cols in TABLES.items():
        p=(root/name).resolve()
        if not p.is_relative_to(root) or not p.is_file() or p.stat().st_size>LIMIT or sha(p)!=declared[name]['sha256']: raise ValueError('Changed/unsafe evidence table: '+name)
        rows[name]=read_table(p,cols)
    for name,key in [('candidate_reviews.tsv','candidate_id'),('locus_definitions.tsv','locus_id'),('genome_metadata.tsv','genome_id')]: unique(rows[name],key)
    pairs=[(r['locus_id'],r['genome_id']) for r in rows['locus_reviews.tsv']]
    if len(set(pairs))!=len(pairs): raise ValueError('Duplicate locus review')
    for r in rows['candidate_reviews.tsv']:
        for k in ('truncated','supplemental_approved'):
            if r[k] not in ('true','false'): raise ValueError('Expected explicit boolean: '+k)
        if r['review_status'] not in ('reviewed','pending','rejected'): raise ValueError('Invalid review status')
    return rows,sha(root/'manifest.json')


def reviewed(r):
    return (r.get('review_status')=='reviewed' and bool(r.get('reviewed_by')) and
            bool(re.fullmatch(r'\d{4}-\d{2}-\d{2}',r.get('reviewed_on',''))) and
            bool(r.get('evidence_id')) and bool(r.get('evidence_source')) and
            bool(re.fullmatch('[a-f0-9]{64}',r.get('evidence_sha256',''))))


def build(rows,seqs,pairs,directions,qc,cfg,evidence,bundle_hash,out,old_members,master_fields=None):
    from core_votu import cluster, rank_candidate
    out=Path(out); by_id=unique(rows,'candidate_id'); reviews=unique(evidence['candidate_reviews.tsv'],'candidate_id')
    if not set(reviews)<=set(by_id): raise ValueError('Review refers to unknown candidate')
    meta=unique(evidence['genome_metadata.tsv'],'genome_id')
    if not set(meta)<=set(qc): raise ValueError('Metadata refers to unknown genome')
    primary=[]; supplemental=[]; eligible=[]; assessments=[]
    for r in rows:
        cid=r['candidate_id']; v=reviews.get(cid,{})
        if v and (v['sequence_sha256']!=r['sequence_sha256'] or v['boundary_version']!=r['boundary_version']): raise ValueError('Stale review identity: '+cid)
        host=str(qc.get(r['genome_id'],{}).get('eligibility','unknown')).lower() in ('true','yes','pass','eligible','1')
        complete=reviewed(v) and v['viral_identity_support']=='strong' and v['boundary_support'] in ('dual_sided_supported','experimentally_validated') and v['truncated']=='false'
        tier='unresolved'; reason='identity_or_boundary_not_reviewed'
        if not host: tier='excluded';reason='host_qc_ineligible_or_unresolved'
        elif r.get('element_class') and r['element_class']!='prophage':tier='excluded';reason='element_class_'+r['element_class']
        else:
            eligible.append(r)
            if complete and r.get('checkv_quality') in ('Complete','High-quality'):
                primary.append(r);tier='primary';reason='reviewed_identity_boundary_and_complete_high'
            elif complete and r.get('checkv_quality') in ('Medium-quality','Not-determined') and v['supplemental_approved']=='true':
                supplemental.append(r);tier='supplemental';reason='reviewed_supplemental_approval'
            else: tier='fragment';reason='not_qualified_as_anchor'
        assessments.append(dict(candidate_id=cid,genome_id=r['genome_id'],anchor_tier=tier,reason=reason,host_qc_eligible=str(host).lower(),complete_element_eligible=str(bool(complete and tier in ('primary','supplemental'))).lower(),
            checkv_quality=r.get('checkv_quality',''),viral_identity_support=v.get('viral_identity_support','unknown'),boundary_support=v.get('boundary_support','unknown'),integration_support=v.get('integration_support','unknown'),
            review_status=v.get('review_status','not_assessed'),reviewed_by=v.get('reviewed_by',''),reviewed_on=v.get('reviewed_on',''),evidence_id=v.get('evidence_id',''),evidence_source=v.get('evidence_source',''),evidence_sha256=v.get('evidence_sha256',''),
            raw_viral_identity=r.get('viral_identity_confidence',''),raw_boundary=r.get('boundary_confidence',''),raw_integration=r.get('integration_evidence',''),boundary_version=r.get('boundary_version',''),sequence_sha256=r['sequence_sha256']))
    assessment={r['candidate_id']:r for r in assessments}
    # Cluster only reviewed anchors. No fragment can define or merge a cluster.
    ranked=[]
    for r in primary:
        ranked.append(dict(r,boundary_confidence='experimentally_validated' if reviews[r['candidate_id']]['boundary_support']=='experimentally_validated' else 'independently_supported'))
    reps,_=cluster(ranked,pairs)
    adjacency=defaultdict(set)
    for p in pairs:
        if p['passes']=='true':adjacency[p['candidate_a']].add(p['candidate_b']);adjacency[p['candidate_b']].add(p['candidate_a'])
    # Supplemental anchors that already match primary remain associated, not a second reference.
    remaining=[r for r in supplemental if not any(p in adjacency[r['candidate_id']] for p in reps)]
    supp_reps,_=cluster(remaining,pairs)
    all_reps=reps+supp_reps
    version='research_'+stable(dict(policy=cfg,evidence={k:v for k,v in evidence.items() if k!='locus_reviews.tsv'},rows=rows,qc=sorted(qc.values(),key=lambda x:x['genome_id']),pairs=pairs,primary=reps,supplemental=supp_reps))[:20]
    ids={r:('votu_p_' if r in reps else 'votu_s_')+stable([version,r])[:20] for r in all_reps}
    members=[]; matches_long=[]
    for r in rows:
        cid=r['candidate_id']; a=assessment[cid]
        matches=([cid] if cid in all_reps else [p for p in all_reps if p in adjacency[cid]]) if a['host_qc_eligible']=='true' and not a['reason'].startswith('element_class_') else []
        single=len(matches)==1
        status='representative' if cid in all_reps else 'member' if single else 'ambiguous_multiple_representatives' if matches else 'excluded_element_class' if a['reason'].startswith('element_class_') else 'excluded_host_qc' if a['host_qc_eligible']=='false' else 'unresolved_no_anchor'
        m=dict(catalog_version=version,votu_id=ids[matches[0]] if single else '',candidate_id=cid,genome_id=r['genome_id'],representative_id=matches[0] if single else '',membership_status=status,matching_representative_ids=';'.join(sorted(matches)),catalog_tier=('primary' if matches[0] in reps else 'supplemental') if single else 'ambiguous' if matches else 'unresolved',member_role='anchor' if a['anchor_tier'] in ('primary','supplemental') else 'fragment',complete_element_eligible=a['complete_element_eligible'],eligibility_reason=a['reason'])
        members.append(m)
        for p in matches: matches_long.append(dict(candidate_id=cid,representative_id=p,votu_id=ids[p],catalog_tier='primary' if p in reps else 'supplemental',assignment_status=status))
    members.sort(key=lambda x:x['candidate_id'])
    for a in assessments:a['catalog_version']=version
    write_table(out/'vOTU_members.tsv',members,MEMBER)
    write_table(out/'vOTU_representatives.tsv',[m for m in members if m['membership_status']=='representative'],MEMBER)
    write_table(out/'fragment_assignments.tsv',[m for m in members if m['member_role']=='fragment'],MEMBER)
    write_table(out/'anchor_assignments_all.tsv',matches_long,['candidate_id','representative_id','votu_id','catalog_tier','assignment_status'])
    af=['catalog_version','candidate_id','genome_id','anchor_tier','reason','host_qc_eligible','complete_element_eligible','checkv_quality','viral_identity_support','boundary_support','integration_support','review_status','reviewed_by','reviewed_on','evidence_id','evidence_source','evidence_sha256','raw_viral_identity','raw_boundary','raw_integration','boundary_version','sequence_sha256']
    write_table(out/'catalog_eligibility.tsv',assessments,af)
    for tier,rr in [('primary',reps),('supplemental',supp_reps),('all',all_reps)]:
        name='vOTU_representatives.fna' if tier=='all' else tier+'_representatives.fna'
        (out/name).write_text(''.join('>'+p+'\n'+seqs[p]+'\n' for p in rr))
    updated=[dict(r,main_catalog_eligible=str(assessment[r['candidate_id']]['anchor_tier']=='primary').lower(),main_catalog_reason=assessment[r['candidate_id']]['reason'],catalog_version=version) for r in rows]
    fields=list(rows[0]) if rows else (master_fields or ['candidate_id','genome_id','contig_id','start0','end0','length','boundary_version','sequence_sha256'])
    write_table(out/'prophage_master.tsv',updated,list(dict.fromkeys(fields+['main_catalog_eligible','main_catalog_reason','catalog_version'])))
    old={m['candidate_id']:m for m in old_members}
    write_table(out/'catalog_crosswalk.tsv',[dict(candidate_id=m['candidate_id'],old_votu_id=old.get(m['candidate_id'],{}).get('votu_id',''),new_votu_id=m['votu_id'],catalog_version=version,old_boundary_version=by_id[m['candidate_id']]['boundary_version'],new_boundary_version=by_id[m['candidate_id']]['boundary_version'],relationship='same_sequence_new_explicit_policy') for m in members],['candidate_id','old_votu_id','new_votu_id','catalog_version','old_boundary_version','new_boundary_version','relationship'])
    write_table(out/'catalog_filter_reasons.tsv',[dict(a,element_class=by_id[a['candidate_id']].get('element_class',''),eligible=str(a['host_qc_eligible']=='true' and not a['reason'].startswith('element_class_')).lower()) for a in assessments],['catalog_version','candidate_id','genome_id','element_class','eligible','reason'])
    write_table(out/'votu_pairwise.tsv',pairs,['candidate_a','candidate_b','ani_percent','af_shorter_percent','ani_threshold_percent','af_threshold_percent','passes','status'])
    write_table(out/'votu_directional_ani.tsv',directions,['query_id','subject_id','ani_percent','af_shorter_percent','aligned_columns','identical_columns','query_covered_bp','subject_covered_bp','accepted_hsps'])
    counts=Counter()
    for a in assessments:
        md=meta.get(a['genome_id'],{})
        for axis in ('species','lineage','study_id','checkv_quality'):
            value=a['checkv_quality'] if axis=='checkv_quality' else md.get(axis) or 'unknown'
            counts[(axis,value,a['anchor_tier'])]+=1
    write_table(out/'catalog_denominators.tsv',[dict(axis=k[0],stratum=k[1],eligibility_tier=k[2],instance_count=n,denominator=sum(v for q,v in counts.items() if q[:2]==k[:2])) for k,n in sorted(counts.items())],['axis','stratum','eligibility_tier','instance_count','denominator'])
    for name,rr in [('primary',reps),('supplemental',supp_reps)]:
        write_table(out/(name+'_representatives.tsv'),[m for m in members if m['candidate_id'] in rr],MEMBER)
    write_table(out/'cross_host_clusters.tsv',[dict(votu_id=v,source_host_genomes=';'.join(sorted({m['genome_id'] for m in members if m['votu_id']==v})),source_species=';'.join(sorted({meta.get(m['genome_id'],{}).get('species') or 'unknown' for m in members if m['votu_id']==v}))) for v in sorted(ids.values())],['votu_id','source_host_genomes','source_species'])
    from mainline_common import write_fasta
    write_fasta(out/'candidates.fna',seqs)
    save(out/'catalog_manifest.json',dict(status='completed',catalog_version=version,policy=cfg,evidence_bundle_sha256=bundle_hash,primary_representatives=len(reps),supplemental_representatives=len(supp_reps),candidates=len(rows),fragment_bridging=False,scientific_status='calibration_pending' if cfg['status']!='frozen' else 'configured_not_runtime_validated',outputs={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}))


def inspect_policy(args,require_frozen=False):
    p=getattr(args,'catalog_policy',None)
    if not p:return {'status':'not_enabled'}
    b=getattr(args,'catalog_evidence',None)
    if not b:raise ValueError('Catalog policy requires an evidence bundle')
    for x in (p,b):
        if any(c in str(x) for c in "'\n\r\0"):raise ValueError('Unsafe catalog path')
    c=config(p,require_frozen); rows,h=bundle(b)
    return dict(status='configuration_'+c['status'],policy_sha256=sha(p),bundle_sha256=h,table_rows={k:len(v) for k,v in rows.items()},analysis_started=False,calibration='pending' if c['status']!='frozen' else 'declared_reference_only')
