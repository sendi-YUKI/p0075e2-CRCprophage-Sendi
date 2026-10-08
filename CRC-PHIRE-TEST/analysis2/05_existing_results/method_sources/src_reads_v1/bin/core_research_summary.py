#!/usr/bin/env python3
"""Join member functions, locus evidence and AimBeta metadata without inference."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from mainline_common import read_table,write_table,save,sha,unique,new_output
from catalog_policy import stable


def summarize(a):
    out=new_output(a.outdir);cat=Path(a.catalog);loc=Path(a.locus)
    cm=json.loads((cat/'catalog_manifest.json').read_text());lm=json.loads((loc/'locus_manifest.json').read_text())
    if cm['catalog_version']!=lm['catalog_version']:raise ValueError('Mixed catalog versions')
    if json.loads((Path(a.function_summary)/'function_summary_status.json').read_text())['input_votu_sha256']!=sha(cat/'vOTU_members.tsv'):raise ValueError('Stale member function summary')
    version=cm['catalog_version'];master=unique(read_table(cat/'prophage_master.tsv'),'candidate_id');members=unique(read_table(cat/'vOTU_members.tsv'),'candidate_id')
    genomes=unique(read_table(Path(a.qc)/'genome_qc.tsv'),'genome_id')
    genes=unique([r for d in a.gene_dirs for r in read_table(Path(d)/'gene_table.tsv')],'gene_id')
    if not {r['genome_id'] for r in genes.values()}<=set(genomes):raise ValueError('Gene source host missing')
    loci=read_table(loc/'locus_evidence.tsv');features={};links=[];observations=[];proof=[];gene_links=[];functions=[]
    def feature(fid,unit,gid='',candidate='',contig='',start=None,end=None,strand='',artifact='',row=''):
        r=dict(feature_id=fid,feature_unit=unit,genome_id=gid or None,candidate_id=candidate or None,contig_id=contig or None,start=start,end=end,strand=strand or None,coordinate_system='zero_based_half_open',source_artifact=artifact or None,source_row_id=row or None)
        if fid in features and features[fid]!=r:raise ValueError('Feature identity collision')
        features[fid]=r
    def link(parent,member,rel,source):
        r=dict(link_id='link_'+stable([parent,member,rel])[:24],parent_feature_id=parent,member_feature_id=member,relation=rel,evidence_ref=source)
        if r not in links:links.append(r)
    for cid,r in master.items():
        feature(cid,'candidate',r['genome_id'],contig=r['contig_id'],start=int(r['start0']),end=int(r['end0']),strand=r['strand'],artifact='prophage_master.tsv',row=cid)
    for gid,r in genes.items():
        feature(gid,'gene',r['genome_id'],contig=r['contig_id'],start=int(r['start0']),end=int(r['end0']),strand=r['strand'],artifact='gene_table.tsv',row=gid)
    for m in members.values():
        if m['votu_id']:
            feature(m['votu_id'],'votu',artifact='vOTU_members.tsv',row=m['votu_id'])
            link(m['votu_id'],m['candidate_id'],'votu_member','vOTU_members.tsv:'+m['candidate_id'])
    def observation(gid,fid,unit,state,value,call,run='success',rule='',eid=''):
        oid='obs_'+stable([version,gid,fid,'carriage'])[:24]
        observations.append(dict(observation_id=oid,entity_type='genome',entity_id=gid,feature_id=fid,feature_unit=unit,run_status=run,carriage=state,measurement='unassessable',value=value,callability=call,callability_evidence_id=eid or None,rule_version=rule or None,value_axis='carriage_binary',value_unit='binary'))
        return oid
    locus_candidates=defaultdict(set);locus_links=[]
    for r in loci:
        if r['catalog_version']!=version or r['genome_id'] not in genomes:raise ValueError('Invalid locus provenance')
        fid='locus_'+stable([version,r['locus_id'],r['genome_id']])[:24]
        # Empty sites can be zero length; preserve their interval in evidence, not a fake positive-length feature.
        feature(fid,'locus',r['genome_id'],artifact='locus_evidence.tsv',row=r['evidence_id'])
        for cid in filter(None,r['observed_candidate_ids'].split(';')):
            if cid not in master or master[cid]['genome_id']!=r['genome_id']:raise ValueError('Locus/candidate source host mismatch')
            link(fid,cid,'locus_candidate','locus_evidence.tsv:'+r['evidence_id']);locus_candidates[cid].add(fid)
        evidence_id='proof_'+stable([r['evidence_id'],'callability'])[:24]
        oid=observation(r['genome_id'],fid,'locus',r['carriage'],None if r['value']=='NA' else int(r['value']),r['callability'],r['run_status'],r['policy_version'],evidence_id)
        accepted=r['review_status']=='reviewed' and r['carriage']=='absent_assessable'
        proof.append(dict(evidence_id=evidence_id,entity_type='observation',entity_id=oid,field='callability',value=r['callability'],source_uri='locus_evidence.tsv',locator=r['evidence_id'],level='B' if accepted else 'U',review_state='accepted' if accepted else 'pending',reviewer=r['reviewed_by'] if accepted else None,review_date=r['reviewed_on'] if accepted else None))
        locus_links.append(dict(locus_feature_id=fid,locus_group=r['locus_id'],genome_id=r['genome_id'],anchor_votu_id=r['anchor_votu_id'],evidence_id=r['evidence_id']))
    genome_carriage=[];whole={};whole_manifest=None
    if getattr(a,'whole_genome_evidence',None):
        from genome_callability import load_decisions
        rows,whole_manifest=load_decisions(a.whole_genome_evidence,cat)
        whole={(r['genome_id'],r['votu_id']):r for r in rows}
        expected={(g,v) for g in genomes for v in {m['votu_id'] for m in members.values()}-{''}}
        if set(whole)!=expected:raise ValueError('Whole-genome decision scope differs from research catalogue')
        import shutil
        for name in ('genome_callability.tsv','callability_policy_receipt.json'):
            shutil.copyfile(Path(a.whole_genome_evidence)/name,out/name)
    for v in sorted({m['votu_id'] for m in members.values()}-{''}):
        for gid in sorted(genomes):
            observed=[m for m in members.values() if m['votu_id']==v and m['genome_id']==gid]
            qualified=[m for m in observed if m['complete_element_eligible']=='true']
            state='present' if qualified and cm['policy']['status'] in ('synthetic','frozen') else 'unknown'
            value=1 if state=='present' else None
            decision=whole.get((gid,v))
            if decision:
                if state=='present' and decision['carriage']!='present':raise ValueError('Positive catalogue/callability conflict')
                state=decision['carriage'];value=None if decision['value']=='NA' else int(decision['value'])
            eid=decision['evidence_id'] if decision and state=='absent_assessable' else ''
            oid=observation(gid,v,'votu',state,value,'assessable' if state!='unknown' else 'unassessable',
                rule=decision['rule_version'] if decision else cm['policy']['policy_version'],eid=eid)
            if eid:
                proof.append(dict(evidence_id=eid,entity_type='observation',entity_id=oid,field='callability',value='assessable',source_uri='genome_callability.tsv',locator=eid,
                    level='B',review_state='accepted',reviewer=decision['reviewer'],review_date=decision['review_date']))
            genome_carriage.append(dict(catalog_version=version,votu_id=v,genome_id=gid,carriage=state,value=value if value is not None else 'NA',observed_instance_count=len(observed),qualified_instance_count=len(qualified),instance_ids=';'.join(sorted(m['candidate_id'] for m in observed)),reason=decision['reason'] if decision else 'qualified_complete_member' if value else 'fragment_only_or_no_assessable_complete_member',scope='whole_genome_target; local_empty_site_never_propagated'))
    seen=set()
    for d in a.function_dirs:
        dd=Path(d)
        for r in read_table(dd/'gene_prophage_membership.tsv'):
            cid,gid=r['candidate_id'],r['gene_id'];m=members.get(cid);g=genes.get(gid)
            if not m or not g or (cid,gid) in seen or m['genome_id']!=g['genome_id']:raise ValueError('Duplicate/invalid gene/member/host link')
            seen.add((cid,gid));c=master[cid]
            overlap=max(0,min(int(g['end0']),int(c['end0']))-max(int(g['start0']),int(c['start0']))) if g['contig_id']==c['contig_id'] else 0
            relation='fully_inside' if int(c['start0'])<=int(g['start0'])<int(g['end0'])<=int(c['end0']) else 'boundary_spanning'
            if overlap<=0 or relation!=r['relation'] or ('overlap_bp' in r and int(r['overlap_bp'])!=overlap):raise ValueError('Gene/candidate coordinate relation differs')
            gene_links.append(dict(catalog_version=version,candidate_id=cid,gene_id=gid,votu_id=m['votu_id'],source_host_genome=g['genome_id'],membership_status=m['membership_status'],member_role=m['member_role'],locus_feature_ids=';'.join(sorted(locus_candidates[cid])),relation=relation))
            link(cid,gid,'candidate_gene','gene_prophage_membership.tsv:'+cid+':'+gid)
        for name in ('member_functions.tsv','system_prophage_membership.tsv'):
            for r in read_table(dd/name):
                cid=r['candidate_id']
                if cid not in members:raise ValueError('Function candidate absent from current catalog')
                functions.append(dict(catalog_version=version,candidate_id=cid,votu_id=members[cid]['votu_id'],source_host_genome=members[cid]['genome_id'],membership_status=members[cid]['membership_status'],member_role=members[cid]['member_role'],source_artifact=name,source_row_id=stable(r),evidence_json=json.dumps(r,sort_keys=True,separators=(',',':'))))
    if a.kofam_summary:
        ko=Path(a.kofam_summary);receipt=json.loads((ko/'kofam_summary.json').read_text())
        if receipt['membership_source_sha256']!=sha(cat/'vOTU_members.tsv'):raise ValueError('Stale KO catalog membership')
        for r in read_table(ko/'member_ko_evidence.tsv'):
            cid=r['candidate_id']
            if cid not in members or (cid,r['gene_id']) not in seen or r['votu_id']!=members[cid]['votu_id'] or r['genome_id']!=members[cid]['genome_id']:raise ValueError('KO member/gene/host mismatch')
            functions.append(dict(catalog_version=version,candidate_id=cid,votu_id=r['votu_id'],source_host_genome=r['genome_id'],membership_status=members[cid]['membership_status'],member_role=members[cid]['member_role'],source_artifact='member_ko_evidence.tsv',source_row_id=stable(r),evidence_json=json.dumps(r,sort_keys=True,separators=(',',':'))))
    write_table(out/'research_member_gene_links.tsv',gene_links,['catalog_version','candidate_id','gene_id','votu_id','source_host_genome','membership_status','member_role','locus_feature_ids','relation'])
    write_table(out/'research_function_links.tsv',functions,['catalog_version','candidate_id','votu_id','source_host_genome','membership_status','member_role','source_artifact','source_row_id','evidence_json'])
    write_table(out/'genome_votu_carriage.tsv',genome_carriage,['catalog_version','votu_id','genome_id','carriage','value','observed_instance_count','qualified_instance_count','instance_ids','reason','scope'])
    write_table(out/'locus_feature_crosswalk.tsv',locus_links,['locus_feature_id','locus_group','genome_id','anchor_votu_id','evidence_id'])
    schema=json.loads(Path(a.study_schema).read_text())['$defs']['rows']
    for name,rr in [('features',list(features.values())),('feature_links',links),('observations',observations),('evidence',proof)]:
        write_table(out/('research_'+name+'.tsv'),rr,list(schema[name]['properties']))
    save(out/'research_summary.json',dict(status='completed_with_failures' if lm['status']=='completed_with_failures' else 'completed',catalog_version=version,policy_status=cm['policy']['status'],member_count=len(members),genome_count=len(genomes),gene_links=len(gene_links),function_links=len(functions),locus_rows=len(loci),study_merge='Merge these four feature/observation/evidence tables into an explicitly supplied clinical contract; do not invent participants or source labels',whole_genome_absence_implemented=bool(whole_manifest),operational_nondetection_only=True,whole_genome_source_manifest=whole_manifest,calibration='pending',real_validation='not_performed',outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('catalog','locus','qc','function-summary','study-schema','outdir'):p.add_argument('--'+n,required=True)
    for n in ('gene-dirs','function-dirs'):p.add_argument('--'+n,nargs='+',required=True)
    p.add_argument('--kofam-summary');p.add_argument('--whole-genome-evidence');summarize(p.parse_args())
