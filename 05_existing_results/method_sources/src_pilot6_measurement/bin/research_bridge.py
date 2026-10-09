#!/usr/bin/env python3
"""Export a research catalogue for existing V1/B1 consumers without reclustering.
Independent V1 retains its own catalogue and gene set; cross-catalogue links are
evidence, never automatic transfer of source host, functional genes or vOTU IDs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
from collections import defaultdict
import sys
sys.path.insert(0, str(Path(os.environ.get("CRC_SPARK_REPO","/srv/CRC-PHIRE/projects/crc-pipeline"))/"bin"))
from mainline_common import fasta, new_output, read_table, write_table, save, sha, unique, run
from viral_catalog import load_units, make_occurrence, OCC_FIELDS, SEQ_FIELDS, MEMBER_FIELDS
import core_votu as cv

def verify_catalog(root):
    root=Path(root)
    m=json.loads((root/'catalog_manifest.json').read_text())
    inventory=m.get('outputs',m.get('output_files',{}))
    if not inventory:raise ValueError('Catalogue lacks immutable file inventory')
    for name,digest in inventory.items():
        p=(root/name).resolve()
        if not p.is_relative_to(root.resolve()) or not p.is_file() or sha(p)!=digest:
            raise ValueError('Changed catalogue artifact: '+name)
    return m

def execute(a):
    import spark_runtime
    spark_runtime.allocation(require_lock=True)
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Sequence bridge requires a controlled Slurm allocation')
    core=Path(a.research_catalog).resolve(); v1=Path(a.independent_catalog).resolve()
    outpath=Path(a.outdir).resolve()
    for p in [core,v1,Path(a.source_manifest).resolve()]:
        if outpath==p or outpath.is_relative_to(p) or p.is_relative_to(outpath):
            raise ValueError('Output overlaps bridge input')
    cm=verify_catalog(core); vm=verify_catalog(v1)
    if not cm['catalog_version'].startswith('research_'):
        raise ValueError('Dedicated bridge requires research anchor/fragment catalogue')
    members=read_table(core/'vOTU_members.tsv',['candidate_id','votu_id','catalog_tier','member_role'])
    originals=unique(read_table(core/'prophage_master.tsv',['candidate_id','genome_id','sequence_sha256']),'candidate_id')
    eligibility=unique(read_table(core/'catalog_eligibility.tsv',['candidate_id','host_qc_eligible']),'candidate_id')
    seq=fasta(core/'candidates.fna')
    if set(originals)!=set(seq) or {x['candidate_id'] for x in members}!=set(seq):
        raise ValueError('Research sequence/master/member mismatch')
    unit_rows=load_units(a.source_manifest,verify=not bool(a.source_fastas))
    if a.source_fastas:
        if len(a.source_fastas)!=len(unit_rows):raise ValueError('Staged source count mismatch')
        for u,p in zip(unit_rows,a.source_fastas):
            if sha(p)!=u['sha256']:raise ValueError('Staged source hash mismatch')
            u['assembly_fasta']=str(Path(p).resolve())
    units=unique(unit_rows,'assembly_id')
    source_cache={}
    occurrences=[]; sid_by_cid={}; by_sid=defaultdict(list)
    for cid,h in originals.items():
        gid=h['genome_id']
        if gid not in units:raise ValueError('Source genome manifest missing '+gid)
        u=units[gid]
        if u['input_kind']!='bacterial_genome':raise ValueError('Research source must be bacterial genome')
        if gid not in source_cache:
            source_cache[gid]={k:(k,s) for k,s in fasta(u['assembly_fasta']).items()}
        hh=dict(h,source_host_genome=gid,main_catalog_eligible=eligibility[cid]['host_qc_eligible'])
        occ=make_occurrence(hh,seq[cid],u,source_cache[gid],legacy=True)
        occ['source_kind']='bacterial_genome'
        occ['source_genome_eligibility']='eligible' if eligibility[cid]['host_qc_eligible']=='true' else 'unknown'
        occurrences.append(occ);sid_by_cid[cid]=occ['viral_sequence_id'];by_sid[occ['viral_sequence_id']].append(cid)
    out=new_output(outpath)
    mapped=[]; seqs={}; master=[]
    membership=unique(members,'candidate_id')
    for sid,cids in sorted(by_sid.items()):
        values={seq[c] for c in cids}
        if len(values)!=1:raise ValueError('Hash collision or inconsistent sequence')
        seqs[sid]=values.pop()
        memberships=[membership[c] for c in cids]
        complete=[m for m in memberships if m['member_role']=='anchor' and m['complete_element_eligible']=='true' and m['votu_id']]
        scopes={o['virus_scope'] for o in occurrences if o['viral_sequence_id']==sid}-{'unknown'}
        if len(scopes)>1:raise ValueError('Conflicting viral scope for identical research sequence')
        master.append(dict(viral_sequence_id=sid,candidate_id=sid,length=len(seqs[sid]),
            sequence_sha256=hashlib.sha256(seqs[sid].encode()).hexdigest(),virus_scope=next(iter(scopes),'unknown'),
            viral_identity_confidence='reviewed_research_anchor' if complete else 'research_fragment_or_unresolved',
            n_source_occurrences=len(cids),viral_catalog_eligibility='eligible' if complete else 'unknown',
            eligibility_reason='upstream_reviewed_research_policy_not_V1_dev_cutoffs' if complete else 'fragment_not_anchor',
            eligibility_policy=cm['policy']['policy_version'],catalog_scope='host_resolved_prophage',
            completeness=originals[cids[0]].get('completeness','unknown'),
            contamination=originals[cids[0]].get('contamination','unknown'),
            checkv_quality=originals[cids[0]].get('checkv_quality','unknown')))
    for m in members:
        mapped.append(dict(m,source_candidate_id=m['candidate_id'],candidate_id=sid_by_cid[m['candidate_id']],
            viral_sequence_id=sid_by_cid[m['candidate_id']],representative_id=sid_by_cid.get(m['representative_id'],''),
            matching_representative_ids=';'.join(sid_by_cid[r] for r in m.get('matching_representative_ids','').split(';') if r),
            catalog_scope='host_resolved_prophage'))
    reps=[m for m in mapped if m['membership_status']=='representative']
    if len({m['votu_id'] for m in reps})!=len(reps):raise ValueError('Duplicate research vOTU representative')
    fields=list(dict.fromkeys(MEMBER_FIELDS+list(members[0] if members else {})+['source_candidate_id']))
    write_table(out/'vOTU_members.tsv',mapped,fields)
    write_table(out/'vOTU_representatives.tsv',reps,fields)
    write_table(out/'source_occurrences.tsv',occurrences,OCC_FIELDS)
    write_table(out/'viral_sequence_master.tsv',master,SEQ_FIELDS)
    write_table(out/'research_candidate_sequence_crosswalk.tsv',
        [dict(source_candidate_id=cid,viral_sequence_id=sid,sequence_sha256=originals[cid]['sequence_sha256']) for cid,sid in sorted(sid_by_cid.items())],
        ['source_candidate_id','viral_sequence_id','sequence_sha256'])
    from mainline_common import write_fasta
    write_fasta(out/'viral_sequences.fna',seqs)
    write_fasta(out/'vOTU_representatives.fna',{m['candidate_id']:seqs[m['candidate_id']] for m in reps})
    # Stable catalogue membership is not collapsed even for sequence-equivalent instances.
    om=[dict(o,**{k:v for k,v in membership[o['candidate_id']].items() if k not in o},
             source_candidate_id=o['candidate_id']) for o in occurrences]
    write_table(out/'occurrence_votu_membership.tsv',om,list(dict.fromkeys(OCC_FIELDS+list(members[0] if members else {})+['source_candidate_id'])))
    queries=fasta(v1/'viral_sequences.fna')
    q={f'q_{i}':s for i,s in enumerate(queries.values())}
    qids=dict(zip(q,queries))
    targets={f'r_{i}':seqs[m['candidate_id']] for i,m in enumerate(reps)}
    tids={f'r_{i}':m for i,m in enumerate(reps)}
    write_fasta(out/'bridge_queries.fna',q);write_fasta(out/'bridge_anchors.fna',targets)
    grouped={}
    if q and targets:
        for tag,query,subject in [('forward','bridge_queries.fna','bridge_anchors.fna'),('reverse','bridge_anchors.fna','bridge_queries.fna')]:
            run(['blastn','-task','blastn','-query',out/query,'-subject',out/subject,
                 '-num_threads',str(a.cpus),'-evalue','1e-5','-max_target_seqs',str(max(len(q),len(targets))),
                 '-outfmt','6 '+cv.BLAST_FIELDS,'-out',out/(tag+'.tsv')],out,tag,timeout=a.timeout)
            grouped.update(cv.parse_blast(out/(tag+'.tsv'),dict(q,**targets)))
    links=[]
    for qid,sid in qids.items():
        matched=[]
        for rid,rep in tids.items():
            exact=q[qid]==targets[rid]
            f=cv.aggregate_hsps(grouped.get((qid,rid),[]));r=cv.aggregate_hsps(grouped.get((rid,qid),[]))
            ani=100 if exact else min(f['ani_percent'],r['ani_percent'])
            af=100 if exact else min(f['af_shorter_percent'],r['af_shorter_percent'])
            if ani>=95 and af>=85:matched.append((rep,ani,af,exact))
        if not matched:
            links.append(dict(independent_catalog_version=vm['catalog_version'],independent_viral_sequence_id=sid,
                research_catalog_version=cm['catalog_version'],relationship='unresolved_no_anchor',annotation_transfer='false',host_transfer='false'))
        for rep,ani,af,exact in matched:
            links.append(dict(independent_catalog_version=vm['catalog_version'],independent_viral_sequence_id=sid,
                research_catalog_version=cm['catalog_version'],research_votu_id=rep['votu_id'],
                research_representative_id=rep['source_candidate_id'],catalog_tier=rep['catalog_tier'],
                relationship='ambiguous_multiple_anchors' if len(matched)>1 else 'exact_sequence' if exact else 'reciprocal_ANI95_AF85',
                ani_percent=ani,af_shorter_percent=af,annotation_transfer='false',host_transfer='false'))
    write_table(out/'independent_V1_links.tsv',links,['independent_catalog_version','independent_viral_sequence_id','research_catalog_version',
        'research_votu_id','research_representative_id','catalog_tier','relationship','ani_percent','af_shorter_percent','annotation_transfer','host_transfer'])
    save(out/'catalog_manifest.json',dict(schema_version='research_V1_bridge_v1',status='completed',
        catalog_version=cm['catalog_version'],catalog_scope='host_resolved_prophage',discovery_samples=sorted({u['biological_sample_id'] for u in unit_rows}-{'unknown'}),
        discovery_provenance_status='provided' if all(u['biological_sample_id']!='unknown' for u in unit_rows) else 'unknown',
        research_manifest_sha256=sha(core/'catalog_manifest.json'),independent_manifest_sha256=sha(v1/'catalog_manifest.json'),
        source_manifest_sha256=sha(a.source_manifest),source_gene_namespace='canonical_preserved_no_transfer',
        independent_gene_namespace='independent_V1_unchanged',reclustered=False,policy=cm['policy'],
        scientific_calibration=False,outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['research-catalog','independent-catalog','source-manifest','outdir']:p.add_argument('--'+name,required=True)
    p.add_argument('--source-fastas',nargs='*')
    p.add_argument('--cpus',type=int,default=1);p.add_argument('--timeout',type=int,default=3600)
    a=p.parse_args()
    if not 1<=a.cpus<=6 or not 1<=a.timeout<=3600:p.error('Resource cap exceeded')
    execute(a)
if __name__=='__main__':main()
