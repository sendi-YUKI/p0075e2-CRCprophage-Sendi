"""Bounded complete-start SE/PE reference probes using the actual fragment kernel.

Perfect reference probes support a layout-specific accessible interval mask, not
biological calibration or a guarantee against missing references/variant reads.
"""
from collections import defaultdict,Counter
from pathlib import Path
from cohort_common import read_json,write_json,write_table,fasta,sha256,require,digest
from cohort_measure import sam_stream,sam_groups,classify_fragment

def ranges(values):
    out=[]
    for v in sorted(values):
        if out and out[-1][1]==v:out[-1][1]=v+1
        else:out.append([v,v+1])
    return out

def build(index_dir,profile_file,outdir,threads):
    from cohort_tools import map_reads
    index=Path(index_dir);out=Path(outdir);out.mkdir(parents=True,exist_ok=False)
    p=read_json(profile_file);policy=p['informative_policy'];ref=fasta(index/'reference.fna');m=read_json(index/'reference_manifest.json')
    require(policy.get('method')=='complete_start_reference_probes_v1','Unknown informative method')
    require(policy.get('status') in ('synthetic','proposed'),'Uncalibrated probe policy must remain proposed/synthetic')
    require(policy.get('layouts') and set(policy['layouts'])<= {'SE','PE'},'Invalid probe layouts')
    require(len(set(policy['layouts']))==len(policy['layouts']),'Duplicate probe layout')
    maximum=policy['max_fragments_total'];require(isinstance(maximum,int) and 1<=maximum<=50000,'Probe preparation cap1..50000')
    L=p['mappability_read_length'];F=policy.get('pe_fragment_length')
    if 'PE' in policy['layouts']:
        require(p['pair_orientation']=='fr' and isinstance(F,int) and max(L,p['insert_min'])<=F<=p['insert_max'],'PE probe domain requires explicit FR fragment length')
    targets=[r for r in m['references'] if r['reference_kind']=='target'];rows={r['reference_id']:{'reference_id':r['reference_id'],'layouts':{}} for r in targets}
    total=0;details=[]
    rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    for layout in policy['layouts']:
        span=L if layout=='SE' else F;folder=out/layout;folder.mkdir();meta={};positions=defaultdict(set);good=Counter();tested=Counter()
        f1=folder/'probe_R1.fastq';f2=folder/'probe_R2.fastq'
        with f1.open('w') as h1,f2.open('w') as h2:
            for t in targets:
                name=t['reference_id'];ss=ref[name];potential=max(0,len(ss)-span+1)*2
                state='complete' if potential and total+potential<=maximum else 'not_assessed_too_short_or_probe_cap'
                rr=dict(status=state,read_length=L,fragment_length=span if layout=='PE' else None,potential_fragments=potential,
                    tested_fragments=0,unique_tile_fraction=None,informative_intervals=[],informative_length=0)
                rows[name]['layouts'][layout]=rr
                if state!='complete':continue
                total+=potential
                for start in range(len(ss)-span+1):
                    frag=ss[start:start+span]
                    for direction in ('+','-'):
                        if not set(frag)<=set('ACGT'):continue
                        oriented=frag if direction=='+' else rc(frag)
                        ident='probe_'+digest([name,start,layout,direction,L,span])[:28]
                        meta[ident]=(name,start,direction)
                        a=oriented[:L];h1.write('@'+ident+'\n'+a+'\n+\n'+'I'*L+'\n')
                        if layout=='PE':h2.write('@'+ident+'\n'+rc(oriented[-L:])+'\n+\n'+'I'*L+'\n')
        if meta:
            unit=dict(unit_id='probe_'+layout,layout=layout,reads_1=str(f1),reads_2=str(f2) if layout=='PE' else None,counts={'fragments':len(meta)},checksums={str(f1):sha256(f1)})
            if layout=='PE':unit['checksums'][str(f2)]=sha256(f2)
            write_json(folder/'unit.json',unit)
            bam=map_reads(index,folder/'unit.json',profile_file,f1,f2 if layout=='PE' else None,folder/'mapping',threads)
            from itertools import zip_longest
            with sam_stream(bam) as audit,sam_stream(folder/'mapping/alignments.primary.qname.bam') as primary:
                observed=set()
                for ag,pg in zip_longest(sam_groups(audit),sam_groups(primary)):
                    require(ag is not None and pg is not None and ag[0]==pg[0] and ag[0] in meta,'Probe SAM identity mismatch')
                    require(ag[0] not in observed,'Repeated probe QNAME');observed.add(ag[0])
                    name,start,direction=meta[ag[0]]
                    state,target,covered,alt=classify_fragment(ag[1],ref,p,layout,pg[1]);tested[name]+=1
                    is_self=state=='assigned' and target==name
                    if is_self:good[name]+=1;positions[name].update(covered)
                    details.append(dict(probe_id=ag[0],reference_id=name,layout=layout,start0=start,direction=direction,
                        assigned_to_origin=is_self,reason=state,alternative_reference_ids=';'.join(alt)))
                require(observed==set(meta),'Missing probe SAM records')
        for name,row in rows.items():
            rr=row['layouts'][layout];rr.update(tested_fragments=tested[name],unique_tile_fraction=good[name]/tested[name] if tested[name] else None,
                informative_intervals=ranges(positions[name]),informative_length=len(positions[name]))
            if not tested[name] and rr['status']=='complete':rr['status']='not_assessed_no_unambiguous_probe_sequence'
    fields=['probe_id','reference_id','layout','start0','direction','assigned_to_origin','reason','alternative_reference_ids']
    write_table(out/'probe_evidence.tsv',fields,details)
    result=dict(schema_version=2,algorithm='complete_start_reference_probes_v1',measurement_profile=p,
        mapping_reference_version=m['mapping_reference_version'],reference_fasta_sha256=sha256(index/'reference.fna'),
        references=list(rows.values()),fragments_budgeted=total,fragments_actually_mapped=len(details),
        scope='declared_reference_perfect_reads_fixed_length_FR_insert_domain; heuristic_mapper_not_exhaustive_alignment_search',
        grouping='exact_equivalence_only; nonidentical_ambiguous_targets_remain_unresolved_no_transitive_merge',
        scientific_calibration=False,outputs={'probe_evidence.tsv':sha256(out/'probe_evidence.tsv')})
    write_json(out/'mappability.json',result)
    return result
