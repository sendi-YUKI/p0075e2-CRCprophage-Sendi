import pathlib,json,csv,collections,random,subprocess,sys,os,hashlib,zipfile
import numpy as np
R=pathlib.Path('/srv/CRC-PHIRE/analysis1');A=R/'14_postpilot_adjustment';M=R/'09_measurement/pilot10';P=R/'12_reports/pilot10';OUT=A/'04_measurement_audit/technical_calibration';W=R/'work/postpilot_adjustment/calibration'
sys.path.insert(0,str(R/'02_pipeline/src_pilot6_measurement/bin'))
from cohort_measure import classify_fragment,sam_groups
from cohort_common import fasta
def table(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def save(p,d):p.write_text(json.dumps(d,indent=2))
def main():
    OUT.mkdir(exist_ok=True);W.mkdir(parents=True,exist_ok=True)
    rows=table(P/'donor_votu_measurement.tsv');refs=table(M/'prepared/reference/reference_map.tsv');byref={r['reference_id']:r for r in refs};sequences=fasta(M/'prepared/reference/reference.fna')
    agg=collections.defaultdict(lambda:collections.Counter())
    for r in rows:
        z=agg[r['reference_id']];z['detected']+=r['detection_status']=='detected';z['amb']+=float(r['ambiguous_fragment_count']);z['partial']+=0<float(r['breadth_full'])<.75;z['breadth']+=float(r['breadth_full'])
    taxa={r['votu_id']:r for r in table(A/'04_measurement_audit/votu_reference_taxa.tsv')};selected=[]
    categories=[('widely_detected',lambda r:agg[r]['detected']>=2,lambda r:-agg[r]['detected']),('ambiguous',lambda r:agg[r]['amb']>0,lambda r:-agg[r]['amb']),('partial_coverage',lambda r:agg[r]['partial']>0,lambda r:-agg[r]['partial']),('reference_only',lambda r:taxa[byref[r]['votu_id']]['has_metagenome_predicted_provirus']=='false' and taxa[byref[r]['votu_id']]['has_reference_predicted_provirus']=='true',lambda r:r),('F_animalis_link',lambda r:'animalis' in taxa[byref[r]['votu_id']]['reference_taxa'],lambda r:r),('P_micra_link',lambda r:'micra' in taxa[byref[r]['votu_id']]['reference_taxa'],lambda r:r)]
    targets=[r['reference_id'] for r in refs if r['reference_kind']=='target']
    for category,eligible,key in categories:
        choices=sorted([r for r in targets if eligible(r)],key=lambda r:(key(r),r))[:2]
        for r in choices:selected.append({'reference_id':r,'votu_id':byref[r]['votu_id'],'technical_category':category,'length':len(sequences[r])})
    chosen=sorted({r['reference_id'] for r in selected});assert len(chosen)>=6
    write(OUT/'selected_technical_examples.tsv',selected)
    # Frozen choices precede simulation. No disease P/effect enters selection.
    pairs=table(A/'04_measurement_audit/target_ambiguity_pairs.tsv');negative_targets=[]
    for r in pairs:
        x,y=r['target_1'],r['target_2']
        other=y if x in chosen and y not in chosen else x if y in chosen and x not in chosen else None
        if other and other not in negative_targets:negative_targets.append(other)
        if len(negative_targets)>=6:break
    hosts=sorted([r for r in refs if r['reference_kind']=='host' and len(sequences[r['reference_id']])>10000],key=lambda r:(-int(r['masked_bases']),r['reference_id']))
    seen=set();negative_hosts=[]
    for r in hosts:
        if r['host_id'] not in seen:
            negative_hosts.append(r['reference_id']);seen.add(r['host_id'])
        if len(negative_hosts)==6:break
    save(OUT/'known_sources.json',{'positive_targets':chosen,'negative_sources_for_selected_targets':negative_targets+negative_hosts,'negative_scope':'selected_targets_absent;other_catalog_viruses_and_masked_host_background_present','synthetic':True,'read_length':100,'PE_fragment_length':300,'substitution_rate':.005,'seed':20261008})
    profile=json.loads((M/'prepared/measurement_profile.json').read_text());bt='/srv/CRC-PHIRE/projects/crc-pipeline/envs/spark-reads-assembly/.pixi/envs/default/bin/bowtie2';index=M/'reference/mapping_index/index';summary=[];curves=[];commands=[];rng=random.Random(20261008)
    complement=str.maketrans('ACGT','TGCA')
    def err(s):return ''.join(rng.choice([x for x in 'ACGT' if x!=b]) if rng.random()<.005 else b for b in s)
    for layout,scenario in [('SE','positive'),('SE','negative'),('PE','positive'),('PE','negative')]:
        name=layout+'_'+scenario;f1=W/(name+'_1.fq');f2=W/(name+'_2.fq');sources=chosen if scenario=='positive' else negative_targets+negative_hosts;truth={};n=0
        with open(f1,'w') as a,open(f2,'w') as b:
            for source in sources:
                seq=sequences[source];span=100 if layout=='SE' else 300
                starts=list(range(0,len(seq)-span+1,33 if layout=='SE' else 67))
                if scenario=='negative' and source.startswith('host_') and len(starts)>5000:starts=sorted(rng.sample(starts,5000))
                for lo in starts:
                    segment=seq[lo:lo+span]
                    if any(x not in 'ACGT' for x in segment):continue
                    q=f'{name}_{n:09d}';n+=1;truth[q]=source;read=err(segment[:100]);a.write(f'@{q}\n{read}\n+\n'+('I'*100)+'\n')
                    if layout=='PE':read=err(segment[-100:].translate(complement)[::-1]);b.write(f'@{q}\n{read}\n+\n'+('I'*100)+'\n')
        for audit in [False,True]:
            dest=W/(name+('.audit.sam' if audit else '.primary.sam'));cmd=[bt,'--end-to-end','--very-sensitive','--reorder','-p','2','-x',str(index)]
            cmd+=['-U',str(f1)] if layout=='SE' else ['-1',str(f1),'-2',str(f2),'--fr','-I','0','-X','1000']
            if audit:cmd+=['-k','20']
            cmd+=['-S',str(dest)];commands.append(cmd)
            with open(OUT/(name+('.audit.log' if audit else '.primary.log')),'w') as log:subprocess.run(cmd,stderr=log,stdout=log,check=True)
        depth={r:np.zeros(len(sequences[r]),dtype=np.uint32) for r in chosen};counts=collections.Counter();wrong=collections.Counter();reasons=collections.Counter();amb=collections.Counter()
        with open(W/(name+'.primary.sam')) as pf,open(W/(name+'.audit.sam')) as af:
            for (q,primary),(q2,audit) in zip(sam_groups(pf),sam_groups(af),strict=True):
                assert q==q2
                reason,target,positions,alternatives=classify_fragment(audit,sequences,profile,layout,primary);reasons[reason]+=1
                if target in depth:
                    counts[target]+=1;wrong[target]+=truth[q]!=target
                    for pos in positions:depth[target][pos]+=1
                for target in alternatives:
                    if target in depth:amb[target]+=1
        assert sum(reasons.values())==n
        for target in chosen:
            d=depth[target];breadth=float(np.mean(d>=1));summary.append({'scenario':scenario,'layout':layout,'target':target,'true_present':scenario=='positive','source_fragments':sum(s==target for s in truth.values()),'assigned_fragments':counts[target],'wrong_source_assigned':wrong[target],'ambiguous_mentions':amb[target],'breadth':breadth,'detected_075':breadth>=.75 and counts[target]>0,'detected_030':breadth>=.30 and counts[target]>0,'synthetic':True})
            for lo in range(0,len(d),250):curves.append({'scenario':scenario,'layout':layout,'target':target,'start0':lo,'end0':min(lo+250,len(d)),'mean_depth':float(d[lo:lo+250].mean())})
        save(OUT/(name+'_filters.json'),dict(reasons));print(name,n,'completed',flush=True)
    write(OUT/'calibration_results.tsv',summary);write(OUT/'synthetic_coverage_bins.tsv',curves);save(OUT/'commands.json',commands)
    # Real donor curves for the SAME technically chosen objects, assembled from saved depth arrays.
    units=json.loads((R/'01_manifests/pilot10_clean_reads.json').read_text())['units'];real=collections.defaultdict(lambda:{r:np.zeros(len(sequences[r]),dtype=np.uint64) for r in chosen})
    for u in units:
        with zipfile.ZipFile(M/'units'/('unit_'+u['unit_id'])/'target_depths.zip') as z:
            for target in chosen:real[u['biological_sample_id']][target]+=np.frombuffer(z.read(target+'.u32le'),dtype='<u4')
    realcurves=[]
    for sid,targetsdepth in real.items():
        for target,d in targetsdepth.items():
            original=next(r for r in rows if r['reference_id']==target and r['biological_sample_id']==sid)
            assert np.isclose(np.mean(d>=1),float(original['breadth_full']))
            for lo in range(0,len(d),250):realcurves.append({'sample':sid,'target':target,'votu_id':byref[target]['votu_id'],'start0':lo,'end0':min(lo+250,len(d)),'mean_depth':float(d[lo:lo+250].mean())})
    write(OUT/'real_coverage_bins.tsv',realcurves);save(OUT/'completed.json',{'status':'completed','selected_targets':len(chosen),'synthetic_scenarios':4,'negative_interpretation':'bounded_challenge_not_global_specificity_or_real_disease_control','new_mapping_of_patients':False})
if __name__=='__main__':main()
