"""Use a second free slot after a verified first donor; publish completed tail samples for legacy driver's reuse.
No writes into a running sample directory. Temporary output is isolated; publication is atomic.
"""
import pathlib,json,os,subprocess,time,collections,sys
from host_profile import R,A,SRC,sha,save,parse_metaphlan,metaphlan_preflight
def main():
    root=A/'03_host_profile'
    # User authorized all currently free resources. The first control can finish
    # in the original lane while independent remaining donors run here.
    for sid in ['SID31866']:
        rec=json.loads((root/sid/'completed.json').read_text());assert rec['status']=='completed' and sha(root/sid/'profile.tsv')==rec['profile_sha256']
    units=json.loads((R/'01_manifests/pilot10_clean_reads.json').read_text())['units'];donors=collections.defaultdict(list)
    for u in units:donors[u['biological_sample_id']].append(u)
    order=['SID31866','SID31160']+sorted(set(donors)-{'SID31866','SID31160'})
    cfg=json.loads((A/'01_manifests/metaphlan_database_adapter.json').read_text());db=cfg['database_dir'];index=cfg['index']
    # Reuse completed first-donor database preflight and freeze, verify identity/size instead of rehashing 50 GB.
    for r in cfg['files']:assert pathlib.Path(r['path']).stat().st_size==r['bytes']
    prefix='/srv/CRC-PHIRE/projects/crc-pipeline/envs/spark-reads-metaphlan/.pixi/envs/default';env={**os.environ,'PATH':prefix+'/bin:'+os.environ['PATH'],'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    base=R/'work/postpilot_adjustment/host_tail';base.mkdir(parents=True,exist_ok=True);records=[]
    for sid in reversed(order[2:]):
        canonical=root/sid
        if canonical.exists():continue
        front=[i for i,s in enumerate(order) if not (root/s/'completed.json').exists()]
        if not front or order.index(sid)<min(front)+2:break
        out=base/sid
        if out.exists():raise RuntimeError('Inspect preserved partial helper output before retry: '+str(out))
        out.mkdir();begin=time.time();files=[];checks=[]
        for u in donors[sid]:
            for key in ['reads_1','reads_2']:
                if u.get(key):
                    p=pathlib.Path(u[key]);assert p.is_file() and str(p.resolve()) not in [x['resolved'] for x in checks]
                    h=sha(p);assert h==u['checksums'][str(p)];files.append(str(p));checks.append({'path':str(p),'resolved':str(p.resolve()),'sha256':h,'bytes':p.stat().st_size})
        if canonical.exists():raise RuntimeError('Front driver reached candidate; preserve input checks and stop without duplicate mapping')
        cmd=[prefix+'/bin/metaphlan',','.join(files),'--input_type','fastq','--db_dir',db,'--index',index,'--offline','--nproc',os.environ['SLURM_CPUS_PER_TASK'],'--bt2_ps','very-sensitive','--min_mapq_val','5','--mapout',str(out/'markers.mapout.bz2'),'-t','rel_ab','-o',str(out/'profile.tsv')]
        save(out/'running.json',{'command':cmd,'input_checks':checks,'started_epoch':begin,'job_id':os.environ.get('SLURM_JOB_ID'),'worker':'isolated_tail'})
        with open(out/'metaphlan.log','w') as log:
            proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            while proc.poll() is None:
                if canonical.exists():
                    import signal
                    os.killpg(proc.pid,signal.SIGTERM)
                    try:proc.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                    save(out/'interrupted.json',{'reason':'front_driver_reached_donor','sample':sid,'job_id':os.environ.get('SLURM_JOB_ID')})
                    raise RuntimeError('Front reached helper donor; helper stopped and preserved partial output')
                time.sleep(5)
            if proc.returncode:raise subprocess.CalledProcessError(proc.returncode,cmd)
        unit={'unit_id':sid,'biological_sample_id':sid,'library_id':'donor_combined_disjoint_clean_reads','eligibility':{'host_adjustment':{'status':'provisional'}}};assert parse_metaphlan(out/'profile.tsv',unit)
        rec={'status':'completed','sample':sid,'group':donors[sid][0]['metadata']['disease_group'],'command':cmd,'input_checks':checks,'input_readsets':len(donors[sid]),'clean_reads':sum(u['counts']['reads'] for u in donors[sid]),'elapsed_seconds':time.time()-begin,'profile_sha256':sha(out/'profile.tsv'),'synthetic_donor_unit':unit,'job_id':os.environ.get('SLURM_JOB_ID'),'worker':'isolated_tail','original_output_path':str(out),'published_output_path':str(canonical),'database_preflight_reuse':'same_fixed_database_as_two_verified_initial_donors_plus_size_check'}
        rec['database_preflight_reuse']='same_fixed_database_as_verified_initial_CRC_donor_plus_size_check'
        save(out/'completed.json',rec)
        # The front must still be >=2 samples away; otherwise leave a separate verified result for explicit reconciliation.
        front=[i for i,s in enumerate(order) if not (root/s/'completed.json').exists()]
        if canonical.exists() or (front and order.index(sid)<min(front)+2):
            rec['publication_status']='preserved_separate_front_approaching';records.append(rec);save(A/'09_handoff/host_tail_progress.json',records);break
        out.rename(canonical);rec['publication_status']='published_for_front_driver_hash_verified_reuse';records.append(rec);save(A/'09_handoff/host_tail_progress.json',records);print(sid,'published',round(time.time()-begin),flush=True)
    save(A/'09_handoff/host_tail_completed.json',{'status':'completed','records':records})
if __name__=='__main__':main()
