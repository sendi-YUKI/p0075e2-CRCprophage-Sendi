import pathlib,json,hashlib,subprocess,os,time,csv,sys,collections
R=pathlib.Path('/srv/CRC-PHIRE/analysis1'); A=R/'14_postpilot_adjustment'
SRC=R/'02_pipeline/src_pilot6_measurement'
sys.path.insert(0,str(SRC/'bin'))
from cohort_tools import metaphlan_preflight,parse_metaphlan
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def save(p,d): p.write_text(json.dumps(d,indent=2))
def main():
    start=time.time(); config=json.loads((SRC/'assets/spark_metaphlan_database.json').read_text())
    db=config['database_dir']; index=config['index']
    adapted={**config,'database_files':{pathlib.Path(r['path']).name:r['sha256'] for r in config['files']}}
    dbmanifest=A/'01_manifests/metaphlan_database_adapter.json';save(dbmanifest,adapted)
    metaphlan_preflight(db,index,dbmanifest)
    units=json.loads((R/'01_manifests/pilot10_clean_reads.json').read_text())['units']
    donors=collections.defaultdict(list)
    for u in units: donors[u['biological_sample_id']].append(u)
    order=['SID31866','SID31160']+sorted(set(donors)-{'SID31866','SID31160'})
    prefix='/srv/CRC-PHIRE/projects/crc-pipeline/envs/spark-reads-metaphlan/.pixi/envs/default'
    env={**os.environ,'PATH':prefix+'/bin:'+os.environ['PATH'],'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    rows=[]; manifest=[]
    for sid in order:
        out=A/'03_host_profile'/sid;out.mkdir(exist_ok=True)
        if (out/'completed.json').exists():
            rec=json.loads((out/'completed.json').read_text());assert sha(out/'profile.tsv')==rec['profile_sha256']
            rows.extend(parse_metaphlan(out/'profile.tsv',rec['synthetic_donor_unit']));manifest.append(rec);continue
        begin=time.time(); files=[];checks=[]
        for u in donors[sid]:
            for key in ['reads_1','reads_2']:
                if u.get(key):
                    p=pathlib.Path(u[key]);assert p.is_file() and str(p.resolve()) not in [x['resolved'] for x in checks], 'Duplicate physical clean read input'
                    observed=sha(p);assert observed==u['checksums'][str(p)],'Clean FASTQ hash changed'
                    files.append(str(p));checks.append({'path':str(p),'resolved':str(p.resolve()),'sha256':observed,'bytes':p.stat().st_size})
        command=[prefix+'/bin/metaphlan',','.join(files),'--input_type','fastq','--db_dir',db,'--index',index,'--offline','--nproc',os.environ.get('SLURM_CPUS_PER_TASK','6'),'--bt2_ps','very-sensitive','--min_mapq_val','5','--mapout',str(out/'markers.mapout.bz2'),'-t','rel_ab','-o',str(out/'profile.tsv')]
        save(out/'running.json',{'command':command,'input_checks':checks,'started_epoch':begin,'job_id':os.environ.get('SLURM_JOB_ID')})
        with open(out/'metaphlan.log','w') as log: subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
        unit={'unit_id':sid,'biological_sample_id':sid,'library_id':'donor_combined_disjoint_clean_reads','eligibility':{'host_adjustment':{'status':'provisional'}}}
        parsed=parse_metaphlan(out/'profile.tsv',unit);assert parsed
        rows.extend(parsed)
        rec={'status':'completed','sample':sid,'group':donors[sid][0]['metadata']['disease_group'],'command':command,'input_checks':checks,'input_readsets':len(donors[sid]),'clean_reads':sum(u['counts']['reads'] for u in donors[sid]),'elapsed_seconds':time.time()-begin,'profile_sha256':sha(out/'profile.tsv'),'synthetic_donor_unit':unit,'job_id':os.environ.get('SLURM_JOB_ID')}
        save(out/'completed.json',rec);manifest.append(rec);save(A/'03_host_profile/progress.json',manifest)
        print(sid,'completed',round(rec['elapsed_seconds']),flush=True)
    with open(A/'03_host_profile/donor_species.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    save(A/'03_host_profile/completed.json',{'status':'completed','donors':manifest,'elapsed_seconds':time.time()-start,'output_sha256':sha(A/'03_host_profile/donor_species.tsv'),'database_index':index})
if __name__=='__main__':main()
