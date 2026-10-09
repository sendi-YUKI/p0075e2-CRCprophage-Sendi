"""Validate converged forward/reverse lanes, aggregate once, diagnose host measurability.
Scheduled after BOTH lanes terminate. Never starts mapping or overwrites running outputs.
"""
import pathlib,json,csv,collections,re,os,subprocess,time,sys
from host_profile import R,A,sha,save,parse_metaphlan

def table(p):
    with open(p,encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    assert rows
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    root=A/'03_host_profile'; identity=table(A/'01_manifests/donor_input_identity.tsv')
    ids={r['sample']:r for r in identity};assert len(ids)==10
    live=json.loads(subprocess.check_output(['spark-status','--json']))
    assert not any(str(j['id']) in ['1133','1139'] for j in live['my_jobs']), 'Wait until both producers leave queue'
    units=json.loads((R/'01_manifests/pilot10_clean_reads.json').read_text())['units']
    expected=collections.defaultdict(dict)
    for u in units:
        for k in ['reads_1','reads_2']:
            if u.get(k):
                assert u[k] not in expected[u['biological_sample_id']]
                expected[u['biological_sample_id']][u[k]]=u['checksums'][u[k]]
    allrows=[];receipts=[];qc=[];recovered=[]
    for sid in sorted(ids):
        out=root/sid; recpath=out/'completed.json'
        if not recpath.exists():
            tail=R/'work/postpilot_adjustment/host_tail'/sid
            if not out.exists() and (tail/'completed.json').exists():
                r=json.loads((tail/'completed.json').read_text());assert sha(tail/'profile.tsv')==r['profile_sha256']
                tail.rename(out);recovered.append(sid)
        assert recpath.is_file(),f'Missing complete receipt: {sid}; retain all partial outputs and inspect'
        rec=json.loads(recpath.read_text());profile=out/'profile.tsv'
        assert rec['status']=='completed' and rec['sample']==sid and rec['group']==ids[sid]['group']
        assert sha(profile)==rec['profile_sha256']
        assert {r['path']:r['sha256'] for r in rec['input_checks']}==expected[sid]
        assert len(rec['input_checks'])==len(expected[sid])
        for r in rec['input_checks']:assert pathlib.Path(r['path']).stat().st_size==r['bytes']
        assert rec['clean_reads']==int(ids[sid]['clean_reads']) and rec['input_readsets']==int(ids[sid]['readsets'])
        parsed=parse_metaphlan(profile,rec['synthetic_donor_unit']);assert parsed
        assert all(r['biological_sample_id']==sid for r in parsed)
        keys=[r['clade_name'] for r in parsed];assert len(keys)==len(set(keys)), 'Duplicate clade records'
        text=profile.read_text();processed=re.search(r'^#(\d+) reads processed$',text,re.M);assert processed
        n=int(processed.group(1));assert 0<n<=rec['clean_reads']
        allrows.extend(parsed);receipts.append(rec)
        qc.append({'sample':sid,'group':rec['group'],'input_clean_reads':rec['clean_reads'],'profile_processed_reads':n,'processed_fraction':n/rec['clean_reads'],'species_rows':sum(r['taxon_level']=='species' for r in parsed),'job_id':rec['job_id'],'nproc':rec['command'][rec['command'].index('--nproc')+1],'elapsed_seconds':rec['elapsed_seconds'],'profile_sha256':rec['profile_sha256'],'input_identity_verified':True,'read_count_note':'MetaPhlAn default minimum read length 70; original clean input count includes shorter reads'})
    assert collections.Counter(r['group'] for r in qc)=={'CRC':5,'control':5}
    write(root/'donor_species.tsv',allrows);write(root/'host_profile_qc.tsv',qc)
    save(root/'completed.json',{'status':'completed','donors':receipts,'output_sha256':sha(root/'donor_species.tsv'),'database_index':'mpa_vJan25_CHOCOPhlAnSGB_202503','finalized_by_job':os.environ.get('SLURM_JOB_ID'),'recovered_tail_samples':recovered,'all_10_identity_and_profile_checks_passed':True})
    db=table(root/'database_taxonomy_coverage.tsv');values={}
    for r in allrows:
        if r['taxon_level']=='species':values[(r['biological_sample_id'],r['clade_name'].split('|')[-1][3:])]=float(r['relative_abundance_percent'])
    host=[];summary=[]
    for tax in db:
        name=tax['requested_taxon'];supported=tax['exact_species_in_database']=='True'
        for sid in sorted(ids):
            value=values.get((sid,name),0) if supported else 'NA'
            status=('reported_species_signal' if value>0 else 'not_reported_at_profiler_resolution') if supported else 'exact_species_label_unresolved_in_database'
            host.append({'sample':sid,'group':ids[sid]['group'],'taxon':name,'relative_abundance_percent':value,'measurement_status':status,'database_species_label_supported':supported,'SGB_taxonomy_records':tax['SGB_taxonomy_records'],'patient_phage_carriage_inference':'requires_linked_phage_evidence','strain_resolution':'not_established'})
        for group in ['CRC','control']:
            vals=[r['relative_abundance_percent'] for r in host if r['taxon']==name and r['group']==group]
            summary.append({'taxon':name,'group':group,'n_donors':5,'n_profile_signal':sum(v>0 for v in vals) if supported else 'NA','minimum_percent':min(vals) if supported else 'NA','maximum_percent':max(vals) if supported else 'NA','database_resolution':'species_label_available' if supported else 'exact_label_unresolved_not_biological_zero'})
    write(root/'priority_host_by_donor.tsv',host);write(root/'priority_host_summary.tsv',summary)
    # Header consistency is an independent postflight on the selected native BAMs.
    import pysam
    loci=table(A/'05_replication_activity/locus_diagnostics.tsv');audit=[]
    for sid in ids:
        bam=R/'work/postpilot_adjustment/replication'/sid/'native_selected_contigs.bam'
        with pysam.AlignmentFile(str(bam),'rb') as b:
            assert b.has_index() and b.header.to_dict().get('HD',{}).get('SO')=='coordinate'
            for r in loci:
                if r['biological_sample_id']!=sid:continue
                assert b.get_reference_length(r['original_contig_id'])==int(r['contig_length'])
                assert 0<=int(r['start0'])<int(r['end0'])<=int(r['contig_length'])
                audit.append({'occurrence_id':r['source_occurrence_id'],'sample':sid,'bam_contig_length_verified':True,'coordinates_valid':True})
    assert len(audit)==382;write(A/'05_replication_activity/coordinate_bam_postflight.tsv',audit)
    save(A/'09_handoff/host_convergence_acceptance.json',{'status':'passed','donors':10,'readsets':len(units),'fastq_inputs':sum(map(len,expected.values())),'coordinate_bam_loci':len(audit),'job_id':os.environ.get('SLURM_JOB_ID'),'host_table_sha256':sha(root/'donor_species.tsv'),'recovered_tail_samples':recovered,'large_fastq_sha_reused_from_per_donor_receipts_not_rehashed':True})
    subprocess.run([sys.executable,str(A/'02_pipeline/make_figures.py'),'--part','host'],check=True)
    print('Ten-donor convergence validated; host summaries and figure generated.',flush=True)
if __name__=='__main__':
    try:main()
    except Exception as e:
        save(A/'09_handoff/host_convergence_failure.json',{'status':'failed','error':repr(e),'job_id':os.environ.get('SLURM_JOB_ID'),'time_epoch':time.time()});raise
