import os,sys,json,csv,pathlib,hashlib,subprocess,time,collections
import numpy as np
import pysam
R=pathlib.Path('/srv/CRC-PHIRE/analysis1');A=R/'14_postpilot_adjustment';W=R/'work/postpilot_adjustment/replication'
SAM='/srv/CRC-PHIRE/projects/crc-pipeline/envs/spark-reads-assembly/.pixi/envs/default/bin/samtools'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def table(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def save(p,d):p.write_text(json.dumps(d,indent=2))
def seqs(path,wanted):
    result={};name=None;buf=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name in wanted:result[name]=''.join(buf)
                name=line[1:].split()[0];buf=[]
            elif name in wanted:buf.append(line.strip())
    if name in wanted:result[name]=''.join(buf)
    return result
def main():
    W.mkdir(parents=True,exist_ok=True); output=A/'05_replication_activity'
    ev=R/'12_reports/pilot10/prophage_evidence/occurrence_evidence.tsv'
    rows=[r for r in table(ev) if r['source_kind']=='metagenome_assembly' and r['candidate_type']=='provirus_locus']
    assert len(rows)==382 and len({r['source_occurrence_id'] for r in rows})==382
    frozen=json.loads((R/'01_manifests/pilot10_frozen_assemblies.json').read_text());units={u['biological_sample_id']:u for u in frozen['units']}
    env={**os.environ,'PATH':str(pathlib.Path(SAM).parent)+':'+os.environ['PATH'],'NUMBA_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    results=[];bins=[];receipts=[]
    for sid in sorted(units):
        t=time.time();u=units[sid];loci=[r for r in rows if r['biological_sample_id']==sid];out=output/sid;out.mkdir(exist_ok=True);work=W/sid;work.mkdir(exist_ok=True)
        if (out/'completed.json').exists():
            rec=json.loads((out/'completed.json').read_text());assert sha(out/'loci.tsv')==rec['loci_sha256']
            results.extend(table(out/'loci.tsv'));bins.extend(table(out/'coverage_bins.tsv'));receipts.append(rec);continue
        fa=pathlib.Path(u['assembly_fasta']);assert sha(fa)==u['sha256']
        bam=fa.with_name('readback.bam');arc=json.loads((R/f'13_handoff/archive_{sid}.json').read_text());ar=[r for r in arc['files'] if r['path']==str(bam)]
        assert len(ar)==1 and ar[0]['status']=='archived_verified_link' and bam.stat().st_size==ar[0]['bytes']
        assert sha(bam)==ar[0]['sha256']
        index=R/f'12_reports/pilot10/candidate_review/bam_indexes/{sid}.bai'
        if not index.is_file():
            index=work/'original.bai';subprocess.run([SAM,'index','-@','1','-o',str(index),str(bam)],check=True)
        wanted={r['original_contig_id'] for r in loci};sequences=seqs(fa,wanted);assert set(sequences)==wanted
        usable=[]
        for r in loci:
            a,b=int(r['start0']),int(r['end0']);s=sequences[r['original_contig_id']]
            assert 0<=a<b<=len(s) and len(s)==int(r['contig_length'])
            assert hashlib.sha256(s[a:b].encode()).hexdigest()==r['sequence_sha256']
            if b-a>=1000:usable.append(r)
        # Full original contigs, original alignments; subset solely avoids scanning unrelated contigs.
        contigs={r['original_contig_id'] for r in usable};subset=work/'native_selected_contigs.bam'
        bfile=pysam.AlignmentFile(str(bam),'rb',index_filename=str(index));order=sorted(contigs,key=bfile.get_tid)
        strict={};native={};alignment_audit=[]
        with pysam.AlignmentFile(str(subset),'wb',template=bfile) as dst:
            for c in order:
                n=len(sequences[c]);d=np.zeros(n+1);q=np.zeros(n+1);nread=0;nm_missing=0;strict_reads=0
                for read in bfile.fetch(c):
                    dst.write(read);nread+=1
                    if not read.query_length or read.reference_end is None:continue
                    if not read.has_tag('NM'):nm_missing+=1
                    nm=read.get_tag('NM') if read.has_tag('NM') else 0
                    if nm/read.query_length>0.03:continue
                    a,b=read.reference_start,read.reference_end;d[a]+=1;d[b]-=1
                    if read.flag & 3844 or read.mapping_quality<20 or not read.has_tag('NM'):continue
                    strict_reads+=1
                    for a,b in read.get_blocks():q[a]+=1;q[b]-=1
                native[c]=np.cumsum(d[:-1]);strict[c]=np.cumsum(q[:-1]);alignment_audit.append({'contig':c,'alignment_records':nread,'missing_NM':nm_missing,'strict_records':strict_reads})
        bfile.close();subprocess.run([SAM,'index',str(subset)],check=True)
        selectedfa=work/'native_selected_contigs.fna'
        with open(selectedfa,'w') as f:
            for c in order:f.write('>'+c+'\n'+sequences[c]+'\n')
        coords=out/'coordinates.tsv';write(coords,[{'scaffold':r['original_contig_id'],'fragment':r['source_occurrence_id'],'start':int(r['start0'])+1,'stop':int(r['end0'])} for r in usable])
        cmd=[sys.executable,str(A/'02_pipeline/PropagAtE/Propagate/Propagate'),'-f',str(selectedfa),'-b',str(subset),'-v',str(coords),'-o',str(out/'propagate'),'-t','2','-p','0.97','-c','2.0','-e','0.70','--mask','150','--min','1.0','--breadth','0.5']
        save(out/'running.json',{'command':cmd,'source_bam':str(bam),'source_sha256':ar[0]['sha256']})
        with open(out/'propagate_console.log','w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,env=env)
        candidates=list((out/'propagate').glob('*.tsv'));assert len(candidates)==1,candidates
        tool={r['prophage']:r for r in table(candidates[0])};assert set(tool)=={r['source_occurrence_id'] for r in usable}
        sample=[];samplebins=[]
        for r in loci:
            oid=r['source_occurrence_id'];c=r['original_contig_id'];a,b=int(r['start0']),int(r['end0']);n=len(sequences[c]);tr=tool.get(oid,{})
            rec={k:r[k] for k in ['source_occurrence_id','biological_sample_id','votu_id','membership_status','original_contig_id','start0','end0','contig_length','boundary_confidence','two_sided_host_marker_evidence']};rec['group']=frozen['groups'][sid];rec['tool_active']=tr.get('active','not_run_short_provirus')
            for k in ['CohenD','prophage-host_ratio','prophage_mean_cov','prophage_cov_breadth','host_mean_cov','host_len']:rec['tool_'+k]=tr.get(k,'NA')
            for k in ['host_usable_bp','host_breadth','strict_phage_mean','strict_host_mean','strict_ratio','local_host_bp','local_ratio','sensitivity_agreement']:rec[k]='NA'
            rec['status']='not_evaluable';rec['reason']='provirus_under_1000bp';rec['replication_call']='not_evaluable'
            if oid in tool:
                d=native[c].copy();q=strict[c].copy();d[:150]=np.nan;d[-150:]=np.nan;q[:150]=np.nan;q[-150:]=np.nan
                hm=np.isfinite(d)
                for other in loci:
                    if other['original_contig_id']==c:hm[int(other['start0']):int(other['end0'])]=False
                hm_len=int(hm.sum());hb=float(np.mean(d[hm]>=1)) if hm_len else 0
                pm=float(np.nanmean(d[a:b]));host=float(np.nanmean(d[hm])) if hm_len else 0
                strictpm=float(np.nanmean(q[a:b]));stricth=float(np.nanmean(q[hm])) if hm_len else 0
                local=hm.copy();local[:max(0,a-5000)]=False;local[min(n,b+5000):]=False
                lh=float(np.nanmean(d[local])) if local.any() else 0
                rec.update(host_usable_bp=hm_len,host_breadth=hb,strict_phage_mean=strictpm,strict_host_mean=stricth,strict_ratio=strictpm/stricth if stricth>0 else 'NA',local_host_bp=int(local.sum()),local_ratio=pm/lh if lh>0 else 'NA')
                # Tool-result agreement checks identify format/missing-last-contig errors.
                if tr['prophage_mean_cov']!='NA':assert np.isclose(pm,float(tr['prophage_mean_cov']),rtol=1e-6,atol=1e-6),(sid,oid,pm,tr)
                reasons=[]
                if hm_len<1000:reasons.append('host_under_1000bp')
                if host<1 or hb<0.5:reasons.append('insufficient_host_coverage')
                if pm<1 or float(tr['prophage_cov_breadth'])<0.5:reasons.append('insufficient_provirus_coverage')
                if not reasons:
                    rec['status']='evaluable';rec['reason']='coverage_qualified_caller_boundary';rec['replication_call']='replication_compatible' if tr['active']=='active' else 'no_elevated_replication_evidence'
                    if stricth>0 and int(local.sum())>=1000 and lh>0:
                        sd=np.nanstd(q[a:b]);hs=np.nanstd(q[hm]);pool=((sd*sd+hs*hs)/2)**.5
                        strictactive=strictpm/stricth>=2 and pool>0 and abs(strictpm-stricth)/pool>=.7 and strictpm>=1 and np.mean(q[a:b]>=1)>=.5
                        rec['sensitivity_agreement']=bool((tr['active']=='active')==strictactive and (pm/lh>=2)==(float(tr['prophage-host_ratio'])>=2))
                else:rec['reason']=';'.join(reasons)
                for lo in range(max(0,a-5000),min(n,b+5000),250):
                    hi=min(lo+250,n,b+5000);samplebins.append({'source_occurrence_id':oid,'sample':sid,'start0':lo,'end0':hi,'relative_start':lo-a,'native_depth':float(np.nanmean(d[lo:hi])) if np.isfinite(d[lo:hi]).any() else 'NA','strict_depth':float(np.nanmean(q[lo:hi])) if np.isfinite(q[lo:hi]).any() else 'NA','provirus_start0':a,'provirus_end0':b})
            sample.append(rec)
        write(out/'loci.tsv',sample);write(out/'coverage_bins.tsv',samplebins);write(out/'alignment_audit.tsv',alignment_audit)
        rec={'status':'completed','sample':sid,'loci':len(sample),'loci_sha256':sha(out/'loci.tsv'),'elapsed_seconds':time.time()-t,'source_bam_sha256':ar[0]['sha256'],'source_fasta_sha256':u['sha256'],'command':cmd,'job_id':os.environ.get('SLURM_JOB_ID')};save(out/'completed.json',rec)
        receipts.append(rec);results.extend(sample);bins.extend(samplebins);save(output/'progress.json',receipts);print(sid,len(sample),round(time.time()-t),flush=True)
    write(output/'locus_diagnostics.tsv',results);write(output/'coverage_bins.tsv',bins)
    save(output/'completed.json',{'status':'completed','loci':len(results),'donors':receipts,'input_evidence_sha256':sha(ev),'result_sha256':sha(output/'locus_diagnostics.tsv'),'status_counts':dict(collections.Counter(r['status'] for r in results)),'call_counts':dict(collections.Counter(r['replication_call'] for r in results))})
if __name__=='__main__':main()
