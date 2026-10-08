#!/usr/bin/env python3
"""R1 shared clean reads, explicit human fragment filtering and qualification."""
import argparse,gzip,hashlib,itertools,json,os,re,shutil,sqlite3,sys
from pathlib import Path
from reads_io import dump,sha,run,fastq,qname,fq_counts,write_tsv

def credible_sam(fields,p):
    if len(fields)<11:raise ValueError('Malformed SAM record')
    for key in ('min_identity','min_aligned_query_fraction'):
        if not 0<=p[key]<=1:raise ValueError('Human filtering thresholds use fractions in [0,1]')
    if not 0<=p['min_mapq']<=254:raise ValueError('Invalid human MAPQ threshold')
    flag=int(fields[1])
    if flag&4 or fields[5]=='*':return False
    cigar=re.findall(r'(\d+)([MIDNSHP=X])',fields[5]);lengths={c:sum(int(n) for n,k in cigar if k==c) for c in 'MIDNSHP=X'}
    if not cigar or ''.join(n+c for n,c in cigar)!=fields[5] or any(int(n)<=0 for n,c in cigar):raise ValueError('Invalid CIGAR')
    if any(lengths[k] for k in ('N','H','P')):raise ValueError('Unsupported decontamination CIGAR; cannot verify original read fraction')
    tags={}
    for item in fields[11:]:
        parts=item.split(':',2)
        if len(parts)!=3:raise ValueError('Malformed SAM tag')
        tags[parts[0]]=parts[2]
    if 'NM' not in tags:raise ValueError('Human decontamination requires NM edit distances')
    aligned=lengths['M']+lengths['=']+lengths['X'];den=aligned+lengths['I']+lengths['D']
    qlen=aligned+lengths['I']+lengths['S'];nm=int(tags['NM'])
    if nm<lengths['I']+lengths['D'] or nm>den:raise ValueError('NM edit distance contradicts CIGAR')
    if fields[9]!='*' and len(fields[9])!=qlen:raise ValueError('SAM sequence length contradicts CIGAR')
    identity=(den-nm)/den if den else 0
    af=(aligned+lengths['I'])/qlen if qlen else 0
    mapq=int(fields[4])
    mapq_pass=p['min_mapq']==0 or (mapq!=255 and mapq>=p['min_mapq'])
    return mapq_pass and identity>=p['min_identity'] and af>=p['min_aligned_query_fraction']

def remove_human(r1,r2,sam,out,p):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    # Validate uniqueness before filtering; two distinct fragments sharing a name
    # cannot be distinguished by SAM and must never be silently co-filtered.
    original=fq_counts(r1,r2)
    db=sqlite3.connect(out/'human_names.sqlite')
    db.execute('DROP TABLE IF EXISTS fragments')
    db.execute('CREATE TABLE fragments (name TEXT PRIMARY KEY, seen INTEGER NOT NULL DEFAULT 0, hit INTEGER NOT NULL DEFAULT 0)')
    try:
        for rec in fastq(r1):db.execute('INSERT INTO fragments(name) VALUES (?)',(qname(rec),))
        with open(sam) as handle:
            for line in handle:
                if line.startswith('@'):continue
                fields=line.rstrip('\r\n').split('\t')
                if len(fields)<11:raise ValueError('Malformed SAM record')
                name=re.sub(r'/[12]$','',fields[0]);flag=int(fields[1])
                seen=(1 if flag&64 else 0)|(2 if flag&128 else 0) if r2 else 1
                if r2 and seen not in (1,2):raise ValueError('Paired SAM lacks a unique mate flag')
                if not db.execute('SELECT 1 FROM fragments WHERE name=?',(name,)).fetchone():
                    raise ValueError('SAM QNAME does not belong to input FASTQ: '+name)
                db.execute('UPDATE fragments SET seen=seen|?,hit=MAX(hit,?) WHERE name=?',
                           (seen,int(credible_sam(fields,p)),name))
        required=3 if r2 else 1
        if db.execute('SELECT COUNT(*) FROM fragments WHERE seen != ?',(required,)).fetchone()[0]:
            raise ValueError('SAM lacks records for one or more input reads/mates')
        db.commit();a=out/'clean_1.fastq.gz';b=out/'clean_2.fastq.gz' if r2 else None
        removed=kept=removed_bases=retained_bases=0
        with gzip.open(a,'wt',compresslevel=1) as fa:
            fb=gzip.open(b,'wt',compresslevel=1) if b else None
            try:
                pairs=itertools.zip_longest(fastq(r1),fastq(r2)) if r2 else ((x,None) for x in fastq(r1))
                for x,y in pairs:
                    if x is None or (r2 and (y is None or qname(x)!=qname(y))):raise ValueError('Pairs lost synchronization')
                    bases=len(x[1])+(len(y[1]) if y else 0)
                    if db.execute('SELECT hit FROM fragments WHERE name=?',(qname(x),)).fetchone()[0]:
                        removed+=1;removed_bases+=bases;continue
                    fa.write('\n'.join(x)+'\n')
                    if fb:fb.write('\n'.join(y)+'\n')
                    kept+=1;retained_bases+=bases
            finally:
                if fb:fb.close()
    finally:db.close()
    return str(a),str(b) if b else None,{'removed_fragments':removed,'retained_fragments':kept,
        'removed_reads':removed*(2 if r2 else 1),'retained_reads':kept*(2 if r2 else 1),
        'removed_bases':removed_bases,'retained_bases':retained_bases,'input_counts':original,
        'rule':'any_credible_mate_removes_pair','sam_coverage_verified':True}

def prepare_viromeqc(script,local,virome_db):
    script=Path(script).resolve();local=Path(local);db=Path(virome_db).resolve()
    # These are the specific companion files consumed by pinned ViromeQC 1.0.2.
    required=['viromeQC.py','fastq_len_filter.py','medians.csv']
    for name in required:
        if not (script.parent/name).is_file():raise ValueError('ViromeQC companion missing: '+name)
    helper=script.parent/'cmseq/cmseq/filter.py'
    if not helper.is_file():raise ValueError('ViromeQC cmseq filter is missing')
    for prefix in ('SILVA_132_LSURef_tax_silva.clean','SILVA_132_SSURef_Nr99_tax_silva.clean'):
        for suffix in ('1.bt2','2.bt2','3.bt2','4.bt2','rev.1.bt2','rev.2.bt2'):
            path=db/(prefix+'.'+suffix)
            if not path.is_file() or path.stat().st_size==0:raise ValueError('ViromeQC database incomplete: '+str(path))
    protein=db/'amphora_bacteria_294.dmnd'
    if not protein.is_file() or protein.stat().st_size==0:raise ValueError('ViromeQC DIAMOND database missing')
    local.mkdir(parents=True,exist_ok=True)
    for name in required:shutil.copy2(script.parent/name,local/name)
    shutil.copytree(script.parent/'cmseq',local/'cmseq',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','.git'))
    idx=local/'index'
    if idx.exists() or idx.is_symlink():
        if not idx.is_symlink() or idx.resolve()!=db:raise ValueError('Unexpected ViromeQC runtime index path')
    else:idx.symlink_to(db,target_is_directory=True)
    # Upstream 1.0.2 does not check all children in marker pipelines. Guard the
    # fixed script copy so a failed aligner cannot be reported as zero contamination.
    runtime=local/script.name;source=runtime.read_text()
    replacements={
        'SSU_reads = int(p4.communicate()[0])':'SSU_reads = int(p4.communicate()[0])\n\tif any(proc.returncode != 0 for proc in (p1,p2,p3,p4)): raise RuntimeError("SSU marker subprocess failed")',
        'LSU_reads = int(p4.communicate()[0])':'LSU_reads = int(p4.communicate()[0])\n\tif any(proc.returncode != 0 for proc in (p1,p2,p3,p4)): raise RuntimeError("LSU marker subprocess failed")',
        'singleCopyMarkers_reads = int(p2.communicate()[0])':'singleCopyMarkers_reads = int(p2.communicate()[0])\n\tp1.wait()\n\tif any(proc.returncode != 0 for proc in (p1,p2)): raise RuntimeError("AMPhORA marker subprocess failed")'}
    for old,new in replacements.items():
        if source.count(old)!=1:raise ValueError('Pinned ViromeQC source differs from reviewed pipeline guards')
        source=source.replace(old,new)
    # pandas 3 removed float(single-element Series); select exactly one median
    # with Series.item(), which fails if the reviewed table is not scalar.
    for parameter in ('rRNA_SSU','rRNA_LSU','AMPHORA2'):
        old="float(medians.loc[medians['parameter']=='"+parameter+"',args.enrichment_preset])"
        new="float(medians.loc[medians['parameter']=='"+parameter+"',args.enrichment_preset].item())"
        if source.count(old)!=1:raise ValueError('Pinned ViromeQC median expression differs from reviewed scalar adapter')
        source=source.replace(old,new)
    old="'-k','0','--quiet']"
    if source.count(old)!=1:raise ValueError('Pinned ViromeQC DIAMOND command differs from reviewed tempdir adapter')
    source=source.replace(old,"'-k','0','--quiet','--tmpdir',args.tempdir]")
    runtime.write_text(source)
    dump(local/'adapter_manifest.json',{'upstream_script_sha256':sha(script),'runtime_script_sha256':sha(runtime),
        'adapter':'viromeqc_1_0_2_exit_guards_scalar_medians_tmpdir_v3','database_path':str(db),
        'copied_files':{str(f.relative_to(local)):sha(f) for f in local.rglob('*') if f.is_file() and 'index' not in f.parts},
        'database_read_only':True,'implicit_database_download_allowed':False})
    return runtime

def qualification(u,counts,virome_status):
    meta=u['metadata'];reasons=[]
    for field in ('mda','enrichment','negative_controls'):
        if not meta.get(field) or meta[field].lower() in ('unknown','na'):reasons.append(field+'_unknown')
    if meta.get('mda','').lower() in ('yes','true'):reasons.append('MDA_distorts_quantitative_abundance')
    if meta.get('run_accession') and meta.get('public_library_source','unknown').lower()=='unknown':reasons.append('public_library_source_unknown')
    if not counts['fragments']:q={'status':'ineligible','reason':'zero_clean_fragments'}
    else:q={'status':'provisional' if reasons else 'eligible','reason':';'.join(reasons) or 'documented_preprocessing'}
    if u['material_type']=='vlp' and virome_status!='assessed':q={'status':'not_assessed','reason':'VLP_enrichment_not_assessed'}
    return {'quantification':q,'host_adjustment':{'status':q['status'] if u['material_type']=='bulk_metagenome' and counts['fragments'] else 'ineligible','reason':q['reason']+';bulk_only;host_reference_required'},'activity':{'status':'not_assessed','reason':'activity_deferred'}}
def reference_evidence(unit,human_index,virome_db,reference_manifest,software_manifest):
    if not software_manifest:raise ValueError('R1 requires a pinned software manifest')
    tools=json.loads(Path(software_manifest).read_text())
    qc=tools.get('reads-qc',{})
    if qc.get('runtime')=='native-locked-prefix':
        import os
        if qc.get('platform')!='linux-aarch64' or qc.get('container') is not None or qc.get('container_executed') is not False:
            raise ValueError('Native R1 software manifest must explicitly distinguish native execution from a container')
        if qc.get('native_environment')!='reads-qc' or os.environ.get('CRC_SPARK_ENVIRONMENT')!='reads-qc':
            raise ValueError('Native R1 environment identity differs from task activation')
        if qc.get('lock_sha256')!=os.environ.get('CRC_SPARK_ENV_LOCK_SHA256'):
            raise ValueError('Native R1 software lock differs from activated environment')
        if qc.get('native_runtime_fingerprint')!=os.environ.get('CRC_SPARK_ENV_FINGERPRINT'):
            raise ValueError('Native R1 full software identity differs from activated environment')
        if not re.fullmatch(r'[a-f0-9]{64}',qc.get('native_manifest_sha256','')) or not re.fullmatch(r'[a-f0-9]{64}',qc.get('write_guard_sha256','')):
            raise ValueError('Native R1 manifest lacks fixed manifest/guard identity')
    elif not re.fullmatch(r'sha256:[a-f0-9]{64}',qc.get('container','')) and '@sha256:' not in qc.get('container',''):
        raise ValueError('R1 software container must be fixed by digest')
    if not qc.get('lock_sha256'):raise ValueError('R1 software lock hash is missing')
    enabled=[]
    if unit['metadata'].get('human_remove','true').lower()=='true':enabled.append('human')
    if unit['material_type']=='vlp':enabled.append('viromeqc')
    if enabled and not reference_manifest:raise ValueError('Enabled reference branch requires --reference-manifest')
    inventory=json.loads(Path(reference_manifest).read_text()) if reference_manifest else {}
    evidence={}
    for name in enabled:
        item=inventory.get(name,{})
        if not item.get('files') or not item.get('version'):raise ValueError('Reference unavailable/unindexed: '+name)
        base=Path(item['path']).resolve()
        if name=='human' and str(Path(human_index or '').resolve())!=str(Path(item.get('index_prefix','')).resolve()):
            raise ValueError('Human index prefix differs from pinned reference manifest')
        if name=='viromeqc' and (not virome_db or Path(virome_db).resolve()!=base):
            raise ValueError('ViromeQC DB differs from pinned reference manifest')
        for entry in item['files']:
            file=Path(entry['path']);file=file if file.is_absolute() else base/file
            if not file.resolve().is_relative_to(base) or not file.is_file() or not file.stat().st_size:
                raise ValueError('Reference inventory path missing/unsafe: '+str(file))
            if not re.fullmatch('[a-f0-9]{64}',entry.get('sha256','')):raise ValueError('Reference SHA256 missing')
        evidence[name]=item
    return {'software':qc,'references':evidence,
            'software_manifest_sha256':sha(software_manifest),
            'reference_manifest_sha256':sha(reference_manifest) if reference_manifest else None,
            'reference_integrity_contract':'full content hashes verified once by launcher preflight; per-unit checks inventory and paths',
            'adapter_sha256':sha(__file__)}

def human_command(profile,threads,index,sam,r1,r2=None):
    k=int(profile['reported_alignments'])
    if k<1:raise ValueError('reported_alignments must be positive')
    if profile['min_mapq']>0 and k>1:
        raise ValueError('Positive human MAPQ threshold requires reported_alignments=1 (best-hit mode); paired -k can report unknown MAPQ 255')
    cmd=['bowtie2','--very-sensitive-local','-p',str(threads),'-x',str(index),'-S',str(sam)]
    if k>1:cmd+=['-k',str(k)]
    return cmd+(['-1',str(r1),'-2',str(r2)] if r2 else ['-U',str(r1)])

def named_fastqc(reads,directory,prefix,threads,log):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    links=directory/'_inputs';links.mkdir(exist_ok=True);paths=[]
    try:
        for mate,source in enumerate(reads,1):
            suffix='.fastq.gz' if str(source).endswith('.gz') else '.fastq'
            target=links/(prefix+'_R'+str(mate)+suffix)
            if target.exists() or target.is_symlink():
                if not target.is_symlink():raise ValueError('QC alias unexpectedly contains a real file')
                target.unlink()
            target.symlink_to(Path(source).resolve());paths.append(target)
        # FastQC's Java image writer ignores TMPDIR unless --dir supplies
        # java.io.tmpdir. Keep only this invocation's cache in the task output.
        import tempfile
        with tempfile.TemporaryDirectory(prefix='_fastqc_tmp_',dir=directory.resolve()) as temp:
            run(['fastqc','--threads',str(min(threads,2)),'--dir',temp,'--outdir',directory]+paths,log)
    finally:
        for target in paths:
            if target.is_symlink():target.unlink()
        if links.exists() and not any(links.iterdir()):links.rmdir()

def virome_parameters(profile):
    supplied=profile.get('vlp',{})
    preset=supplied.get('enrichment_preset','human')
    if preset not in ('human','environmental'):raise ValueError('Invalid ViromeQC enrichment_preset')
    length=supplied.get('minlen',75);quality=supplied.get('minqual',20)
    if isinstance(length,bool) or not isinstance(length,int) or length<1:raise ValueError('Invalid ViromeQC minlen')
    if isinstance(quality,bool) or not isinstance(quality,int) or not 0<=quality<=93:raise ValueError('Invalid ViromeQC minqual')
    return {'enrichment_preset':preset,'minlen':length,'minqual':quality,'diagnostic_filter_only':True}

def process(unit_file,profile_file,out,threads=2,human_index=None,virome_db=None,reference_manifest=None,software_manifest=None):
    u=json.load(open(unit_file));p=json.load(open(profile_file));out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    if p.get('quality_encoding','phred33')!='phred33' or u.get('quality_encoding','phred33')!='phred33':raise ValueError('This R1 profile requires declared phred33 encoding')
    evidence=reference_evidence(u,human_index,virome_db,reference_manifest,software_manifest)
    vp=virome_parameters(p) if u['material_type']=='vlp' else None
    evidence['viromeqc_parameters']=vp
    dump(out/'preprocessing_evidence.json',evidence)
    for f,h in {**u['checksums'],**u.get('source_singleton_checksums',{})}.items():
        if sha(f)!=h:raise ValueError('Input checksum changed since R0')
    raw=out/'fastqc_raw';raw.mkdir(exist_ok=True)
    named_fastqc([u[k] for k in ('reads_1','reads_2') if u[k]],raw,u['unit_id']+'_raw',threads,out/'commands.log')
    a=out/'trimmed_1.fastq.gz';b=out/'trimmed_2.fastq.gz' if u['layout']=='PE' else None
    cmd=['fastp','--thread',str(threads),'--in1',u['reads_1'],'--out1',a,'--json',out/'fastp.json','--html',out/'fastp.html','--length_required',str(p['min_length']),'--qualified_quality_phred',str(p['quality_phred']),'--unqualified_percent_limit',str(p['unqualified_percent_limit'])]
    singleton=[]
    if b:
        cmd+=['--in2',u['reads_2'],'--out2',b,'--detect_adapter_for_pe','--unpaired1',out/'singleton_1.fastq.gz','--unpaired2',out/'singleton_2.fastq.gz']
        singleton=[str(out/f'singleton_{i}.fastq.gz') for i in (1,2)]
    for mate in (1,2):
        if u['metadata'].get(f'adapter_r{mate}'):cmd += ['--adapter_sequence' if mate==1 else '--adapter_sequence_r2',u['metadata'][f'adapter_r{mate}']]
    run(cmd,out/'commands.log');trimmed_counts=fq_counts(a,b)
    # Acquisition orphans are separately cleaned and never added to PE quantification.
    for i,source in enumerate(u.get('source_singletons',[])):
        target=out/f'acquisition_singleton_{i}.fastq.gz'
        scmd=['fastp','--thread',str(threads),'--in1',source,'--out1',target,'--json',out/f'acquisition_singleton_{i}.json','--html',out/f'acquisition_singleton_{i}.html','--length_required',str(p['min_length']),'--qualified_quality_phred',str(p['quality_phred']),'--unqualified_percent_limit',str(p['unqualified_percent_limit'])]
        if fq_counts(source)['fragments']:run(scmd,out/'commands.log');singleton.append(str(target))
    dehuman=u['metadata'].get('human_remove','true').lower()
    if dehuman not in ('true','false'):raise ValueError('human_remove must be true/false')
    human_stats={'status':'not_assessed','reason':'explicitly_disabled'}
    if dehuman=='true':
        if not human_index:raise ValueError('Human removal enabled but no frozen human index')
        sam=out/'human.sam';args=human_command(p['human'],threads,human_index,sam,a,b)
        if trimmed_counts['fragments']:
            run(args,out/'human_alignment.log');a,b,human_stats=remove_human(a,b,sam,out/'dehuman',p['human']);human_stats['status']='assessed'
        else:
            a,b=str(a),str(b) if b else None
            human_stats={'status':'not_assessed','reason':'zero_fragments_after_fastp','removed_fragments':0,'retained_fragments':0}
        cleaned=[]
        for i,s in enumerate(singleton):
            if not Path(s).exists() or not fq_counts(s)['fragments']:continue
            ss=out/f'singleton_{i}.sam';run(human_command(p['human'],threads,human_index,ss,s),out/'human_alignment.log')
            sa,_,st=remove_human(s,None,ss,out/f'dehuman_singleton_{i}',p['human']);cleaned.append(sa)
        singleton=cleaned
    else:
        a,b=str(a),str(b) if b else None
        if not u['metadata'].get('human_removal_reason'):raise ValueError('Disabled human removal needs an explicit reason')
    counts=fq_counts(a,b);singleton=[s for s in singleton if Path(s).exists()];single_counts=[fq_counts(s) for s in singleton]
    cleanqc=out/'fastqc_clean';cleanqc.mkdir(exist_ok=True)
    if counts['fragments']:named_fastqc([a]+([b] if b else []),cleanqc,u['unit_id']+'_clean',threads,out/'commands.log')
    vstatus='not_applicable'
    if u['material_type']=='vlp':
        if not virome_db:raise ValueError('VLP route requires pinned ViromeQC marker database')
        if counts['fragments']:
            # Bioconda's ViromeQC resolves index relative to its script. Copy the tiny code tree,
            # link the immutable externally mounted DB, and leave its quality filtering diagnostic-only.
            executable=shutil.which('viromeQC.py')
            if not executable:raise ValueError('ViromeQC executable unavailable')
            script=Path(executable).resolve();local=out/'viromeqc_tool'
            runtime=prepare_viromeqc(script,local,virome_db)
            run(['python3',runtime,'-i',a]+([b] if b else [])+['-o',out/'viromeqc.tsv','--debug','--enrichment_preset',vp['enrichment_preset'],'--minlen',str(vp['minlen']),'--minqual',str(vp['minqual']),'--bowtie2_threads',str(threads),'--diamond_threads',str(threads),'--tempdir',out,'--sample_name',u['unit_id']],out/'viromeqc.log');vstatus='assessed'
            # Keep the referenced DB out of recursively published result trees.
            if (local/'index').is_symlink():(local/'index').unlink()
        else:vstatus='not_assessed'
    preprocessing_payload={'profile_sha256':sha(profile_file),'software':evidence['software'],'references':evidence['references'],'adapter_sha256':evidence['adapter_sha256'],'material_type':u['material_type'],'human_remove':dehuman,'adapters':{k:v for k,v in u['metadata'].items() if k.startswith('adapter_r')}}
    preprocessing_id='prep_'+hashlib.sha256(json.dumps(preprocessing_payload,sort_keys=True).encode()).hexdigest()[:20]
    result={k:u[k] for k in ('unit_id','biological_sample_id','library_id','subject_id','material_type','molecule','platform','layout','source_readsets','metadata')}
    result.update({'schema_version':1,'quality_encoding':'phred33','quality_encoding_basis':'explicit_stage_contract_not_inferred','preprocessing_id':preprocessing_id,'profile_sha256':sha(profile_file),'preprocessing_evidence':evidence,'reads_1':a,'reads_2':b,'singletons':singleton,'counts':counts,'singleton_counts':single_counts,'checksums':{f:sha(f) for f in [a,b]+singleton if f},'raw_counts':u['counts'],'trimmed_counts':trimmed_counts,'source_singleton_checksums':u.get('source_singleton_checksums',{}),'raw_checksums':u['checksums'],'human_filter':human_stats,'viromeqc_status':vstatus,'eligibility':qualification(u,counts,vstatus),'status':'success' if counts['fragments'] else 'zero_clean_reads'})
    dump(out/'clean_unit.json',result);return result
def collect(dirs,out,profile,raw_manifest):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True);units=[];summary=[]
    for d in dirs:
        d=Path(d);u=json.load(open(d/'clean_unit.json'));target=out/'units'/u['unit_id'];target.parent.mkdir(exist_ok=True)
        if d.resolve()!=target.resolve():shutil.copytree(d,target,dirs_exist_ok=True,symlinks=True,ignore=shutil.ignore_patterns('*.sam','*.sqlite','*.sqlite-*','trimmed_*.fastq.gz','singleton_*.fastq.gz'))
        for field in ('reads_1','reads_2'):
            if u[field]:u[field]=str(target/Path(u[field]).relative_to(d.resolve()))
        u['singletons']=[str(target/Path(f).relative_to(d.resolve())) for f in u['singletons']]
        u['checksums']={str(target/Path(f).relative_to(d.resolve())):h for f,h in u['checksums'].items()}
        dump(target/'clean_unit.json',u);units.append(u)
        summary.append({'unit_id':u['unit_id'],'status':u['status'],'raw_fragments':u['raw_counts']['fragments'],'clean_fragments':u['counts']['fragments'],'clean_reads':u['counts']['reads'],'singleton_reads':sum(s['reads'] for s in u['singleton_counts']),'material_type':u['material_type'],'qualification':u['eligibility']['quantification']['status'],'reasons':u['eligibility']['quantification']['reason']})
    collection_ids=sorted({u['preprocessing_id'] for u in units})
    if any(u['profile_sha256']!=sha(profile) for u in units):raise ValueError('Cannot collect units from a different preprocessing profile')
    if len({u['preprocessing_evidence']['software_manifest_sha256'] for u in units})>1:raise ValueError('Cannot collect units from different software manifests')
    collection_hash=hashlib.sha256(json.dumps(collection_ids).encode()).hexdigest()[:20]
    dump(out/'clean_reads_manifest.json',{'schema_version':1,'kind':'clean_reads_manifest','preprocessing_id':'prep_collection_'+collection_hash,'unit_preprocessing_ids':collection_ids,'profile_sha256':sha(profile),'raw_manifest_sha256':sha(raw_manifest),'units':units})
    write_tsv(out/'sample_status.tsv',summary,['unit_id','material_type','status','raw_fragments','clean_fragments','clean_reads','singleton_reads','qualification','reasons'])
    run(['multiqc',str(out/'units'),'--outdir',out/'multiqc','--force','--dirs','--dirs-depth','2'],out/'multiqc.log')
def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='cmd',required=True)
    a=s.add_parser('unit');a.add_argument('--unit',required=True);a.add_argument('--profile',required=True);a.add_argument('--outdir',required=True);a.add_argument('--threads',type=int,default=2);a.add_argument('--human-index');a.add_argument('--virome-db');a.add_argument('--reference-manifest');a.add_argument('--software-manifest')
    a=s.add_parser('collect');a.add_argument('--dirs',nargs='+',required=True);a.add_argument('--outdir',required=True);a.add_argument('--profile',required=True);a.add_argument('--raw-manifest',required=True)
    a=p.parse_args()
    if a.cmd=='unit':process(a.unit,a.profile,a.outdir,a.threads,a.human_index,a.virome_db,a.reference_manifest,a.software_manifest)
    else:collect(a.dirs,a.outdir,a.profile,a.raw_manifest)
if __name__=='__main__':main()
