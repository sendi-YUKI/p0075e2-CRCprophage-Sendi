#!/usr/bin/env python3
"""Run or resume A2+A3+A4, optionally first executing the established A0+A1 entry."""
import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'bin'))
from nextflow_resume import select_resume_session, verify_nextflow_session
import spark_runtime
from core_genomes import validate_reference_manifest


def sha(p):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([p])
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def save(p,obj):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')


def read(p):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([p])
    return json.loads(Path(p).read_text())


def rows(p):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([p])
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def run(argv,**kw):return subprocess.run(list(map(str,argv)),check=True,**kw)


def code_manifest():
    files=[]
    for folder in ['bin','modules','workflows','conf','assets']:
        files.extend(f for f in (REPO/folder).rglob('*') if f.is_file() and '__pycache__' not in f.parts)
    files.extend(REPO/f for f in ['main.nf','main_core.nf','nextflow.config','pixi.toml','pixi.lock'])
    return dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
                dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=REPO,text=True)),
                branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip(),
                baseline_v01=subprocess.check_output(['git','rev-parse','v0.1^{commit}'],cwd=REPO,text=True).strip(),
                files={str(f.relative_to(REPO)):sha(f) for f in sorted(files)})


def images_preflight(info):
    sources={str(p.relative_to(REPO)):read(p) for p in sorted((REPO/'assets').glob('core_*tools.json'))}
    pins=set()
    def walk(x):
        if isinstance(x,dict):
            for key,value in x.items():
                if key in ['container','image_id'] and isinstance(value,str) and ('@sha256:' in value or re.fullmatch(r'sha256:[a-f0-9]{64}',value)):pins.add(value)
                walk(value)
        elif isinstance(x,list):
            for value in x:walk(value)
    walk(sources);pins.add(read(REPO/'assets/containers.json')['python'])
    inspected=[]
    if spark_runtime.enabled():
        inspected=[spark_runtime.preflight_image(image) for image in sorted(pins)]
        spark_runtime.native_record(info,[row['environment'] for row in inspected])
        save(info/'container_manifest.json',dict(runtime='native-locked-prefix',legacy_tool_sources=sources,images=inspected))
        return inspected
    for image in sorted(pins):
        data=json.loads(subprocess.check_output(['docker','image','inspect',image],text=True))[0]
        if data['Architecture']!='amd64' or data['Os']!='linux':raise ValueError('Expected linux/amd64 '+image)
        run(['docker','run','--rm','--platform','linux/amd64','--cpus','1','--memory','1g','--user',f'{os.getuid()}:{os.getgid()}',image,'sh','-c','ps -e -o pid= -o ppid= >/dev/null && command -v awk && command -v sed && command -v grep'],stdout=subprocess.DEVNULL)
        inspected.append(dict(reference=image,image_id=data['Id'],repo_digests=data.get('RepoDigests',[]),platform='linux/amd64'))
    save(info/'container_manifest.json',dict(images=inspected,tool_sources=sources))
    return inspected


def database_preflight(args,info):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([args.checkm2_db,args.functions_db,args.reference_dir])
    assets=read(REPO/'assets/core_genomes_tools.json')['checkm2']
    if args.checkm2_db.stat().st_size!=assets['database_bytes'] or sha(args.checkm2_db)!=assets['database_sha256']:raise ValueError('CheckM2 database content differs from the pinned full database')
    refs=validate_reference_manifest(read(args.reference_dir/'reference_manifest.json'))
    for ref in refs['references']:
        target=args.reference_dir/ref['fasta']
        if not target.resolve().is_relative_to(args.reference_dir.resolve()):raise ValueError('Reference escapes its catalog directory')
        if sha(target)!=ref['fasta_sha256']:raise ValueError('Reference checksum mismatch '+str(target))
    # Full consumed model inventory is produced by the verified asset installer.
    inventory=args.functions_db/'consumed_files_manifest.json'
    if not inventory.exists():raise ValueError('Functions consumed-file inventory missing; run bootstrap/core_functions_assets.py')
    data=read(inventory); files=data.get('files',[]) if isinstance(data,dict) else data
    if not files:raise ValueError('Empty functions model inventory')
    if isinstance(files,dict):files=[dict(path=k,**v) for k,v in files.items()]
    for item in files:
        f=args.functions_db/item['path']
        if not f.resolve().is_relative_to(args.functions_db.resolve()):raise ValueError('Unsafe model manifest path')
        if f.stat().st_size!=item['bytes'] or sha(f)!=item['sha256']:raise ValueError('Function model modified '+str(f))
    manifests={}
    checkm2_manifest=args.checkm2_db.parent.parent/('database_manifest.spark.json' if spark_runtime.enabled() else 'database_manifest.json')
    if spark_runtime.enabled() and not checkm2_manifest.is_file():raise ValueError('Native CheckM2 server database receipt missing: '+str(checkm2_manifest))
    for f in [args.functions_db/'database_manifest.json',inventory,args.reference_dir/'reference_manifest.json',checkm2_manifest,args.checkm2_db.parent.parent/('runtime_manifest.spark.json' if spark_runtime.enabled() else 'runtime_manifest.json')]:
        if f.exists():manifests[str(f)]=dict(sha256=sha(f),content=read(f))
    result=dict(checkm2_file=str(args.checkm2_db),checkm2_sha256=assets['database_sha256'],reference_catalog=refs,manifests=manifests)
    save(info/'database_manifest.json',result)
    return result


def configured_resource_limits():
    """Configured ceilings; observed Slurm allocation is recorded separately."""
    if spark_runtime.enabled():
        if os.environ.get('CRC_SPARK_RESOURCE_PROFILE')=='exclusive20':return dict(cpus=20,memory_gib=110,task_memory_gib=46,queue_size=4)
        return dict(cpus=6,memory_gib=48,task_memory_gib=46,queue_size=1)
    return dict(cpus=6,memory_gib=12,queue_size=2)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    inputs=parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--input',type=Path);inputs.add_argument('--a1-result',type=Path)
    parser.add_argument('--outdir',type=Path,required=True);parser.add_argument('--resume',action='store_true')
    parser.add_argument('--database-root',type=Path,default=Path.home()/'databases')
    parser.add_argument('--work-dir',type=Path,default=Path.home()/'scratch/crc-phage/core/work')
    parser.add_argument('--checkm2-db',type=Path,default=Path.home()/'databases/core/checkm2/1.1.0/CheckM2_database/uniref100.KO.1.dmnd')
    parser.add_argument('--conjscan-models-dir',type=Path)
    parser.add_argument('--element-pfam-db',type=Path)
    parser.add_argument('--functions-db',type=Path,default=Path.home()/'databases/core/functions')
    parser.add_argument('--reference-dir',type=Path,default=REPO/'assets/reference_catalog_empty')
    parser.add_argument('--metadata',type=Path,default=REPO/'assets/empty_metadata.tsv')
    parser.add_argument('--genbanks',type=Path,default=REPO/'assets/empty_genbanks.tsv')
    for option,default in [('completeness-min',90),('contamination-max',5),('high-completeness-min',95),('high-contamination-max',2),('fragmentation-n50',10000),('taxonomy-ani-min',95),('taxonomy-fragment-fraction-min',0.65),('votu-ani-percent',95),('votu-af-shorter-percent',85)]:parser.add_argument('--'+option,type=int if option=='fragmentation-n50' else float,default=default)
    parser.add_argument('--host-association-config',type=Path)
    parser.add_argument('--host-association-bundle',type=Path,default=REPO/'assets/host_source_empty')
    parser.add_argument('--catalog-policy',type=Path)
    parser.add_argument('--catalog-evidence',type=Path,default=REPO/'assets/catalog_evidence_empty')
    parser.add_argument('--enable-kofam',action='store_true')
    parser.add_argument('--host-recombination-config',type=Path)
    parser.add_argument('--whole-genome-search-config',type=Path)
    parser.add_argument('--kofam-db',type=Path,default=Path('/srv/CRC-PHIRE/databases/kofam/2026-08-02'))
    parser.add_argument('--host-config',type=Path)
    parser.add_argument('--host-mask-bundle',type=Path,default=REPO/'assets/host_mask_empty')
    args=parser.parse_args()
    for k,v in vars(args).items():
        if isinstance(v,Path):setattr(args,k,v.expanduser().resolve())
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([v for v in vars(args).values() if isinstance(v,Path)])
    from mainline_preflight import inspect_extensions
    try:mainline=inspect_extensions(args,require_frozen=True)
    except (ValueError,RuntimeError,KeyError,OSError) as exc:parser.error(str(exc))
    out=args.outdir;info=out/'pipeline_info';state_file=info/'run_status.json' 
    if args.a1_result and out==args.a1_result:parser.error('Core output must be separate from preserved A0+A1 output')
    if out.exists() and any(out.iterdir()) and not args.resume:parser.error('Output is not empty; choose a new directory or --resume')
    if args.resume and not state_file.exists():parser.error('No saved core run status in this output directory')
    out.mkdir(parents=True,exist_ok=True);info.mkdir(exist_ok=True)
    lock=(info/'run.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:parser.error('This output directory already has an active process')
    previous=read(state_file) if state_file.exists() else {}
    previous_workflow=read(info/'workflow_status.json') if (info/'workflow_status.json').exists() else {}
    resume_session=select_resume_session(previous,previous_workflow,REPO/'.nextflow/history') if args.resume else None
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    history=info/'history'/stamp;history.mkdir(parents=True)
    for f in list(info.iterdir()):
        if f.is_file() and f.name!='run.lock':shutil.copyfile(f,history/f.name)
    record=dict(current_target='A2_A3_A4',status='running',attempt=stamp,started_utc=stamp,nextflow_started=False,exit_code=None,
                outdir=str(out),actual_repo_path=str(REPO),resource_limits=configured_resource_limits(),
                resume_target=resume_session,resume_session_id=resume_session,
                parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()})
    save(state_file,record)
    try:
        if spark_runtime.enabled():
            record['resource_allocation']=spark_runtime.allocation()
            save(state_file,record)
        if not 50<=args.votu_ani_percent<=100 or not 1<args.votu_af_shorter_percent<=100:raise ValueError('vOTU thresholds use percent units, e.g. 95 and 85')
        for name in ['completeness_min','contamination_max','high_completeness_min','high_contamination_max']:
            if not 0<=getattr(args,name)<=100:raise ValueError('QC thresholds must be percent 0..100')
        if args.input:
            a1=out/'a0_a1'
            cmd=[REPO/'bin/run','--input',args.input,'--outdir',a1,'--database-root',args.database_root,'--work-dir',args.work_dir/'a0_a1']
            if args.resume and (a1/'pipeline_info').exists():cmd.append('--resume')
            run(cmd,cwd=REPO)
            args.a1_result=a1
        a1=args.a1_result
        success=read(a1/'pipeline_info/workflow_status.json')
        if success.get('success') is not True:raise ValueError('A0+A1 workflow did not succeed')
        inherited=read(a1/'aggregate_manifest.json')
        for name,expected in inherited['files'].items():
            if Path(name).name!=name or sha(a1/name)!=expected:raise ValueError('A0+A1 inherited output checksum mismatch: '+name)
        required=['prophage_master.tsv','candidates.fna','sample_status.tsv','checkv_quality.tsv','coordinate_mapping.tsv','contig_map.tsv','genome_stats.tsv','caller_hits.tsv']
        input_manifest=dict(a1_result=str(a1),files={name:sha(a1/name) for name in required},samples=rows(a1/'sample_status.tsv'),
                            metadata=dict(path=str(args.metadata),sha256=sha(args.metadata)),genbanks=dict(path=str(args.genbanks),sha256=sha(args.genbanks)))
        inherited_data=read(a1/'data_manifests.json')
        for row in input_manifest['samples']:
            prepared=a1/'samples'/row['genome_id']/'input/prepared'
            if spark_runtime.enabled():
                reject_ipc_storage([prepared,*prepared.rglob('*')])
            actual=sha(prepared/'genome.fna')
            evidence=[m for m in inherited_data if m['genome_id']==row['genome_id'] and m['manifest_type']=='input_manifest.json']
            if len(evidence)!=1 or sha(prepared/'input_manifest.json')!=evidence[0]['manifest_sha256'] or actual!=evidence[0]['content']['normalized_fasta_sha256']:
                raise ValueError('Prepared FASTA/inherited geNomad input evidence mismatch '+row['genome_id'])
            for m in inherited_data:
                if m['genome_id']==row['genome_id'] and 'source_normalized_fasta_sha256' in m['content'] and m['content']['source_normalized_fasta_sha256']!=actual:
                    raise ValueError('Stale caller evidence for '+row['genome_id'])
            input_manifest['files'][str(prepared/'genome.fna')]=actual
        sample_ids={r['genome_id'] for r in input_manifest['samples']}
        bank_rows=rows(args.genbanks)
        if len({r['genome_id'] for r in bank_rows})!=len(bank_rows):raise ValueError('Duplicate GenBank sample ID')
        input_manifest['genbanks']['files']=[]
        for r in bank_rows:
            if r['genome_id'] not in sample_ids:raise ValueError('Unknown GenBank genome ID '+r['genome_id'])
            if r['genbank']:
                gb=Path(r['genbank']);gb=gb if gb.is_absolute() else args.genbanks.parent/gb
                actual=sha(gb)
                expected=r.get('sha256',r.get('genbank_sha256',''))
                if expected and expected!=actual:raise ValueError('GenBank checksum mismatch '+str(gb))
                input_manifest['genbanks']['files'].append(dict(genome_id=r['genome_id'],path=str(gb),sha256=actual))
        save(info/'input_manifest.json',input_manifest)
        code=code_manifest();save(info/'code_manifest.json',code)
        images_preflight(info);database_preflight(args,info)
        save(info/'mainline_preflight.json',mainline)
        name='core_'+stamp.replace('Z','')+'_'+str(os.getpid())
        command=['nextflow','-log',info/('nextflow-'+stamp+'.log'),'run',REPO/'main_core.nf','-c',REPO/'conf/core.config','-name',name,'-work-dir',args.work_dir,'--a1_result',a1,'--outdir',out]
        for key in ['checkm2_db','functions_db','reference_dir','metadata','genbanks','completeness_min','contamination_max','high_completeness_min','high_contamination_max','fragmentation_n50','taxonomy_ani_min','taxonomy_fragment_fraction_min','votu_ani_percent','votu_af_shorter_percent']:
            command.extend(['--'+key,str(getattr(args,key))])
        if spark_runtime.enabled():
            command.extend(spark_runtime.nextflow_options())
        for optional in ('conjscan_models_dir','element_pfam_db'):
            value=getattr(args,optional)
            if value:
                if not value.exists() or any(c in str(value) for c in "'\n\r\0"):raise ValueError('Invalid element database path')
                command.extend(['--'+optional,str(value.resolve())])
        if args.host_association_config:command.extend(['--host_association_config',str(args.host_association_config),'--host_association_bundle',str(args.host_association_bundle)])
        if args.catalog_policy:command.extend(['--catalog_policy',str(args.catalog_policy),'--catalog_evidence',str(args.catalog_evidence)])
        if args.host_recombination_config:command.extend(['--host_recombination_config',str(args.host_recombination_config.resolve())])
        if args.whole_genome_search_config:command.extend(['--whole_genome_search_config',str(args.whole_genome_search_config.resolve())])
        if args.enable_kofam:command.extend(['--enable_kofam','true','--kofam_db',str(args.kofam_db)])
        if args.host_config:command.extend(['--host_config',str(args.host_config),'--host_mask_bundle',str(args.host_mask_bundle)])
        if resume_session:
            command.extend(['-resume',resume_session])
        record.update(nextflow_run_name=name,nextflow_started=True,command=list(map(str,command)),code_commit=code['commit'],code_dirty=code['dirty'])
        save(state_file,record);save(info/'core_run_manifest.json',record)
        with (info/('console-'+stamp+'.log')).open('w') as log:
            result=subprocess.run(list(map(str,command)),cwd=REPO,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError('Nextflow failed exit '+str(result.returncode)+'; see '+str(info/('console-'+stamp+'.log')))
        workflow=read(info/'workflow_status.json')
        actual_session=verify_nextflow_session(workflow,name,resume_session)
        if workflow.get('success') is not True:raise RuntimeError('Missing successful current core workflow receipt')
        record.update(nextflow_session_id=actual_session,resume_session_id=actual_session)
        trace=rows(info/'trace.tsv')
        if not any('CORE_REPORT' in r.get('name','') and r.get('status') in ['COMPLETED','CACHED'] for r in trace):
            raise RuntimeError('This attempt did not produce a successful report task')
        final=read(out/'core_data_manifest.json')
        if final.get('status')!='completed':raise RuntimeError('Core aggregation missing or failed')
        from host_source import verify_artifact_recursive
        if args.host_recombination_config:
            verify_artifact_recursive(out/'host/recombination/recombination','recombination_manifest.json')
        if args.whole_genome_search_config:
            verify_artifact_recursive(out/'whole_genome_evidence/whole_genome','whole_genome_evidence_manifest.json')
        record.update(status='completed',exit_code=0,catalog_version=final['catalog_version'],output_sha256={f.name:sha(f) for f in sorted(out.iterdir()) if f.is_file()},completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        save(state_file,record);save(info/'core_run_manifest.json',record)
        print(json.dumps(dict(status='completed',outdir=str(out),catalog_version=final['catalog_version'],report=str(out/'report_core.html')),indent=2))
        return 0
    except Exception as exc:
        record.update(status='failed',exit_code=1,error=str(exc),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        save(state_file,record);save(info/'core_run_manifest.json',record)
        print(str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
