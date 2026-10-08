#!/usr/bin/env python3
"""Preflight, execute, and preserve provenance for a real A0+A1 run."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from pathlib import PurePosixPath
import re
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'bin'))
from nextflow_resume import select_resume_session, verify_nextflow_session
import spark_runtime
SAMPLE_PROCESSES = {'PREPARE_GENOME':'prepare', 'GENOMAD_ENDTOEND':'genomad', 'NORMALIZE_CALLS':'normalize', 'CHECKV_ENDTOEND':'checkv', 'ANNOTATE_QUALITY':'quality', 'EMPTY_QUALITY':'quality'}
SUCCESS = {'COMPLETED', 'CACHED'}
MODEL_FILES = ('decision_forest.ubj', 'nn_classifier.h5', 'provirus_tagger.crfsuite', 'rbs_categories.tsv', 'score_calibration_weights.npz')

def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def read_table(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise RuntimeError(f'Malformed table header: {path}')
        rows=list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise RuntimeError(f'Malformed table row: {path}')
        return rows

def safe_relative(root, relative):
    parts=PurePosixPath(relative)
    if parts.is_absolute() or '..' in parts.parts or '\\' in relative or ':' in relative:
        raise RuntimeError(f'Unsafe relative manifest path: {relative}')
    result=root.joinpath(*parts.parts).resolve()
    if not result.is_relative_to(root.resolve()):
        raise RuntimeError(f'Manifest path escapes database root: {relative}')
    return result

def verify_databases(db_root, catalog, manifest):
    """Verify fetch_assets' actual receipt schema, independent of Docker."""
    if manifest.get('kind') != 'databases' or manifest.get('download_exit_code') != 0:
        raise RuntimeError('Database download manifest is not a successfully completed database fetch')
    databases={}
    for asset in catalog['databases']:
        base=safe_relative(db_root,asset['relative_directory'])
        directory=safe_relative(base,asset['archive_top_directory'])
        record=manifest.get('assets',{}).get(asset['id'],{})
        expected_status='downloaded_index_pending' if asset['requires_generated_index'] else 'ready'
        if record.get('version') != asset['version'] or record.get('status') != expected_status or record.get('url') != asset['url'] or not record.get('files'):
            raise RuntimeError(f'Missing successful verified database receipt: {asset["id"]}')
        extraction_path=base/'extraction_manifest.json'
        if not extraction_path.is_file():
            raise RuntimeError(f'Missing extraction receipt: {extraction_path}')
        extraction=read_json(extraction_path)
        if record['files'] != extraction.get('files') or record.get('archive',{}).get('sha256') != extraction.get('archive_sha256'):
            raise RuntimeError(f'Database/extraction receipt disagreement: {asset["id"]}')
        archive=record.get('archive',{})
        if not re.fullmatch(r'[0-9a-f]{64}',archive.get('sha256','')) or archive.get('bytes') != asset['download_bytes']:
            raise RuntimeError(f'Invalid archive checksum/size receipt: {asset["id"]}')
        for algorithm in ('md5','sha256'):
            if asset.get('expected_'+algorithm) and archive.get(algorithm) != asset['expected_'+algorithm]:
                raise RuntimeError(f'Archive receipt differs from publisher {algorithm}: {asset["id"]}')
        for item in record['files']:
            file=safe_relative(base,item['path'])
            if not file.is_file() or file.stat().st_size != item['bytes'] or sha(file) != item['sha256']:
                raise RuntimeError(f'Database content changed or missing: {file}')
        for rel in asset['required_paths']:
            required=safe_relative(directory,rel)
            if not (required.is_file() or required.is_dir()):
                raise RuntimeError(f'Missing database path: {directory/rel}')
        if asset['id']=='genomad' and (directory/'version.txt').read_text(encoding='utf-8').strip() != asset['version']:
            raise RuntimeError('geNomad database version.txt differs from its pinned version')
        databases[asset['id']]=directory
    return databases

def verify_index_receipt(index, receipt, version, container, source):
    if not receipt.is_file():
        raise RuntimeError('CheckV index exists without provenance; inspect it before reuse.')
    previous=read_json(receipt)
    expected={'sha256':sha(index),'bytes':index.stat().st_size,'diamond_version':version,'container':container,'source_sha256':sha(source)}
    if any(previous.get(key) != value for key,value in expected.items()):
        raise RuntimeError('CheckV index provenance/content/source differs from pinned runtime')

def resume_target(previous, workflow_status, history_path=None):
    """Return a full output-scoped UUID; never pass a run name to Nextflow."""
    return select_resume_session(previous, workflow_status, history_path or REPO/'.nextflow/history')

def trace_sample_task(record):
    match=re.search(r'(?:^|:)([A-Z_]+)\s+\(([^()]*)\)\s*$',record.get('name',''))
    if not match or match.group(1) not in SAMPLE_PROCESSES:
        return None
    process, sample=match.groups()
    if process=='EMPTY_QUALITY' and sample.endswith(':zero_candidates'):
        sample=sample.removesuffix(':zero_candidates')
    return process,sample

def execution_rows(genome_ids, trace_rows, status_by_sample, *, nextflow_started, trace_available):
    """Use this attempt's trace; missing/failed calls never imply zero viruses."""
    grouped={genome_id:[] for genome_id in genome_ids}
    for record in trace_rows:
        parsed=trace_sample_task(record)
        if parsed:
            process,sample=parsed
            if sample not in grouped:
                raise RuntimeError(f'Trace refers to sample absent from current inputs: {sample}')
            grouped[sample].append((process,record))
    results=[]
    for genome_id, tasks in grouped.items():
        row=dict(genome_id=genome_id,execution_status='not_run',genome_qc_status='not_assessed',candidate_count='unknown',failed_tasks='',task_exit_codes='',evidence='no_task_in_current_attempt_trace')
        row.update({stage+'_status':'not_run' for stage in SAMPLE_PROCESSES.values()})
        failures=[]
        for process,record in tasks:
            raw_status=record.get('status','UNKNOWN').upper()
            stage=SAMPLE_PROCESSES[process]
            exit_code=record.get('exit','-')
            failed=raw_status in {'FAILED','ABORTED','KILLED'} or (exit_code not in ('0','-','',None) and raw_status!='CACHED')
            row[stage+'_status']='failed' if failed else 'cached' if raw_status=='CACHED' else 'completed' if raw_status=='COMPLETED' else 'unfinished'
            if failed:
                failures.append(process)
            if process=='EMPTY_QUALITY' and raw_status in SUCCESS and not failed:
                row['checkv_status']='skipped_empty'
        row['failed_tasks']=';'.join(sorted(set(failures)))
        row['task_exit_codes']=';'.join(f'{process}={record.get("exit","unknown")}' for process,record in tasks)
        if tasks:
            row['evidence']='current_attempt_nextflow_trace'
            row['execution_status']='failed' if failures else 'completed' if row['quality_status'] in ('completed','cached') else 'incomplete'
        elif nextflow_started and not trace_available:
            row['execution_status']='unknown'
            row['evidence']='nextflow_started_but_trace_unavailable'
        # Only a completed normalization in this attempt (or an explicitly
        # cached one) can establish candidate counts. Ignore stale final files.
        published=status_by_sample.get(genome_id)
        if row['normalize_status'] in ('completed','cached') and published and published.get('genomad_status')=='completed':
            count=published.get('candidate_count','unknown')
            if not str(count).isdigit():
                raise RuntimeError(f'Invalid published candidate count for {genome_id}: {count}')
            row['candidate_count']=count
            if count=='0' and row['execution_status']=='completed':
                row['execution_status']='completed_zero_candidates'
        results.append(row)
    return results

def save_execution_status(info, out, inputs, trace_path, nextflow_started):
    available=trace_path.is_file()
    traces=read_table(trace_path) if available else []
    published={}
    for sample in inputs:
        sample_id=sample['genome_id']
        successful={parsed[0] for trace in traces if (parsed:=trace_sample_task(trace)) and parsed[1]==sample_id and trace.get('status','').upper() in SUCCESS and trace.get('exit','-') in ('0','-','')}
        for stage in ('final','normalized'):
            if stage=='final' and not successful.intersection({'ANNOTATE_QUALITY','EMPTY_QUALITY'}):
                continue
            if stage=='normalized' and 'NORMALIZE_CALLS' not in successful:
                continue
            candidates=list((out/'samples'/sample_id/stage).rglob('sample_status.tsv'))
            if len(candidates)>1:
                raise RuntimeError(f'Multiple published {stage} sample status files for {sample_id}')
            if candidates:
                rows=read_table(candidates[0])
                if len(rows)==1 and rows[0].get('genome_id')==sample_id:
                    published[sample_id]=rows[0]
                    break
    rows=execution_rows([r['genome_id'] for r in inputs],traces,published,nextflow_started=nextflow_started,trace_available=available)
    fields=['genome_id','execution_status','prepare_status','genomad_status','normalize_status','checkv_status','quality_status','candidate_count','genome_qc_status','failed_tasks','task_exit_codes','evidence']
    with (info/'sample_execution_status.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    return rows

def validate_final_outputs(out, genome_ids):
    expected=['prophage_master.tsv','caller_hits.tsv','checkv_quality.tsv','coordinate_mapping.tsv','contig_map.tsv','genome_stats.tsv','sample_status.tsv','stage_counts.tsv','report.html','candidates.fna','data_manifests.json','aggregate_manifest.json']
    for name in expected:
        if not (out/name).is_file():
            raise RuntimeError(f'Missing final output: {name}')
    statuses=read_table(out/'sample_status.tsv')
    if len(statuses)!=len(genome_ids) or {r.get('genome_id') for r in statuses} != set(genome_ids):
        raise RuntimeError('Final sample-status IDs differ from current input samples')
    if any(r.get('sample_status') not in {'completed_with_candidates','zero_candidates'} or r.get('genome_qc_status')!='not_assessed' for r in statuses):
        raise RuntimeError('Final sample statuses are incomplete or claim unassessed bacterial QC passed')
    master=read_table(out/'prophage_master.tsv')
    counts={sample:0 for sample in genome_ids}
    ids=set()
    for row in master:
        sample=row.get('genome_id')
        if sample not in counts or not row.get('candidate_id') or row['candidate_id'] in ids:
            raise RuntimeError('Final master contains an unknown sample or duplicate candidate ID')
        counts[sample]+=1;ids.add(row['candidate_id'])
    for status in statuses:
        if str(counts[status['genome_id']])!=status.get('candidate_count'):
            raise RuntimeError('Final master and sample status candidate counts differ')
    for name in ('caller_hits.tsv','checkv_quality.tsv','coordinate_mapping.tsv'):
        rows=read_table(out/name)
        if len(rows)!=len(ids) or {r.get('candidate_id') for r in rows}!=ids:
            raise RuntimeError(f'Final candidate ID sets differ in {name}')
    manifest=read_json(out/'aggregate_manifest.json')
    if manifest.get('samples')!=sorted(genome_ids):
        raise RuntimeError('Aggregate manifest belongs to different input samples')
    for name,digest in manifest.get('files',{}).items():
        path=safe_relative(out,name)
        if not path.is_file() or sha(path)!=digest:
            raise RuntimeError(f'Final output changed after aggregation: {name}')
    return {name:sha(out/name) for name in expected}

def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()

def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    tmp.replace(path)

def capture(command, *, check=True):
    p=subprocess.run([str(x) for x in command],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if check and p.returncode:
        raise RuntimeError(f'Command failed ({p.returncode}): {command}\n{p.stdout}')
    return p.stdout.strip()

def stream(command):
    print('+ '+' '.join(map(str,command)),flush=True)
    subprocess.run([str(x) for x in command],cwd=REPO,check=True)

def docker_run(image, *command, mounts=(), memory='12g', cpus=6):
    cmd=['docker','run','--rm','--platform','linux/amd64','--cpus',str(cpus),'--memory',memory,
         '--user',f'{os.getuid()}:{os.getgid()}']
    for host,guest,mode in mounts:
        cmd+=['--mount',f'type=bind,src={host},dst={guest}'+(',readonly' if mode=='ro' else '')]
    return cmd+[image,*command]

def prepare_runtime(db_root, info):
    if spark_runtime.enabled():
        return spark_runtime.prepare_a0_runtime(db_root, info)
    if sys.platform!='linux':
        raise RuntimeError('Run this entrypoint inside the initialized WSL Ubuntu project.')
    if capture(['docker','info','--format','{{.OSType}}'])!='linux':
        raise RuntimeError('Docker is not serving Linux containers; enable Ubuntu WSL integration.')
    nf=capture(['nextflow','-version'])
    if '25.10.4' not in nf:
        raise RuntimeError('Expected Nextflow25.10.4; source bootstrap/environment.sh')
    containers=read_json(REPO/'assets/containers.json')
    images={}
    for name,image in containers.items():
        if '@sha256:' not in image: raise RuntimeError(f'Unpinned container: {name}')
        try:
            inspection=json.loads(capture(['docker','image','inspect',image]))[0]
        except RuntimeError:
            stream(['docker','pull','--platform','linux/amd64',image])
            inspection=json.loads(capture(['docker','image','inspect',image]))[0]
        if inspection.get('Architecture')!='amd64' or inspection.get('Os')!='linux':
            raise RuntimeError(f'Unexpected image platform: {name}')
        images[name]={'requested':image,'id':inspection['Id'],'repo_digests':inspection.get('RepoDigests',[]),'architecture':inspection['Architecture'],'os':inspection['Os']}
        # Nextflow task metrics require ps even for short validation processes.
        images[name]['nextflow_task_tools'] = capture(docker_run(
            image, 'bash', '-euc',
            'command -v bash; command -v ps; command -v awk; command -v sed; command -v grep; ps -e -o pid= -o ppid='))
    save(info/'container_manifest.json',images)
    catalog=read_json(REPO/'assets/required_assets.json')
    manifest_path=db_root/'database_manifest.json'
    if not manifest_path.is_file(): raise RuntimeError(f'Missing required database manifest: {manifest_path}')
    manifest=read_json(manifest_path)
    databases=verify_databases(db_root,catalog,manifest)
    db=databases['checkv']; index=db/'genome_db/checkv_reps.dmnd'
    # A lock prevents simultaneous writes when two run wrappers start together.
    import fcntl
    with (db_root/'.checkv-index.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        receipt=db_root/'checkv_generated_index.json'
        pending_receipt=db_root/'checkv_generated_index.pending.json'
        source=db/'genome_db/checkv_reps.faa'
        version=capture(docker_run(containers['checkv'],'diamond','version'))
        if index.exists():
            # Recover the narrow interruption window between the atomic index
            # rename and receipt rename, without trusting an unrecorded index.
            candidate_receipt=receipt if receipt.is_file() else pending_receipt
            verify_index_receipt(index,candidate_receipt,version,containers['checkv'],source)
            if candidate_receipt==pending_receipt:
                pending_receipt.replace(receipt)
        else:
            temporary_index=db/'genome_db/.checkv_reps.building.dmnd'
            if temporary_index.is_symlink():
                raise RuntimeError('Refusing a symlink at the temporary DIAMOND index path')
            temporary_index.unlink(missing_ok=True)
            stream(docker_run(containers['checkv'],'diamond','makedb','--in','/db/genome_db/checkv_reps.faa','--db','/db/genome_db/.checkv_reps.building','--threads','6',mounts=[(db,'/db','rw')]))
            if not temporary_index.is_file(): raise RuntimeError('DIAMOND did not create CheckV index')
            save(pending_receipt,{'sha256':sha(temporary_index),'bytes':temporary_index.stat().st_size,'diamond_version':version,'container':containers['checkv'],'source_sha256':sha(source)})
            temporary_index.replace(index)
            pending_receipt.replace(receipt)
    shutil.copy2(receipt,info/receipt.name)
    shutil.copy2(manifest_path,info/'database_manifest.json')
    versions={'nextflow':nf,'java':capture(['java','-version']),'docker':capture(['docker','version']),
              'genomad':capture(docker_run(containers['genomad'],'genomad','--version')),
              'checkv':capture(docker_run(containers['checkv'],'python','-c',"import importlib.metadata; print(importlib.metadata.version('checkv'))")),
              'python_validation':capture(docker_run(containers['python'],'python3','--version')),
              'pixi':capture(['pixi','--version']), 'nf_core':capture(['pixi','run','--locked','nf-core','--version']),
              'resource_budget':{'cpus':6,'memory_gb':12,'max_tasks':2}}
    save(info/'software_versions.json',versions)
    for name in ('genomad','checkv'):
        if not re.search(r'(?<![0-9.])'+re.escape(catalog['tools'][name]['version'])+r'(?![0-9.])',versions[name]):
            raise RuntimeError(f'Actual {name} version differs from pinned tool catalog: {versions[name]}')
    # Official v1.12.0 GenomadData defines package/data and these five assets.
    # Streaming hashing works even if the pinned container's Python is <3.11.
    model_code="""import hashlib,json
from genomad._paths import GenomadData
p=GenomadData.data_dir
required=%r
for name in required:
    if not (p/name).is_file(): raise RuntimeError('Missing bundled geNomad asset: '+name)
files=[]
for f in sorted(p.rglob('*')):
    if f.is_file():
        h=hashlib.sha256()
        with f.open('rb') as handle:
            for block in iter(lambda:handle.read(1048576),b''): h.update(block)
        files.append({'path':str(f.relative_to(p)),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
print(json.dumps(files))
""" % (MODEL_FILES,)
    models=json.loads(capture(docker_run(containers['genomad'],'python','-c',model_code)))
    if not models: raise RuntimeError('No bundled geNomad data/models found; inspect package layout before proceeding.')
    save(info/'genomad_model_manifest.json',{'container':containers['genomad'],'source':'https://github.com/apcamargo/genomad/blob/v1.12.0/genomad/_paths.py','data_directory':'GenomadData.data_dir','required_assets':MODEL_FILES,'files':models})
    return databases

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--database-root',type=Path,default=Path.home()/'databases')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--work-dir',type=Path,default=Path.home()/'scratch/crc-phage/work')
    args=p.parse_args(argv)
    out=args.outdir.expanduser().resolve(); sheet=args.input.expanduser().resolve()
    info=out/'pipeline_info'; info.mkdir(parents=True,exist_ok=True)
    previous_state=read_json(info/'run_status.json') if (info/'run_status.json').is_file() else {}
    previous_workflow=read_json(info/'workflow_status.json') if (info/'workflow_status.json').is_file() else {}
    target=resume_target(previous_state,previous_workflow) if args.resume else None
    history=info/'history'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    history.mkdir(parents=True)
    prior=history/'previous';prior.mkdir()
    for previous in info.iterdir():
        if previous.is_file(): shutil.copy2(previous,prior/previous.name)
    if (info/'input').is_dir():
        shutil.copytree(info/'input',prior/'input')
    trace_path=history/'trace.tsv'
    state={'status':'preflight','current_target':'A0_A1','input':str(sheet),'outdir':str(out),'resume_requested':args.resume,'resume_target':target,'resume_session_id':target,'nextflow_started':False,'started_utc':datetime.now(timezone.utc).isoformat(),'attempt_directory':str(history),'trace_path':str(trace_path)}
    inputs=[]
    result=1
    try:
        save(info/'run_status.json',state)
        if (out/'prophage_master.tsv').exists() and not args.resume:
            raise RuntimeError('This output already has results. Use --resume for this run or choose a new output directory.')
        if args.resume and not target and ((out/'prophage_master.tsv').exists() or (out/'samples').exists()):
            raise RuntimeError('Cannot identify this output directory\'s previous Nextflow session; do not use bare -resume. Inspect pipeline_info/history first.')
        # All previous aliases were archived above. Removing only these known
        # generated aliases prevents old evidence masquerading as this attempt.
        for name in ('trace.tsv','report.html','timeline.html','dag.html','workflow_status.json','output_manifest.json'):
            (info/name).unlink(missing_ok=True)
        stream([sys.executable,REPO/'bin/phageflow.py','validate-sheet','--input',sheet,'--outdir',info/'input'])
        resolved=info/'input/normalized_samples.tsv'
        rows=read_table(resolved)
        inputs=[dict(r,sha256=sha(r['fasta']),bytes=Path(r['fasta']).stat().st_size) for r in rows]
        if spark_runtime.enabled():spark_runtime.validate_output_inputs(out,[r['fasta'] for r in rows]+[sheet,args.database_root,args.work_dir])
        save(info/'input_manifest.json',{'source_samplesheet':str(sheet),'sheet_sha256':sha(sheet),'genomes':inputs})
        code={str(f.relative_to(REPO)):sha(f) for pattern in ('*.nf','*.config','bin/*.py','bootstrap/*.sh','assets/*.json','pixi.toml','pixi.lock') for f in REPO.glob(pattern) if f.is_file()}
        for directory in ['modules','workflows','conf']:
            code.update({str(f.relative_to(REPO)):sha(f) for f in (REPO/directory).rglob('*') if f.is_file()})
        git_commit=capture(['git','rev-parse','HEAD'],check=False)
        git_status=capture(['git','status','--porcelain=v1'],check=False)
        valid_commit=bool(re.fullmatch('[0-9a-f]{40,64}',git_commit))
        save(info/'code_manifest.json',{'git_commit':git_commit if valid_commit else None,'git_dirty':bool(git_status) if valid_commit else None,'git_status':git_status,'git_commit_error':None if valid_commit else git_commit,'files':code})
        state['preflight_phase']='containers_databases_models';save(info/'run_status.json',state)
        databases=prepare_runtime(args.database_root.expanduser().resolve(),info)
        run_name='crcphage_'+datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
        cmd=['nextflow','-log',str(history/'nextflow.log'),'run',str(REPO/'main.nf'),'-profile','docker','-name',run_name,'--input',str(resolved),'--outdir',str(out),
             '--genomad_db',str(databases['genomad']),'--checkv_db',str(databases['checkv']),'-work-dir',str(args.work_dir.expanduser().resolve()),
             '-with-trace',str(trace_path),'-with-report',str(history/'report.html'),'-with-timeline',str(history/'timeline.html'),'-with-dag',str(history/'dag.html')]
        if spark_runtime.enabled():
            cmd[cmd.index('-profile') + 1] = 'spark_native'
            cmd.extend(spark_runtime.nextflow_options())
        if target: cmd.extend(['-resume',target])
        state.update(status='running',command=cmd,nextflow_run_name=run_name,effective_resume=bool(target))
        save(info/'run_status.json',state)
        with (history/'console.log').open('w',encoding='utf-8') as log:
            proc=subprocess.Popen(cmd,cwd=REPO,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            state['nextflow_started']=True;save(info/'run_status.json',state)
            for line in proc.stdout:
                print(line,end='',flush=True);log.write(line);log.flush()
            rc=proc.wait()
        state.update(exit_code=rc,status='validating_outputs' if rc==0 else 'failed')
        if rc: raise RuntimeError(f'Nextflow exited {rc}; inspect trace and {history}/console.log')
        workflow=read_json(info/'workflow_status.json')
        actual_session=verify_nextflow_session(workflow,run_name,target)
        if workflow.get('success') is not True:
            raise RuntimeError('Missing successful current A0+A1 workflow receipt')
        state.update(nextflow_session_id=actual_session,resume_session_id=actual_session)
        save(info/'output_manifest.json',validate_final_outputs(out,[r['genome_id'] for r in inputs]))
        state['status']='completed'
        result=0
    except Exception as error:
        state.update(status='failed',error=str(error))
        print(f'ERROR: {error}',file=sys.stderr)
    finally:
        for name in ('trace.tsv','report.html','timeline.html','dag.html','nextflow.log'):
            if (history/name).is_file():
                shutil.copy2(history/name,info/name)
        try:
            execution=save_execution_status(info,out,inputs,trace_path,state['nextflow_started'])
            shutil.copy2(info/'sample_execution_status.tsv',out/'sample_execution_status.tsv')
            if result==0 and any(row['execution_status'] not in {'completed','completed_zero_candidates'} for row in execution):
                raise RuntimeError('Exit-zero Nextflow run lacks successful final tasks for every input in its current trace')
            state['sample_execution_status']=str(info/'sample_execution_status.tsv')
        except Exception as error:
            state.update(status='failed',execution_status_error=str(error))
            print(f'ERROR recording execution status: {error}',file=sys.stderr)
            result=1
        state['finished_utc']=datetime.now(timezone.utc).isoformat()
        save(info/'run_status.json',state)
        for artifact in info.iterdir():
            if artifact.is_file(): shutil.copy2(artifact,history/artifact.name)
    if state['status']=='completed':
        print(f'Real run completed. Report: {out}/report.html')
    return 0 if state['status']=='completed' else 1

if __name__=='__main__': raise SystemExit(main())
