#!/usr/bin/env python3
"""Output-bound, UUID-only launchers for development reads/V1/B1/B2 workflows."""
import argparse,datetime,fcntl,json,math,os,re,shutil,subprocess,sys,uuid,shlex
from pathlib import Path
R=Path(__file__).resolve().parents[1]
from reads_io import dump,sha,prepare
from nextflow_resume import select_resume_session,history_session
import spark_runtime
import contextlib
def read(p):return json.loads(Path(p).read_text())
def breadth_depth_feasible(label, path, entry):
    """Reject a breadth/depth pair that random coverage can never satisfy.

    breadth = 1 - exp(-depth) under Poisson coverage, so breadth_min must sit
    strictly below the ceiling its own mean_depth_min implies. breadth 0.80
    needs 1.61x; 0.75 needs 1.39x; 0.55 needs 0.80x. Pairs are frozen together
    per tier (protocol section 9.6b) precisely so this cannot drift apart.
    """
    b, d = entry.get('breadth_min'), entry.get('mean_depth_min')
    pid = entry.get('profile_id', '?')
    if b is None or d is None:
        return []
    try:
        b, d = float(b), float(d)
    except (TypeError, ValueError):
        return ['%s profile %s has non-numeric breadth_min/mean_depth_min' % (label, pid)]
    if d <= 0:
        # depth 0 imposes no floor, so any breadth is reachable given enough reads
        return []
    ceiling = 1.0 - math.exp(-d)
    if b >= ceiling:
        need = float('inf') if b >= 1 else -math.log(1 - b)
        return ['%s profile %s asks breadth_min=%.3f at mean_depth_min=%.3gx, but %.3gx '
                'caps breadth at %.3f (breadth=1-exp(-depth)); this pair can never '
                'produce a detection. breadth %.3f needs depth >= %.2fx, or lower '
                'breadth_min below %.3f. Freeze breadth and depth as a pair '
                '(protocol 9.6b), file %s'
                % (label, pid, b, d, d, ceiling, b, need, ceiling, path)]
    return []


def research_profile_guard(measurement_profile, detection_profiles):
    """Reject development profiles in study runs.

    Rationale: cohort_detection_dev_v1 uses breadth_field=breadth_full with
    mean_depth_min=0 (dev_breadth_10 sits at breadth_min=0.1), and
    cohort_measurement_dev_v1 uses mapq_min=10. Together that reproduces the
    legacy rule the study protocol forbids in section 12: claiming detection
    from 10 percent breadth plus a single MAPQ value. Engineering runs may keep
    using them; study runs must not.
    """
    problems=[]
    for label,path in (('measurement',measurement_profile),('detection',detection_profiles)):
        if not path:continue
        doc=read(check_path(path))
        entries=doc.get('profiles') if isinstance(doc.get('profiles'),list) else [doc]
        eligible=[e for e in entries if e.get('research_eligible') is True]
        if not eligible:
            ids=', '.join(str(e.get('profile_id','?')) for e in entries)
            problems.append(f"{label} profile {path} has no research_eligible entry (found: {ids})")
            continue
        unset=[e.get('profile_id','?') for e in eligible
               if e.get('breadth_min') is None and e.get('mean_depth_min') is None
               and e.get('status')=='pilot_pending']
        if unset:
            problems.append(f"{label} profile entries {unset} are pilot_pending with null thresholds; "
                            "freeze them after ground-truth simulation before a study run")
        for e in eligible:
            problems.extend(breadth_depth_feasible(label,path,e))
    if problems:
        raise ValueError('research run rejected: '+'; '.join(problems))
def check_path(value):
    p=Path(value).expanduser().resolve()
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([p])
    if any(c in str(p) for c in "'\n\r"):raise ValueError('Paths containing quotes/newlines unsupported')
    return p
def preflight_image(image):
    if spark_runtime.enabled():return spark_runtime.preflight_image(image)
    if not re.fullmatch(r'(?:.+@)?sha256:[a-f0-9]{64}',image or ''):raise ValueError('Container must be pinned by digest')
    j=json.loads(subprocess.check_output(['docker','image','inspect',image],text=True))[0]
    if j['Os']!='linux' or j['Architecture']!='amd64':raise ValueError('Require linux/amd64 image')
    return {'image':image,'image_id':j['Id']}
def db_verify_symlinks(root, records):
    """An explicit symlinks field is the complete non-following link inventory."""
    if not isinstance(records,list):raise ValueError('Database symlinks must be a complete list')
    expected={}
    for item in records:
        if not isinstance(item,dict) or not isinstance(item.get('path'),str) or not isinstance(item.get('target'),str):
            raise ValueError('Invalid database symlink record')
        name,target=item['path'],item['target'];relative=Path(name)
        if not name or relative.is_absolute() or '..' in relative.parts or name!=relative.as_posix() or name=='.' or not target or '\x00' in name+target:
            raise ValueError('Unsafe database symlink path/target: '+name)
        if name in expected:raise ValueError('Duplicate database symlink declaration: '+name)
        path=root/relative
        if not path.is_symlink():raise ValueError('Declared database symlink missing or replaced: '+name)
        if os.readlink(path)!=target:raise ValueError('Database symlink target changed: '+name)
        try:resolved=path.resolve(strict=True)
        except (OSError,RuntimeError) as exc:raise ValueError('Database symlink is dangling or cyclic: '+name) from exc
        if not resolved.is_relative_to(root) or not resolved.exists():
            raise ValueError('Database symlink target escapes root or is missing: '+name)
        expected[name]=target
    actual={}
    def traversal_error(error):raise error
    for parent,dirs,files in os.walk(root,followlinks=False,onerror=traversal_error):
        for name in dirs+files:
            path=Path(parent)/name
            if path.is_symlink():actual[path.relative_to(root).as_posix()]=os.readlink(path)
    if actual!=expected:
        missing=sorted(set(expected)-set(actual));extra=sorted(set(actual)-set(expected))
        raise ValueError('Database symlink inventory differs; missing='+repr(missing)+'; extra='+repr(extra))


def db_verify(d):
    """Validate files, and exact symlinks when the optional complete list is present.

    Legacy manifests without symlinks retain their file-checksum contract; they
    do not make a link-inventory completeness claim. An explicit [] means no links.
    """
    if spark_runtime.enabled() or os.environ.get('CRC_SPARK_ENVIRONMENT'):
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([d.get('path'),d.get('inventory'),d.get('index_prefix')])
    if d.get('status')!='ready':raise ValueError('Required database not ready: '+json.dumps(d))
    files=d.get('files')
    if d.get('inventory'):
        inventory=Path(d['inventory']).expanduser().resolve()
        expected=d.get('inventory_sha256')
        if not isinstance(expected,str) or not re.fullmatch(r'[a-f0-9]{64}',expected):
            raise ValueError('External database inventory requires a SHA256')
        if not inventory.is_file():
            raise ValueError('External database inventory checksum mismatch: '+str(inventory))
        payload=inventory.read_bytes()
        if __import__('hashlib').sha256(payload).hexdigest()!=expected:
            raise ValueError('External database inventory checksum mismatch: '+str(inventory))
        files=json.loads(payload)
    if not isinstance(files,list) or not files:raise ValueError('Ready database lacks nonempty file checksum inventory')
    root=Path(d['path']).expanduser().resolve()
    if not root.is_dir():raise ValueError('Database root is missing: '+str(root))
    if 'symlinks' in d:db_verify_symlinks(root,d['symlinks'])
    for item in files:
        if not isinstance(item,dict) or not isinstance(item.get('path'),str) or not item['path']:
            raise ValueError('Invalid database inventory file record')
        expected=item.get('sha256')
        if not isinstance(expected,str) or not re.fullmatch(r'[a-f0-9]{64}',expected):
            raise ValueError('Invalid database file SHA256: '+item['path'])
        path=(root/item['path']).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError('Database inventory path escapes root or is missing: '+str(path))
        if 'bytes' in item:
            size=item['bytes']
            if type(size) is not int or size<0 or path.stat().st_size!=size:
                raise ValueError('Database file bytes mismatch: '+str(path))
        if sha(path)!=expected:raise ValueError('Read-only database checksum mismatch: '+str(path))

def canonical_gene_inputs(values,out,previous=None):
    """Adapt CLI file arguments to V1's existing comma-tokenized parameter."""
    paths=[]
    for value in values or []:
        if ',' in value:raise ValueError('Canonical gene table paths must not contain commas; pass separate file arguments')
        path=check_path(value)
        if ',' in str(path):raise ValueError('Canonical gene table paths must not contain commas')
        if path==out or out in path.parents:raise ValueError('Canonical gene table inputs must not be inside the output directory')
        if not path.is_file():raise ValueError('Canonical gene table file missing: '+str(path))
        if path in paths:raise ValueError('Duplicate canonical gene table path: '+str(path))
        paths.append(path)
    entries=[dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size) for path in sorted(paths)]
    before={entry['path']:entry['sha256'] for entry in (previous or {}).get('inputs',[])}
    for entry in entries:
        if entry['path'] in before and before[entry['path']]!=entry['sha256']:
            raise ValueError('Canonical gene table changed at a previously recorded path; restore it or use a new outdir: '+entry['path'])
    return ','.join(str(path) for path in sorted(paths)) or None,dict(schema_version=1,kind='canonical_gene_tables',inputs=entries)

def stable_cohort_prepare(command, destination, info, attempt):
    """Preserve unchanged prepared-file bytes, mtimes and inodes across resume."""
    staging=info/'prepare_staging'/attempt
    staged=list(command)
    staged[staged.index('--outdir')+1]=str(staging)
    subprocess.run(staged,check=True)
    expected=set()
    for source in sorted(staging.rglob('*')):
        if not source.is_file():continue
        relative=source.relative_to(staging);expected.add(relative)
        target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists() or sha(target)!=sha(source):
            shutil.copy2(source,target)
    if destination.exists():
        for stale in sorted(destination.rglob('*')):
            if stale.is_file() and stale.relative_to(destination) not in expected:
                archived=info/'attempts'/attempt/'obsolete_prepared'/stale.relative_to(destination)
                archived.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(stale),str(archived))
    dump(info/'prepared_inventory.json',{str(k):sha(destination/k) for k in sorted(expected)})
    # Only files just generated in this private attempt are removed. No broad deletion.
    for source in sorted(staging.rglob('*'),key=lambda q:len(q.parts),reverse=True):
        if source.is_file():source.unlink()
        elif source.is_dir():source.rmdir()
    staging.rmdir()

def verify_viral_assets(dbroot):
    if spark_runtime.enabled():return spark_runtime.verify_viral_assets(dbroot)
    from run_pipeline import verify_databases,verify_index_receipt
    from core_functions import verify_databases as verify_functions
    catalog=read(R/'assets/required_assets.json')
    manifest=read(dbroot/'database_manifest.json')
    paths=verify_databases(dbroot,catalog,manifest)
    containers=read(R/'assets/containers.json')
    cv=paths['checkv']
    observed=subprocess.check_output(['docker','run','--rm','--platform','linux/amd64',containers['checkv'],'diamond','version'],text=True).strip()
    verify_index_receipt(cv/'genome_db/checkv_reps.dmnd',dbroot/'checkv_generated_index.json',
                         observed,containers['checkv'],cv/'genome_db/checkv_reps.faa')
    functions=dbroot/'core/functions'
    function_hash=verify_functions(functions)
    return {'genomad_db':str(paths['genomad']),'checkv_db':str(paths['checkv']),'functions_db':str(functions)}, {
        'a1_manifest_sha256':sha(dbroot/'database_manifest.json'),
        'checkv_index_manifest_sha256':sha(dbroot/'checkv_generated_index.json'),
        'functions_consumed_manifest_sha256':function_hash}

def groovy_literal(value):
    return "'" + value.replace(chr(92),chr(92)*2).replace("'",chr(92)+"'") + "'"

def readonly_mount_options(mounts):
    args=['--platform','linux/amd64','--user',str(os.getuid())+':'+str(os.getgid())]
    for mount in sorted(mounts):
        if ':' in mount or '\n' in mount or '\r' in mount:raise ValueError('Unsafe container bind path')
        args+=['-v',mount+':'+mount+':ro']
    return ' '.join(shlex.quote(value) for value in args)


def verify_terminal_outputs(workflow,out,evidence_tools=()):
    """A zero engine exit is insufficient when an empty channel skipped aggregation."""
    paths={
        'reads':['clean/clean_reads_manifest.json','clean/sample_status.tsv'],
        'assembly':['assemblies/assembly_manifest.json','assemblies/viral_assembly_manifest.json','assemblies/genomes.tsv','assemblies/sample_status.tsv'],
        'viral':['catalog/catalog_manifest.json','catalog/viral_sequence_master.tsv','catalog/source_occurrences.tsv','catalog/vOTU_members.tsv','catalog/vOTU_representatives.tsv','catalog/vOTU_representatives.fna','catalog/occurrence_votu_membership.tsv','genes/gene_calling_status.json','crosswalk/gene_set_crosswalk.tsv','crosswalk/crosswalk_manifest.json','functions/viral_function_status.json','summary/summary_manifest.json','summary/votu_function_summary.tsv'],
        'cohort':['summary/cohort_manifest.json','summary/measurement_long.tsv','summary/detection_long.tsv','summary/unit_metadata.tsv'],
        'association':['association/effects.tsv','association/analysis_sample_audit.tsv','association/design_report.tsv','association/analysis_spec.json','association/input_manifest.json','association/association_status.json','association/association_report.html','association/output_manifest.json'],
    }[workflow]
    for relative in paths:
        if not (out/relative).is_file():raise RuntimeError('Missing mandatory '+workflow+' terminal output: '+relative)
    if workflow=='viral':
        from viral_evidence_summary import verify_receipt
        summary_path=out/'evidence_summary/summary_manifest.json'
        summary=read(summary_path)
        if summary.get('status')!='completed' or set(summary.get('enabled_tools',[]))!=set(evidence_tools):
            raise RuntimeError('Viral evidence summary/request mismatch')
        if summary.get('catalog_manifest_sha256')!=sha(out/'catalog/catalog_manifest.json'):
            raise RuntimeError('Evidence summary uses another catalog')
        if not summary.get('output_files'):raise RuntimeError('Evidence summary has no output hashes')
        paths.append('evidence_summary/summary_manifest.json')
        for relative,checksum in summary['output_files'].items():
            base=summary_path.parent.resolve();path=(base/relative).resolve()
            if not path.is_relative_to(base) or not path.is_file() or sha(path)!=checksum:
                raise RuntimeError('Evidence summary output checksum mismatch')
            paths.append('evidence_summary/'+relative)
        for tool in evidence_tools:
            directory=out/'evidence'/tool/'evidence'
            verify_receipt(directory,tool,sha(out/'catalog/catalog_manifest.json'))
            if summary.get('evidence_receipt_sha256',{}).get(tool)!=sha(directory/'status.json'):
                raise RuntimeError('Evidence summary/tool receipt mismatch')
            paths.append('evidence/'+tool+'/evidence/status.json')
        for name in ('catalog/catalog_manifest.json','genes/gene_calling_status.json',
                     'crosswalk/crosswalk_manifest.json','functions/viral_function_status.json','summary/summary_manifest.json'):
            receipt=read(out/name)
            if receipt.get('status')!='completed':raise RuntimeError('V1 terminal status is not completed: '+name)
            base=(out/name).parent.resolve()
            files=receipt.get('output_files',{})
            if name.startswith('crosswalk/'):
                files={'gene_set_crosswalk.tsv':receipt.get('output_sha256')}
            if not files:raise RuntimeError('V1 terminal receipt has no output hashes: '+name)
            for relative,checksum in files.items():
                path=(base/relative).resolve()
                if not path.is_relative_to(base) or not path.is_file() or sha(path)!=checksum:
                    raise RuntimeError('V1 terminal output checksum mismatch: '+name+':'+relative)
    if workflow=='association':
        if read(out/'association/association_status.json').get('execution_status')!='success':
            raise RuntimeError('B2 did not report successful execution')
        files=read(out/'association/output_manifest.json').get('files',{})
        if not files:raise RuntimeError('B2 output manifest has no checksums')
        base=(out/'association').resolve()
        for relative,checksum in files.items():
            path=(base/relative).resolve()
            if not path.is_relative_to(base) or not path.is_file() or sha(path)!=checksum:
                raise RuntimeError('B2 output checksum mismatch: '+relative)
    return {relative:sha(out/relative) for relative in paths}

def code_manifest():
    files=[p for root in ('bin','modules','workflows','conf','assets') for p in (R/root).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    files += list(R.glob('main*.nf'))+[R/'nextflow.config']
    origin=read(R/'analysis1_source_manifest.json')
    return {'commit':None,'source_commit':origin['source_commit'],'snapshot':True,'source_manifest_sha256':sha(R/'analysis1_source_manifest.json'),'files':{str(p.relative_to(R)):sha(p) for p in files}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('workflow',choices=['reads','assembly','viral','cohort','association']);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');p.add_argument('--preview',action='store_true')
    p.add_argument('--samplesheet');p.add_argument('--merge-map');p.add_argument('--clean-reads-manifest');p.add_argument('--assembly-manifest');p.add_argument('--catalog-dir');p.add_argument('--catalog-scope',choices=['host_resolved_prophage','virome_discovery','combined_exploratory']);p.add_argument('--legacy-catalog');p.add_argument('--legacy-source-manifest')
    p.add_argument('--research-reference-config');p.add_argument('--enable-viral-kofam',action='store_true');p.add_argument('--research-catalog');p.add_argument('--research-source-manifest')
    p.add_argument('--viral-evidence-tools',default='',help='Opt-in V1 evidence: comma-separated bacphlip,vibrant,vcontact3; default disabled')
    p.add_argument('--canonical-gene-tables',nargs='+',help='Viral workflow: separate canonical A2 gene_table.tsv files; commas are forbidden')
    p.add_argument('--preprocessing-profile',default=str(R/'assets/profiles/reads_preprocessing_dev_v1.json'));p.add_argument('--measurement-profile',default=str(R/'assets/profiles/cohort_measurement_dev_v1.json'));p.add_argument('--detection-profiles',default=str(R/'assets/profiles/cohort_detection_dev_v1.json'));p.add_argument('--host-references');p.add_argument('--decoy-fasta');p.add_argument('--enable-phispy',action='store_true');p.add_argument('--research-run',action='store_true',help='Study run: reject profiles without research_eligible=true (protocol section 12)')
    p.add_argument('--assembly-threads',type=int,default=16,choices=range(1,21));p.add_argument('--assembly-memory-gib',type=int,default=40,choices=range(1,101));p.add_argument('--enable-metaphlan',action='store_true');p.add_argument('--association-long');p.add_argument('--association-metadata');p.add_argument('--association-spec');a=p.parse_args()
    if a.enable_viral_kofam and a.workflow!='viral':raise ValueError('--enable-viral-kofam requires viral mode')
    if a.viral_evidence_tools and a.workflow!='viral':raise ValueError('--viral-evidence-tools is only valid for viral')
    if a.canonical_gene_tables and a.workflow!='viral':raise ValueError('--canonical-gene-tables is only valid for viral')
    if (a.research_catalog or a.research_source_manifest) and a.workflow!='viral':raise ValueError('Research bridge flags require viral mode')
    if bool(a.research_catalog)!=bool(a.research_source_manifest):raise ValueError('Research bridge requires both catalogue and source manifest')
    if getattr(a,'research_run',False):research_profile_guard(getattr(a,'measurement_profile',None),getattr(a,'detection_profiles',None))
    out=check_path(a.outdir);dbroot=Path.home()/'databases';workroot=Path('/srv/CRC-PHIRE/analysis1/work/reads_v1')
    if out.is_relative_to(dbroot) or out.is_relative_to(workroot):raise ValueError('Outputs must be separate from database and Nextflow work directories')
    for key,value in vars(a).items():
        if key in ('workflow','outdir','catalog_scope') or not value:continue
        values=value if isinstance(value,list) else [value]
        for item in values:
            if not isinstance(item,str):continue
            inp=check_path(item)
            if inp==out or out in inp.parents:raise ValueError('Inputs must not be inside this new output directory')
    unbound_payload=[q for q in out.iterdir() if q.name!='.run.lock' and not (q.name=='pipeline_info' and q.is_dir() and not any(q.iterdir()))] if out.exists() else []
    if unbound_payload and not (out/'pipeline_info/stage_state.json').exists():
        raise ValueError('Existing nonempty directory has no matching stage state; choose a new output to protect historical results')
    out.mkdir(parents=True,exist_ok=True);info=out/'pipeline_info';info.mkdir(exist_ok=True)
    statefile=info/'stage_state.json';old=read(statefile) if statefile.exists() else {}
    if old and old.get('workflow')!=a.workflow:raise ValueError('Output directory bound to another workflow')
    if a.resume and not old:raise ValueError('Resume requires an existing output-bound stage state')
    if old and not a.resume:raise ValueError('Output has prior stage state; use --resume or new outdir')
    lock=open(out/'.run.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    history=R/'.nextflow/history';resume=select_resume_session(old,{},history) if a.resume else None
    name='readsbuild_'+a.workflow+'_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:6]
    params={'outdir':str(out),'assembly_threads':a.assembly_threads,'assembly_memory_gib':a.assembly_memory_gib};mounts={str(R),str(dbroot)};tools=read(R/'assets/reads_tools.json')
    database_manifest=spark_runtime.reads_database_manifest(info) if spark_runtime.enabled() else R/'assets/reads_databases.json'
    dbs=read(database_manifest) if database_manifest.exists() else {}
    if spark_runtime.enabled():params['database_manifest']=str(database_manifest)
    state={'schema_version':1,'target':'READS_B1_B2_BUILD','workflow':a.workflow,'status':'preflight','nextflow_started':False,'nextflow_run_name':name,'resume_session_id':resume,'args':vars(a),'real_data_validated':False,'scientifically_calibrated':False}
    if old:dump(info/'attempts'/(name+'.previous.json'),old)
    dump(statefile,state)
    try:
        images=[]
        used_database_evidence={}
        if a.research_catalog:
            from research_bridge import verify_catalog
            bridge_source=check_path(a.research_catalog);verify_catalog(bridge_source)
            if spark_runtime.enabled():spark_runtime.validate_output_inputs(out,[bridge_source])
        if a.workflow=='reads':
            if not a.samplesheet:raise ValueError('reads requires --samplesheet')
            if spark_runtime.enabled():
                from spark_acquire import validate_declared_reads_inputs, validate_prepared_reads
                validate_declared_reads_inputs(a.samplesheet,out)
            raw=prepare(a.samplesheet,out/'r0',a.merge_map);params.update(raw_reads_manifest=str(out/'r0/raw_reads_manifest.json'),preprocessing_profile=str(check_path(a.preprocessing_profile)))
            if spark_runtime.enabled():validate_prepared_reads(a.samplesheet,out,raw,a.merge_map)
            images.append(tools['reads-qc']['container'])
            if any(str(u['metadata'].get('human_remove','true')).lower()=='true' for u in raw['units']):db_verify(dbs.get('human',{}));params['human_index']=dbs['human']['index_prefix']
            if any(u['material_type']=='vlp' for u in raw['units']):db_verify(dbs.get('viromeqc',{}));params['virome_db']=dbs['viromeqc']['path']
            for u in raw['units']:
                for f in u['checksums']:mounts.add(str(Path(f).parent))
            config='reads'
        elif a.workflow=='assembly':
            if not a.clean_reads_manifest:raise ValueError('assembly requires --clean-reads-manifest')
            params['clean_reads_manifest']=str(check_path(a.clean_reads_manifest));params['enable_phispy']=a.enable_phispy;images.append(tools['reads-assembly']['container']);config='reads'
            for u in read(a.clean_reads_manifest)['units']:
                if spark_runtime.enabled():spark_runtime.validate_output_inputs(out,u['checksums'])
                for f,h in u['checksums'].items():
                    if sha(f)!=h:raise ValueError('Clean reads changed')
                    mounts.add(str(Path(f).parent))
            if a.enable_phispy:
                annotation_units=[u for u in read(a.clean_reads_manifest)['units'] if u['material_type']=='bacterial_isolate']
                if annotation_units:images.append(tools['reads-bakta']['container'])
                missing_gbff=[u for u in annotation_units if not u.get('metadata',{}).get('genbank')]
                if missing_gbff:
                    db_verify(dbs.get('bakta',{}));params['bakta_db']=dbs['bakta']['path']
                for u in annotation_units:
                    provided=u.get('metadata',{}).get('genbank')
                    if provided:
                        provided=check_path(provided)
                        if not provided.is_file():raise ValueError('Provided GBFF missing: '+str(provided))
                        mounts.add(str(provided.parent))
                        # Exact FASTA/GBFF sequence agreement is checked after assembly by the annotation adapter.

        elif a.workflow=='viral':
            if a.enable_viral_kofam:
                if not spark_runtime.enabled():raise ValueError('Viral KOfam requires Spark native runtime')
                spark_runtime.environment('kofam')
                deps=read(R/'assets/spark_databases.json')
                kd=deps['databases']['kofam']
                if kd.get('status')!='ready':raise ValueError('KOfam database not ready')
                kp=Path(kd['path'])
                params.update(enable_viral_kofam=True,kofam_db=str(kp),kofam_db_manifest_sha256=sha(kp/'database_manifest.json'))
            evidence_tools=a.viral_evidence_tools.split(',') if a.viral_evidence_tools else []
            if len(evidence_tools)!=len(set(evidence_tools)) or any(t not in ('bacphlip','vibrant','vcontact3') for t in evidence_tools):raise ValueError('Invalid/duplicate evidence tool')
            if evidence_tools:
                if not spark_runtime.enabled():raise ValueError('Viral evidence currently requires Spark native runtime')
                from viral_evidence import verify_assets
                evidence_contract=read(R/'assets/viral_evidence_tools.json')
                for tool in evidence_tools:
                    item=evidence_contract['tools'][tool]
                    if item['state']!='prepared':raise ValueError('Evidence tool not prepared: '+tool)
                    spark_runtime.environment(tool)
                    verify_assets(item)
            params['viral_evidence_tools']=','.join(evidence_tools)

            if not a.assembly_manifest or not a.catalog_scope:raise ValueError('viral requires --assembly-manifest and --catalog-scope')
            params.update(assembly_manifest=str(check_path(a.assembly_manifest)),catalog_scope=a.catalog_scope,viral_gene_container=tools['viral-gene']['container']);images.append(tools['viral-gene']['container']);config='viral'
            for k in ('legacy_catalog','legacy_source_manifest','research_catalog','research_source_manifest'):
                if getattr(a,k):params[k]=str(check_path(getattr(a,k)))
            if spark_runtime.enabled() and a.legacy_catalog:
                from spark_fs_guard import reject_ipc_storage
                legacy=check_path(a.legacy_catalog)
                reject_ipc_storage([legacy,*legacy.rglob('*')])
            canonical_receipt=info/'canonical_gene_tables_manifest.json'
            previous=read(canonical_receipt) if a.resume and canonical_receipt.is_file() else None
            canonical_value,canonical_manifest=canonical_gene_inputs(a.canonical_gene_tables,out,previous)
            params['canonical_gene_tables']=canonical_value
            dump(canonical_receipt,canonical_manifest)
            for entry in canonical_manifest['inputs']:mounts.add(str(Path(entry['path']).parent))
            verified,used_database_evidence=verify_viral_assets(dbroot);params.update(verified)
            images += list(read(R/'assets/containers.json').values())
            images += re.findall(r"container\s*=\s*'([^']+)'",(R/'conf/viral.config').read_text())
            for manifest_path in [a.assembly_manifest,a.legacy_source_manifest,a.research_source_manifest]:
                if not manifest_path:continue
                manifest_path=check_path(manifest_path)
                for u in read(manifest_path)['units']:
                    if u.get('assembly_fasta'):
                        path=Path(u['assembly_fasta'])
                        path=path if path.is_absolute() else manifest_path.parent/path
                        if spark_runtime.enabled():spark_runtime.validate_output_inputs(out,[path])
                        if sha(path)!=u.get('sha256'):raise ValueError('Assembly SHA256 changed: '+str(path))
                        mounts.add(str(path.resolve().parent))
                    if u.get('contig_map'):
                        path=Path(u['contig_map']);path=path if path.is_absolute() else manifest_path.parent/path
                        if spark_runtime.enabled():spark_runtime.validate_output_inputs(out,[path])
                        mounts.add(str(path.resolve().parent))
        elif a.workflow=='cohort':
            if not a.clean_reads_manifest or not a.catalog_dir or not a.catalog_scope:raise ValueError('cohort requires --clean-reads-manifest --catalog-dir --catalog-scope')
            cmd=[sys.executable,str(R/'bin/cohort_cli.py'),'prepare','--clean-reads-manifest',a.clean_reads_manifest,'--catalog-dir',a.catalog_dir,'--catalog-scope',a.catalog_scope,'--measurement-profile',a.measurement_profile,'--detection-profiles',a.detection_profiles,'--outdir',str(out/'prepared')]
            for k in ('host_references','decoy_fasta','research_reference_config'):
                if getattr(a,k):cmd+=['--'+k.replace('_','-'),getattr(a,k)]
            stable_cohort_prepare(cmd,out/'prepared',info,name)
            params['cohort_prepared']=str(out/'prepared');config='cohort';images.append(tools['reads-qc']['container'])
            params['cohort_enable_metaphlan']=a.enable_metaphlan
            if a.enable_metaphlan:
                if dbs.get('metaphlan',{}).get('status')!='ready':
                    raise ValueError('MetaPhlAn enabled but unavailable: '+json.dumps(dbs.get('metaphlan',{})))
                db_verify(dbs['metaphlan'])
                params.update(cohort_metaphlan_db=dbs['metaphlan']['path'],cohort_metaphlan_index=dbs['metaphlan']['index'],
                              cohort_metaphlan_database_manifest=dbs['metaphlan']['manifest'])
                images.append(tools['reads-metaphlan']['container'])
            for u in read(a.clean_reads_manifest)['units']:
                if spark_runtime.enabled():spark_runtime.validate_output_inputs(out,u['checksums'])
                for f in u['checksums']:mounts.add(str(Path(f).parent))
            mounts.add(str(check_path(a.catalog_dir)))
        else:
            for k in ('association_long','association_metadata','association_spec'):
                if not getattr(a,k):raise ValueError('association requires --'+k.replace('_','-'))
                params[k]=str(check_path(getattr(a,k)))
            config='association';assoc=read(R/'assets/association_tools.json');images.append(assoc.get('container',assoc.get('image_id')))
        if spark_runtime.enabled() and a.workflow=='reads':
            params['software_manifest']=str(spark_runtime.write_reads_software_manifest(info/'reads_tools.native.json'))
        dump(info/'container_manifest.json',[preflight_image(i) for i in sorted(set(images))]);dump(info/'code_manifest.json',code_manifest());dump(info/'database_manifest.json',{'reads_databases':dbs,'verified_workflow_assets':used_database_evidence});dump(info/'parameters.json',params)
        cfg=info/'mounts.config'
        opts='' if spark_runtime.enabled() else readonly_mount_options(mounts)
        content='// Spark native uses its explicit filesystem guard; no Docker mounts.\n' if spark_runtime.enabled() else 'docker.writableInputMounts = false\ndocker.runOptions = '+groovy_literal(opts)+'\n'
        if not cfg.exists() or cfg.read_text()!=content:cfg.write_text(content)
        cmd=['nextflow','run',str(R/('main_'+a.workflow+'.nf')),'-c',str(R/'conf'/f'{config}.config'),'-c',str(cfg),'-params-file',str(info/'parameters.json'),'-name',name,'-work-dir',str(workroot/a.workflow)]
        if spark_runtime.enabled():
            cmd+=spark_runtime.nextflow_options()
            evidence_envs=params.get('viral_evidence_tools','').split(',') if params.get('viral_evidence_tools') else []
            spark_runtime.native_record(info,sorted(set(['python']+[spark_runtime.image_environment(i) for i in images]+evidence_envs)))
            state['resource_limits']={'allocation':spark_runtime.allocation(),'profile':os.environ.get('CRC_SPARK_RESOURCE_PROFILE','standard'),'assembly_threads':a.assembly_threads,'assembly_memory_gib':a.assembly_memory_gib}
        if resume:cmd+=['-resume',resume]
        if a.preview:cmd+=['-preview']
        state.update(status='running',nextflow_started=True,command=cmd);dump(statefile,state)
        with (contextlib.nullcontext() if spark_runtime.enabled() else open(workroot/'core-heavy.lock','w')) as heavy:
            if heavy is not None:fcntl.flock(heavy,fcntl.LOCK_EX)
            else:spark_runtime.allocation()
            with open(info/(name+'.log'),'w') as log:rc=subprocess.run(cmd,cwd=R,stdout=log,stderr=subprocess.STDOUT).returncode
        actual=history_session(history,name)
        if rc==0 and not a.preview and not actual:raise RuntimeError('Successful Nextflow command lacks exact session UUID evidence')
        if resume and actual and actual!=resume:raise RuntimeError('Resume UUID mismatch')
        state.update(exit_code=rc,nextflow_session_id=actual,status='wiring_checked' if a.preview and rc==0 else 'success' if rc==0 else 'failed');dump(statefile,state)
        if rc==0 and not a.preview:
            terminal=verify_terminal_outputs(a.workflow,out,params.get('viral_evidence_tools','').split(',') if params.get('viral_evidence_tools') else [])
            if a.research_catalog:
                from research_bridge import verify_catalog
                bridge_dir=out/'research_bridge/catalog';verify_catalog(bridge_dir)
                terminal['research_bridge/catalog/catalog_manifest.json']=sha(bridge_dir/'catalog_manifest.json')
            if a.workflow=='viral' and params.get('enable_viral_kofam'):
                from host_source import verify_artifact
                kp=out/'viral_kofam/kofam';verify_artifact(kp,'viral_ko_manifest.json')
                terminal['viral_kofam/kofam/viral_ko_manifest.json']=sha(kp/'viral_ko_manifest.json')
            dump(info/'terminal_output_manifest.json',terminal)
        if rc:raise RuntimeError(f'Nextflow failed exit {rc}; see {info/name}.log')
        print(json.dumps({'workflow':a.workflow,'status':state['status'],'outdir':str(out),'nextflow_session_id':actual},indent=2))
    except Exception as e:
        state.update(status='failed',error=str(e));dump(statefile,state);raise
if __name__=='__main__':
    try:main()
    except Exception as e:print('ERROR: '+str(e),file=sys.stderr);sys.exit(1)
