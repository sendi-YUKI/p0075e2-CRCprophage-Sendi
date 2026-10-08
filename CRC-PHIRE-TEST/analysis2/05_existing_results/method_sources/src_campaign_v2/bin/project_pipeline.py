#!/usr/bin/env python3
"""Sequential project manifest over existing CRC entries; check never executes tools."""
import argparse,datetime,hashlib,html,json,os,re,shutil,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(os.environ.get('CRC_SPARK_REPO','/srv/CRC-PHIRE/projects/crc-pipeline')).resolve()
sys.path.insert(0,str(ROOT/'bin'))
KINDS={
 'aggregate_f1':('python','host_source.py',['aggregate-f1']),
 'detection_calibration':('python','detection_calibration.py',[]),
 'a0':('python','run',[]),'core':('python','run-core',[]),'core-from-reads':('python','run-core-from-reads',[]),
 'reads':('python','run-stage',['reads']),'assembly':('python','run-stage',['assembly']),
 'viral':('python','run-stage',['viral']),'cohort':('python','run-stage',['cohort']),
 'association_legacy':('python','run-stage',['association']),
 'research_bridge':('core-blast','research_bridge.py',[]),'whole_genome_evidence':('core-blast','whole_genome_evidence.py',[]),
 'research_summary':('python','core_research_summary.py',[]),
 'host_core_reference':('core-blast','host_core_reference.py',[]),'host_recombination':('recombination','host_recombination.py',['run']),
 'host_lineage':('recombination','host_lineage.py',[]),'host_source':('python','host_source.py',['run']),
 'community_modules':('python','community_modules.py',[]),'community_analysis':('python','community_analysis.py',['run'])
}
TOKEN=re.compile(r'\{\{([A-Za-z0-9_.-]+)\.outdir\}\}')
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):
 p=Path(p)
 if p.stat().st_size>4*2**20:raise ValueError('Project metadata file exceeds4MiB')
 def pairs(items):
  d={}
  for k,v in items:
   if k in d:raise ValueError('Duplicate JSON key '+k)
   d[k]=v
  return d
 return json.loads(p.read_text(),object_pairs_hook=pairs)
def save(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);temp=p.with_suffix(p.suffix+'.tmp');temp.write_text(json.dumps(x,indent=2)+'\n');temp.replace(p)
def resolve(x,outputs):
 if isinstance(x,str):
  def repl(m):
   if m[1] not in outputs:raise ValueError('Unavailable upstream output '+m[1])
   return str(outputs[m[1]])
  return TOKEN.sub(repl,x)
 if isinstance(x,list):return [resolve(v,outputs) for v in x]
 if isinstance(x,dict):return {k:resolve(v,outputs) for k,v in x.items()}
 return x
def inspect(path,outdir):
 cfg=read(path)
 if set(cfg)!={'schema_version','mode','authorization_manifest','stages'} or cfg['schema_version']!=1 or cfg['mode'] not in ('engineering_smoke','formal_research'):raise ValueError('Project manifest schema mismatch')
 seen=set();nodes=[]
 for n in cfg['stages']:
  if set(n)!={'id','kind','depends_on','parameters','config_templates'}:raise ValueError('Stage keys mismatch')
  if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,63}',n['id']) or n['id'] in seen or n['kind'] not in KINDS:raise ValueError('Stage identity/type invalid')
  if not set(n['depends_on'])<=seen:raise ValueError('Missing/cyclic/unsorted stage dependency')
  for key in n['parameters']:
   if not re.fullmatch('[a-z][a-z0-9-]*',key) or key in ('outdir','resume','work-dir'):raise ValueError('Uncontrolled parameter '+key)
  dependencies=set(TOKEN.findall(json.dumps(n['parameters'])))
  for item in n['config_templates']:
   if set(item)!={'parameter','template'} or not Path(item['template']).is_file():raise ValueError('Config template missing')
   if not re.fullmatch('[a-z][a-z0-9-]*',item['parameter']) or item['parameter'] in ('outdir','resume','work-dir'):raise ValueError('Uncontrolled template parameter')
   value=read(item['template']);dependencies.update(TOKEN.findall(json.dumps(value)))
  if not dependencies<=set(n['depends_on']):raise ValueError('Undeclared output dependency')
  seen.add(n['id']);nodes.append(dict(id=n['id'],kind=n['kind'],depends_on=n['depends_on'],outdir=str(Path(outdir)/n['id']/'attempt_1'),environment=KINDS[n['kind']][0]))
 if not nodes:raise ValueError('Empty project plan')
 return cfg,dict(status='structurally_valid',analysis_started=False,environment_probes=False,scientific_readiness='requires_each_stage_contract',stages=nodes)
def run(a):
 import spark_runtime
 cfg,plan=inspect(a.config,a.outdir);allocation=spark_runtime.allocation()
 dest=Path(a.outdir).resolve();spark_runtime.validate_output_inputs(dest,[a.config,cfg['authorization_manifest']])
 for root in [ROOT/'envs',Path('/srv/CRC-PHIRE/databases'),Path('/srv/CRC-PHIRE/scratch/crc-phage')]:
  if dest==root or dest.is_relative_to(root) or root.is_relative_to(dest):raise ValueError('Output overlaps protected root')
 auth=read(cfg['authorization_manifest'])
 if auth.get('purpose')!=cfg['mode']:raise ValueError('Authorization purpose mismatch')
 if cfg['mode']=='engineering_smoke':
  if auth.get('campaign')!='basic_build_20260930' or auth.get('status')!='registered' or auth.get('input_kind')!='synthetic':raise ValueError('Only registered synthetic campaign inputs accepted')
  if not auth.get('files'):raise ValueError('Register exact fixture files')
  for row in auth['files']:
   p=Path(row['path']).resolve()
   if not any(p.is_relative_to(ROOT/'results'/n) for n in ('basic_build_20260930','methods_prepare_20261001')) or not p.is_file() or sha(p)!=row['sha256']:raise ValueError('Engineering fixture identity mismatch')
  old_budget=ROOT/'results/basic_build_20260930/budget.json'
  continued=ROOT/'results/methods_prepare_20261001/budget.json'
  budget=read(continued if continued.is_file() else old_budget)
  if continued.is_file():
   prior=read(old_budget)
   if budget.get('baseline_budget_sha256')!=sha(old_budget) or budget['limits']!=prior['limits'] or budget['campaign']!=prior['campaign'] or budget['usage']['analysis_task_seconds_total']<prior['usage']['analysis_task_seconds_total']:
    raise ValueError('Invalid cumulative budget continuation')
  remaining=budget['limits']['analysis_task_seconds_total']-budget['usage']['analysis_task_seconds_total']
  if remaining<=0:raise ValueError('Campaign budget exhausted')
 else:
  if not all(auth.get(k) for k in ('approved_by','approved_on','decision_ref')) or auth.get('formal_research_explicitly_authorized') is not True:raise ValueError('Separate formal research authorization required')
  remaining=3600
 registered={Path(r['path']).resolve():r['sha256'] for r in auth.get('files',[])}
 def validate_engineering_input(value,upstream,seen=None):
  if cfg['mode']!='engineering_smoke':return
  seen=set() if seen is None else seen
  if isinstance(value,dict):
   for v in value.values():validate_engineering_input(v,upstream,seen)
  elif isinstance(value,list):
   for v in value:validate_engineering_input(v,upstream,seen)
  elif isinstance(value,str) and Path(value).is_absolute():
   p=Path(value).resolve()
   if any(p==u or p.is_relative_to(u) for u in upstream):return
   if p.is_relative_to(Path('/srv/CRC-PHIRE/databases')) or p.is_relative_to(ROOT/'assets'):return
   if not p.exists():raise ValueError('Missing engineering input '+str(p))
   if p.is_dir():
    for f in p.rglob('*'):
     if f.is_file():validate_engineering_input(str(f),upstream,seen)
   else:
    if p not in registered or sha(p)!=registered[p]:raise ValueError('Input outside registered engineering fixture '+str(p))
    if p.suffix=='.json' and p not in seen:
     seen.add(p);validate_engineering_input(read(p),upstream,seen)
 statefile=dest/'project_state.json'
 if dest.exists() and any(dest.iterdir()) and not statefile.is_file():raise ValueError('Nonempty output has no project checkpoint')
 dest.mkdir(parents=True,exist_ok=True)
 if statefile.exists() and not a.continue_prepared:raise ValueError('Use explicit --continue-prepared for this project checkpoint')
 source_hash=hashlib.sha256(json.dumps({str(p.relative_to(ROOT)):sha(p) for d in ('bin','modules','workflows','conf','assets') for p in (ROOT/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts},sort_keys=True).encode()).hexdigest()
 expected=dict(config_sha256=sha(a.config),authorization_sha256=sha(cfg['authorization_manifest']),source_sha256=source_hash,
  root_workflows={p.name:sha(p) for p in ROOT.glob('main*.nf')},database_registry_sha256=sha(ROOT/'assets/spark_databases.json'),environment_registry_sha256=sha(ROOT/'assets/spark_dependencies.json'),
  templates={x['template']:sha(x['template']) for n in cfg['stages'] for x in n['config_templates']})
 state=read(statefile) if statefile.exists() else dict(schema_version=1,identity=expected,stages={},mode=cfg['mode'])
 if state['identity']!=expected:raise ValueError('Source/config/authorization changed; choose new output')
 outputs={};started=time.monotonic();limit=min(3600,remaining,a.timeout)
 save(dest/'resolved_plan.json',plan);save(statefile,state)
 for node in cfg['stages']:
  ident=node['id'];prior=state['stages'].get(ident,{})
  if prior.get('status')=='completed':
   for name,h in prior['outputs'].items():
    if sha(Path(prior['outdir'])/name)!=h:raise ValueError('Completed output changed')
   outputs[ident]=Path(prior['outdir']);continue
  attempt=prior.get('attempt',0)+1;out=dest/ident/('attempt_'+str(attempt));par=resolve(node['parameters'],outputs)
  configroot=dest/'resolved_configs'/ident/str(attempt);configroot.mkdir(parents=True,exist_ok=False)
  for item in node['config_templates']:
   value=resolve(read(item['template']),outputs);validate_engineering_input(value,list(outputs.values()));fp=configroot/(item['parameter']+'.json');save(fp,value);par[item['parameter']]=str(fp)
  validate_engineering_input(par,[*outputs.values(),configroot])
  key,script,leading=KINDS[node['kind']];item,activation=spark_runtime.environment(key);prefix=Path(item['prefix'])
  cmd=([str(prefix/'bin/python'),str(ROOT/'bin'/script)] if script.endswith('.py') else [str(ROOT/'bin'/script)])+leading
  for name,value in par.items():
   if value is False or value is None:continue
   cmd.append('--'+name)
   if value is not True:cmd.extend(map(str,value if isinstance(value,list) else [value]))
  cmd+=['--outdir',str(out)]
  if node['kind']=='core':cmd+=['--work-dir',str(Path('/srv/CRC-PHIRE/scratch/crc-phage')/('project_'+hashlib.sha256(str(out).encode()).hexdigest()[:20]))]
  env=dict(os.environ,**activation,CRC_PHAGE_RUNTIME='spark-native')
  record=dict(status='running',attempt=attempt,outdir=str(out),argv=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),environment_lock_sha256=item['lock_sha256'])
  state['stages'][ident]=record;save(statefile,state);logs=dest/'logs';logs.mkdir(exist_ok=True)
  with (logs/(ident+'_'+str(attempt)+'.log')).open('w') as log:
   proc=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   try:
    while proc.poll() is None:
     if time.monotonic()-started>limit:raise RuntimeError('Time ceiling reached')
     if shutil.disk_usage(dest).free<300*2**30:raise RuntimeError('Free space below300GiB')
     if cfg['mode']=='engineering_smoke' and int(time.monotonic()-started)%15==0:
      size=sum(p.stat().st_size for p in dest.rglob('*') if p.is_file() and not p.is_symlink())
      if size>20*2**30:raise RuntimeError('Engineering project output exceeds20GiB')
     time.sleep(1)
   except BaseException:
    os.killpg(proc.pid,signal.SIGTERM)
    try:proc.wait(timeout=3)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
    record.update(status='stopped',exit_code=proc.returncode);save(statefile,state);raise
  record['exit_code']=proc.returncode;record['status']='completed' if proc.returncode==0 else 'failed'
  record['outputs']={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink()}
  save(statefile,state)
  if proc.returncode:raise RuntimeError('Stage '+ident+' failed; see '+str(logs))
  outputs[ident]=out
 state['status']='completed';state['elapsed_seconds']=time.monotonic()-started;save(statefile,state)
 cards=[]
 for n in plan['stages']:
  record=state['stages'][n['id']];folder=Path(record['outdir']);links=[]
  artifact_names=sorted(record['outputs'],key=lambda name:(not name.endswith(('.html','.pdf')),name))
  for name in artifact_names[:100]:
   if Path(name).suffix in ('.html','.pdf','.tsv','.json'):
    rel=str((folder/name).relative_to(dest));links.append('<a href="'+html.escape(rel,quote=True)+'">'+html.escape(name)+'</a>')
  cards.append('<section><h2>'+html.escape(n['id']+' / '+n['kind'])+'</h2><p>Status: '+record['status']+'; '+str(len(record['outputs']))+' hash-bound files. Navigation capped at100; full inventory in checkpoint.</p><ul>'+''.join('<li>'+v+'</li>' for v in links)+'</ul></section>')
 (dest/'report_project.html').write_text('<!doctype html><meta charset="utf-8"><title>CRC project run</title><style>body{font:16px system-ui;max-width:1100px;margin:2em auto}section{border-top:1px solid #ccc}</style><h1>CRC project run</h1><p>Mode: '+html.escape(cfg['mode'])+'. Engineering execution is not scientific calibration. Source host, predicted host and sample phasing are distinct. KO/cargo are not validated AMG; unknown is not zero. Only declared stages were run.</p><a href="project_state.json">Commands, exits, identities and full output inventory</a>'+''.join(cards))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['check','run']);p.add_argument('--config',required=True);p.add_argument('--outdir',required=True);p.add_argument('--continue-prepared',action='store_true');p.add_argument('--timeout',type=int,default=3600);a=p.parse_args()
 if a.action=='check':print(json.dumps(inspect(a.config,a.outdir)[1],indent=2))
 else:
  if not 1<=a.timeout<=3600:raise ValueError('Timeout cap3600seconds')
  run(a)
if __name__=='__main__':main()
