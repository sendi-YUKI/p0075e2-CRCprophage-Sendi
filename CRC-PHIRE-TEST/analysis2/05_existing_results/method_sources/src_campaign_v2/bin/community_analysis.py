#!/usr/bin/env python3
"""CRC-outcome community model adapter, clinical identity gate and holdout lock."""
import argparse,datetime,hashlib,json,math,os,subprocess,sys
from pathlib import Path
ROOT=Path(os.environ.get('CRC_SPARK_REPO','/srv/CRC-PHIRE/projects/crc-pipeline'))
sys.path.insert(0,str(ROOT/'bin'))
from mainline_common import read_table,write_table,save,sha,new_output,unique
import study_contract as study
import spark_runtime
from host_source import proof
FIELDS=['unit_id','biological_sample_id','participant_id','cohort_id','recruitment_id','family_id','feature_id','reference_host_id','status','reason','CRC','value','fragment_rpkm_all_clean','host_core_depth','host_core_breadth','clean_fragments_denominator','age','sex','lineage']
MODEL_KEYS=['schema_version','signal','sex_levels','minimum_n','minimum_outcome_count','minimum_host_breadth','firth_on_separation','detection_profile_id','catalog_version','mapping_reference_version']
def config(path,run=False):
 c=study.strict_json(study.read_small(path))
 required=set(MODEL_KEYS+['status','phase','cohorts','family','study_contract','covariates','measurement_directories','approval','holdout_lock'])
 if set(c)-{'module_directories'}!=required or c['schema_version']!=1:raise ValueError('Community config keys/schema mismatch')
 if c['status'] not in ('proposed','synthetic','frozen') or c['phase'] not in ('development','holdout'):raise ValueError('Unknown configuration status/phase')
 if c['signal'] not in ('detection','log1p_RPKM') or len(c['sex_levels'])<2 or len(set(c['sex_levels']))!=len(c['sex_levels']):raise ValueError('Explicit signal/sex reference levels required')
 if len(set(c['cohorts']))!=len(c['cohorts']) or not c['cohorts']:raise ValueError('Declare unique cohorts')
 if c['minimum_n'] is None or c['minimum_outcome_count'] is None or c['minimum_host_breadth'] is None:
  if run:raise ValueError('Scientific model choices pending')
 else:
  if type(c['minimum_n']) is not int or c['minimum_n']<3 or type(c['minimum_outcome_count']) is not int or c['minimum_outcome_count']<1 or not 0<=c['minimum_host_breadth']<=1:raise ValueError('Invalid scientific thresholds')
 if run and c['status']=='proposed':raise ValueError('Proposed configuration cannot execute')
 if c['status']=='frozen' and not all(c['approval'].get(k) for k in ('reviewer','date','decision_ref')):raise ValueError('Frozen config lacks scientific approval')
 return c
def model_hash(c):
 return hashlib.sha256(json.dumps({k:c[k] for k in MODEL_KEYS},sort_keys=True).encode()).hexdigest()
def family(path):
 rows=read_table(path,['family_id','feature_id','reference_host_id','quantification_group','module_id'])
 if not rows:raise ValueError('Declare complete F2/F3 family before fitting')
 unique(rows,'feature_id')
 for r in rows:
  if r['family_id'] not in ('F2','F3') or not r['reference_host_id']:raise ValueError('Family must name reference host')
  if r['family_id']=='F2' and (not r['quantification_group'] or r['module_id']):raise ValueError('F2 needs one explicit quantification group')
  if r['family_id']=='F3' and (not r['module_id'] or r['quantification_group']):raise ValueError('F3 needs one predefined module')
 return rows
def verify_measurements(path,c):
 p=Path(path);m=json.loads((p/'cohort_manifest.json').read_text())
 if m.get('status')!='completed' or m.get('measurement_contract')!='research_quantification_group_v1':raise ValueError('Research group measurement required')
 if (m.get('catalog_version'),m.get('mapping_reference_version'))!=(c['catalog_version'],c['mapping_reference_version']):raise ValueError('Mixed reference/catalogue')
 for n,h in m['outputs'].items():
  f=(p/n).resolve()
  if not f.is_relative_to(p.resolve()) or not f.is_file() or sha(f)!=h:raise ValueError('Changed measurement output: '+n)
 return m
def prepare(c,out):
 rows=family(c['family']);cfg,t,_,_=study.load_contract(c['study_contract'])
 if c['status']=='frozen' and cfg['status']!='frozen':raise ValueError('Clinical contract not frozen')
 if c['status']=='synthetic' and cfg['stage']!='preparation':raise ValueError('Artificial statistics require preparation clinical contract')
 errs,warns=study.relationships(cfg,t)
 if errs:raise ValueError('Clinical relationships: '+json.dumps(errs[:5]))
 samples=unique(t['samples'],'sample_id');people=unique(t['participants'],'participant_id')
 cohorts=unique(t['cohorts'],'cohort_id');studies=unique(t['studies'],'study_id');ev=unique(t['evidence'],'evidence_id')
 covs=unique(read_table(c['covariates'],['sample_id','age','sex','lineage','lineage_evidence_id']),'sample_id')
 role='development' if c['phase']=='development' else 'validation'
 for cid in c['cohorts']:
  co=cohorts.get(cid,{});st=studies.get(co.get('study_id'),{})
  if co.get('role')!=role or st.get('role')!=role:raise ValueError('Cohort/study role mismatch')
 signal=[];seen=set();receipts=[]
 for root in c['measurement_directories']:
  receipts.append(dict(path=root,sha256=sha(Path(root)/'cohort_manifest.json')));verify_measurements(root,c)
  signal.extend(read_table(Path(root)/'community_signals.tsv'))
  # Module producer is separately provenance-bound to this same group output.
  mp=Path(root)/'modules/module_manifest.json'
  external=[]
  for candidate in c.get('module_directories',[]):
   manifest=Path(candidate)/'module_manifest.json'
   mm_candidate=json.loads(manifest.read_text())
   if mm_candidate.get('community_signals_sha256')==sha(Path(root)/'community_signals.tsv'):external.append(manifest)
  if len(external)>1:raise ValueError('Multiple module outputs for the same measurement; choose one')
  if external:mp=external[0]
  if any(h['family_id']=='F3' for h in rows):
   if not mp.is_file():raise ValueError('F3 requested without actual module evidence producer output')
   mm=json.loads(mp.read_text())
   if mm.get('status')!='completed' or mm.get('community_signals_sha256')!=sha(Path(root)/'community_signals.tsv'):raise ValueError('Module scope changed')
   receipts.append(dict(path=str(mp),sha256=sha(mp)))
   if sha(mp.parent/'module_signals.tsv')!=mm['outputs']['module_signals.tsv']:raise ValueError('Module output changed')
   signal.extend(read_table(mp.parent/'module_signals.tsv'))
 normalized=[];audit=[];selected_people={};unitbindings={}
 for s in signal:
  if s['detection_profile_id']!=c['detection_profile_id']:continue
  sid=s['biological_sample_id'];sample=samples.get(sid,{})
  if sample.get('cohort_id') not in c['cohorts']:continue
  cid=sample['cohort_id'];pid=sample.get('participant_id');person=people.get(pid,{})
  identity=(s['unit_id'],sid,cid)
  if s['unit_id'] in unitbindings and unitbindings[s['unit_id']]!=identity:raise ValueError('Unit metadata collision')
  unitbindings[s['unit_id']]=identity
  reason=''
  if sample.get('baseline_selection')!='selected':reason='not_preselected_baseline'
  elif not pid or not proof(ev.get(person.get('identity_evidence_id'),{}),'participant',pid,'identity',pid):reason='independence_identity_not_reviewed'
  elif s.get('subject_id')!=pid:reason='measurement_clinical_participant_mismatch'
  elif not proof(ev.get(sample.get('evidence_id'),{}),'sample',sid,'disease_status',sample.get('disease_status')):reason='clinical_evidence_not_reviewed'
  elif sample.get('disease_status') not in ('CRC_confirmed','control_qualified'):reason='unknown_or_ineligible_clinical_label'
  if not reason:
   if pid in selected_people and selected_people[pid]!=s['unit_id']:raise ValueError('Multiple selected units per independent participant; predeclare selection/merge')
   selected_people[pid]=s['unit_id']
  cv=covs.get(sid,{})
  lineage=cv.get('lineage','unknown')
  if lineage not in ('unknown','','NA') and not proof(ev.get(cv.get('lineage_evidence_id'),{}),'sample',sid,'lineage',lineage):lineage='unknown'
  for h in rows:
   match=(s.get('module_id')==h['module_id']) if h['family_id']=='F3' else (not s.get('module_id') and s.get('quantification_group')==h['quantification_group'])
   if not match or s['reference_host_id']!=h['reference_host_id']:continue
   key=(s['unit_id'],h['feature_id'])
   if key in seen:raise ValueError('Duplicate group/module measurement across inputs')
   seen.add(key)
   normalized.append(dict(s,participant_id=pid,cohort_id=cid,recruitment_id=person.get('recruitment_id'),family_id=h['family_id'],feature_id=h['feature_id'],
       status='unassessable' if reason else s['status'],reason=reason or s['reason'],
       CRC=1 if sample.get('disease_status')=='CRC_confirmed' and not reason else 0 if sample.get('disease_status')=='control_qualified' and not reason else '',
       age=cv.get('age',''),sex=cv.get('sex',''),lineage=lineage))
 # Missing declared hypotheses are retained by R from family, not silently removed.
 if c['status']=='synthetic' and len(normalized)>2000:raise ValueError('Artificial statistics cap is 2000 model rows')
 write_table(out/'model_input.tsv',normalized,FIELDS);write_table(out/'family.tsv',rows,list(rows[0]))
 save(out/'input_manifest.json',dict(status='prepared',clinical_contract_sha256=sha(c['study_contract']),covariates_sha256=sha(c['covariates']),family_sha256=sha(c['family']),
  model_hash=model_hash(c),measurement_receipts=receipts,independent_participants=sorted(selected_people),recruitment_ids=sorted({r['recruitment_id'] for r in normalized if r.get('recruitment_id')}),
  cohort_ids=c['cohorts'],phase=c['phase'],scientific_calibration=False))
 return normalized
def lock(a):
 c=config(a.config,run=True)
 if c['phase']!='development':raise ValueError('Freeze using development config')
 review=study.strict_json(study.read_small(a.review))
 if review.get('independent_unseen_validation_confirmed') is not True or not all(review.get(k) for k in ('reviewer','date','provenance','validation_recruitment_ids','validation_participant_ids','validation_cohort_ids')):raise ValueError('Independent validation provenance incomplete')
 prior=json.loads((Path(a.development)/'community_analysis_manifest.json').read_text())
 for name,h in prior.get('outputs',{}).items():
  f=(Path(a.development)/name).resolve()
  if not f.is_relative_to(Path(a.development).resolve()) or not f.is_file() or sha(f)!=h:raise ValueError('Development output changed before lock')
 if not prior.get('outputs'):raise ValueError('Development output inventory missing')
 inp=json.loads((Path(a.development)/'input_manifest.json').read_text())
 if prior.get('status')!='completed' or inp['model_hash']!=model_hash(c):raise ValueError('Development result/model mismatch')
 for k in ('recruitment_ids','participant_ids','cohort_ids'):
  dev=inp['independent_participants'] if k=='participant_ids' else inp[k]
  if set(dev)&set(review['validation_'+k]):raise ValueError('Development/holdout overlap: '+k)
 dest=Path(a.outfile)
 if dest.exists():raise ValueError('Never overwrite a candidate lock')
 record=dict(schema_version=1,purpose='community_holdout_lock',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    status=c['status'],family_sha256=sha(c['family']),model_hash=model_hash(c),development_manifest_sha256=sha(Path(a.development)/'community_analysis_manifest.json'),
    development_participant_ids=inp['independent_participants'],development_recruitment_ids=inp['recruitment_ids'],development_cohort_ids=inp['cohort_ids'],
    family=family(c['family']),review=review)
 with dest.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
def execute(a):
 allocation=spark_runtime.allocation(require_lock=True);c=config(a.config,run=True)
 spark_runtime.validate_output_inputs(a.outdir,[a.config,c['family'],c['covariates'],Path(c['study_contract']).parent,*c['measurement_directories'],*c.get('module_directories',[]),c.get('holdout_lock')])
 out=new_output(a.outdir);normalized=prepare(c,out)
 if c['phase']=='holdout':
  if not c['holdout_lock']:raise ValueError('Holdout requires prior immutable candidate lock')
  lk=json.loads(Path(c['holdout_lock']).read_text())
  if lk.get('purpose')!='community_holdout_lock' or lk['model_hash']!=model_hash(c) or lk['family_sha256']!=sha(c['family']) or lk['status']!=c['status']:raise ValueError('Holdout lock/model/family mismatch')
  ins=json.loads((out/'input_manifest.json').read_text())
  for key,observed in [('participant_ids',ins['independent_participants']),('recruitment_ids',ins['recruitment_ids']),('cohort_ids',ins['cohort_ids'])]:
   if set(observed)&set(lk['development_'+key]) or not set(observed)<=set(lk['review']['validation_'+key]):raise ValueError('Unregistered/overlapping holdout identity '+key)
  c['holdout_lock_sha256']=sha(c['holdout_lock'])
 save(out/'config.json',c)
 _,item,prefix,_,activation=spark_runtime.environment_definition('r-sensitivity')
 spark_runtime.stable_identity('r-sensitivity',item,spark_runtime.read(spark_runtime.DEPENDENCIES))
 env=dict(os.environ,**activation,CONDA_PREFIX=str(prefix),R_HOME=str(prefix/'lib/R'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 env['PATH']=str(prefix/'bin')+':'+env.get('PATH','')
 cmd=[str(prefix/'bin/Rscript'),'--vanilla',str(ROOT/'bin/community_models.R'),str(out.resolve())]
 with (out/'R.stdout.log').open('w') as stdout,(out/'R.stderr.log').open('w') as stderr:
  r=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=env,timeout=a.timeout)
 save(out/'command.json',dict(argv=cmd,exit_code=r.returncode,timeout=a.timeout))
 if r.returncode:raise RuntimeError('Community R model failed')
 save(out/'community_analysis_manifest.json',dict(status='completed',phase=c['phase'],config_sha256=sha(a.config),
    model_hash=model_hash(c),r_script_sha256=sha(ROOT/'bin/community_models.R'),environment_lock_sha256=item['lock_sha256'],
    scientific_calibration=False,outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
 q=sub.add_parser('check');q.add_argument('--config',required=True)
 q=sub.add_parser('run');q.add_argument('--config',required=True);q.add_argument('--outdir',required=True);q.add_argument('--timeout',type=int,default=900)
 q=sub.add_parser('freeze-holdout');q.add_argument('--config',required=True);q.add_argument('--development',required=True);q.add_argument('--review',required=True);q.add_argument('--outfile',required=True)
 a=p.parse_args()
 if a.action=='check':
  c=config(a.config);f=family(c['family']);print(json.dumps(dict(status=c['status'],phase=c['phase'],family_size=len(f),analysis_started=False)))
 elif a.action=='freeze-holdout':lock(a)
 else:
  if not 1<=a.timeout<=3600:raise ValueError('Timeout bound')
  execute(a)
if __name__=='__main__':main()
