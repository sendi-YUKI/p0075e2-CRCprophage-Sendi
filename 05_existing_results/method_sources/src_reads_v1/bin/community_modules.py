#!/usr/bin/env python3
"""Predeclared KO-set carrier modules; measures annotated carrier groups, not gene activity."""
import argparse,json,os,sys
from pathlib import Path
ROOT=Path(os.environ.get('CRC_SPARK_REPO','/srv/CRC-PHIRE/projects/crc-pipeline'));sys.path.insert(0,str(ROOT/'bin'))
from collections import defaultdict
from mainline_common import read_table,write_table,save,sha,new_output,unique
def execute(a):
 c=json.loads(Path(a.config).read_text())
 import spark_runtime
 spark_runtime.validate_output_inputs(a.outdir,[a.config,a.measurement,*[x['directory'] for x in c.get('ko_sources',[])],*[x['research_bridge'] for x in c.get('ko_sources',[]) if x.get('research_bridge')]])
 if c.get('schema_version')!=1 or c.get('status') not in ('synthetic','frozen') or c.get('exposure')!='sum_group_RPKM_of_own_sequence_KO_set_carriers':
  raise ValueError('Explicit KO carrier module definition required; no generic functional abundance substitution')
 definitions=c['modules'];unique(definitions,'module_id')
 for m in definitions:
  if not m.get('required_kos') or len(m['required_kos'])!=len(set(m['required_kos'])):raise ValueError('Unique required KO set needed')
 p=Path(a.measurement);manifest=json.loads((p/'cohort_manifest.json').read_text())
 if manifest.get('measurement_contract')!='research_quantification_group_v1':raise ValueError('Group measurement required')
 for name,h in manifest['outputs'].items():
  target=(p/name).resolve()
  if not target.is_relative_to(p.resolve()) or sha(target)!=h:raise ValueError('Changed group measurement')
 groups=read_table(p/'quantification_group_members.tsv')
 # Required normalized own-sequence evidence includes producing manifests.
 own=defaultdict(set);assessed=set();source_receipts=[]
 for source in c['ko_sources']:
  root=Path(source['directory'])
  if source['kind']=='viral':
   receipt=root/'viral_ko_manifest.json';r=json.loads(receipt.read_text())
   for name,h in r['outputs'].items():
    if sha(root/name)!=h:raise ValueError('Changed viral KO output')
   for x in read_table(root/'viral_member_ko_evidence.tsv'):
    assessed.add(x['viral_sequence_id'])
    if x['passes_kofamscan_rule']=='true':own[x['viral_sequence_id']].add(x['ko'])
  elif source['kind']=='canonical':
   bridge=Path(source['research_bridge']);bm=json.loads((bridge/'catalog_manifest.json').read_text())
   cross=unique(read_table(bridge/'research_candidate_sequence_crosswalk.tsv'),'source_candidate_id')
   for name,h in bm['outputs'].items():
    if sha(bridge/name)!=h:raise ValueError('Changed research bridge')
   receipt=root/'kofam_summary.json';r=json.loads(receipt.read_text())
   if r.get('status')!='completed':raise ValueError('KOfam summary incomplete')
   # The upstream summary now contains output hashes; older summaries need regeneration.
   if not r.get('outputs'):raise ValueError('Canonical KO summary lacks content inventory; regenerate with current code')
   for name,h in r['outputs'].items():
    if sha(root/name)!=h:raise ValueError('Changed canonical KO output')
   for x in read_table(root/'member_ko_evidence.tsv'):
    if x['candidate_id'] not in cross:continue
    sid=cross[x['candidate_id']]['viral_sequence_id'];assessed.add(sid)
    if x['ko_evidence_status']=='accepted' and x['relation']=='fully_inside':own[sid].add(x['ko'])
  else:raise ValueError('Unsupported KO source')
  source_receipts.append(dict(path=str(receipt),sha256=sha(receipt)))
 bygroup=defaultdict(set)
 for g in groups:bygroup[g['quantification_group']].add(g['viral_sequence_id'])
 membership=[];module_groups=defaultdict(set)
 for m in definitions:
  required=set(m['required_kos'])
  for qg,sids in bygroup.items():
   state='carrier' if all(sid in assessed and required<=own[sid] for sid in sids) else 'not_qualified_or_unknown'
   membership.append(dict(module_id=m['module_id'],quantification_group=qg,status=state,required_kos=';'.join(sorted(required)),scope='own_sequence_KO_set_cooccurrence_not_activity_or_complete_pathway'))
   if state=='carrier':module_groups[m['module_id']].add(qg)
 signal=read_table(p/'community_signals.tsv');buckets=defaultdict(list)
 for r in signal:buckets[(r['unit_id'],r['reference_host_id'],r['detection_profile_id'])].append(r)
 output=[]
 for (unit,host,profile),rr in buckets.items():
  byq=unique(rr,'quantification_group')
  for m in definitions:
   qq=module_groups[m['module_id']] & set(byq)
   if not qq:continue
   ss=[byq[q] for q in sorted(qq)];base=dict(ss[0]);ok=all(s['status']=='assessed' for s in ss)
   if len({(s['host_core_depth'],s['host_core_breadth'],s['clean_fragments_denominator']) for s in ss})!=1:raise ValueError('Module host/denominator mismatch')
   base.update(module_id=m['module_id'],quantification_group='',status='assessed' if ok else 'unassessable',
      reason='' if ok else 'one_or_more_carrier_groups_unassessable',
      value=int(any(float(s['value'])==1 for s in ss)) if ok else '',
      fragment_rpkm_all_clean=sum(float(s['fragment_rpkm_all_clean']) for s in ss) if ok else '',
      member_group_count=len(qq),exposure_definition=c['exposure'])
   output.append(base)
 out=new_output(a.outdir)
 fields=list(signal[0]) if signal else ['unit_id','biological_sample_id','subject_id','reference_host_id','quantification_group','detection_profile_id','status','reason','value','fragment_rpkm_all_clean','host_core_depth','host_core_breadth','clean_fragments_denominator']
 write_table(out/'module_signals.tsv',output,fields+['module_id','member_group_count','exposure_definition'])
 write_table(out/'module_group_membership.tsv',membership,['module_id','quantification_group','status','required_kos','scope'])
 save(out/'module_manifest.json',dict(status='completed',config_sha256=sha(a.config),community_signals_sha256=sha(p/'community_signals.tsv'),
   source_receipts=source_receipts,definitions=definitions,exposure_definition=c['exposure'],scientific_calibration=False,
   limitations='Carrier-group abundance, not KO gene coverage, pathway activity or AMG; groups without confirmed own-sequence set are not assigned to module',
   outputs={f.name:sha(f) for f in out.iterdir() if f.is_file()}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('config','measurement','outdir'):p.add_argument('--'+k,required=True)
 execute(p.parse_args())
