#!/usr/bin/env python3
"""Consume scoped detection truth; no automatic simulations or LOD claims."""
import argparse,json,math
from pathlib import Path
from collections import defaultdict
from mainline_common import read_table,write_table,save,sha,new_output
def wilson(k,n):
 if not n:return ('','')
 z=1.959963984540054;p=k/n;den=1+z*z/n
 center=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0,center-half),min(1,center+half)
def execute(a):
 import spark_runtime
 spark_runtime.validate_output_inputs(a.outdir,[a.manifest,a.truth,a.detection])
 c=json.loads(Path(a.manifest).read_text())
 if c.get('purpose')!='detection_truth_review' or c.get('status') not in ('synthetic','reviewed'):raise ValueError('Explicit reviewed/synthetic truth required')
 if sha(a.truth)!=c['truth_sha256'] or sha(a.detection)!=c['detection_sha256']:raise ValueError('Calibration input identity mismatch')
 if c['status']=='reviewed' and not all(c.get(k) for k in ('reviewer','review_date','truth_source','assay_version')):raise ValueError('Unreviewed truth source')
 keys=['unit_id','quantification_group','detection_profile_id']
 truth={};counts=defaultdict(lambda:dict(TP=0,FN=0,TN=0,FP=0,unassessable=0,unknown_truth=0))
 for r in read_table(a.truth,keys+['truth','depth_condition','truth_source']):
  key=tuple(r[k] for k in keys)
  if key in truth:raise ValueError('Duplicate truth key')
  if r['truth'] not in ('0','1','unknown') or not r['truth_source']:raise ValueError('Truth source/value missing')
  truth[key]=r
 if c['status']=='synthetic' and len(truth)>2000:raise ValueError('Synthetic calibration interface size cap')
 seen=set()
 for r in read_table(a.detection,keys+['detection_status']):
  key=tuple(r[k] for k in keys)
  if key not in truth:continue
  if key in seen:raise ValueError('Duplicate detection result')
  seen.add(key);t=truth[key];group=(r['detection_profile_id'],t['depth_condition']);n=counts[group]
  if t['truth']=='unknown':n['unknown_truth']+=1
  elif r['detection_status'] not in ('detected','not_detected','not_detected_at_this_depth'):n['unassessable']+=1
  else:n[('T' if (r['detection_status']=='detected')==(t['truth']=='1') else 'F')+('P' if r['detection_status']=='detected' else 'N')]+=1
 if seen!=set(truth):raise ValueError('Truth rows lack matching predictions')
 rows=[]
 for (profile,depth),n in sorted(counts.items()):
  pos=n['TP']+n['FN'];neg=n['TN']+n['FP'];sl,sh=wilson(n['TP'],pos);pl,ph=wilson(n['TN'],neg)
  rows.append(dict(detection_profile_id=profile,depth_condition=depth,**n,sensitivity=n['TP']/pos if pos else '',specificity=n['TN']/neg if neg else '',
     sensitivity_CI_low=sl,sensitivity_CI_high=sh,specificity_CI_low=pl,specificity_CI_high=ph,interval='Wilson_95',LOD='not_estimated'))
 out=new_output(a.outdir);fields=['detection_profile_id','depth_condition','TP','FN','TN','FP','unassessable','unknown_truth','sensitivity','specificity','sensitivity_CI_low','sensitivity_CI_high','specificity_CI_low','specificity_CI_high','interval','LOD']
 write_table(out/'calibration_summary.tsv',rows,fields)
 save(out/'calibration_manifest.json',dict(status='completed',input_manifest_sha256=sha(a.manifest),truth_kind=c['status'],scientific_calibration=False,truth_comparison_performed=True,project_calibration_completed=False,
    generalization='Only supplied truth/depth/assay; no biological accuracy claim from synthetic truth',LOD='requires separately approved truth and prespecified depth series; not inferred from this summary',
    outputs={'calibration_summary.tsv':sha(out/'calibration_summary.tsv')}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('manifest','truth','detection','outdir'):p.add_argument('--'+k,required=True)
 execute(p.parse_args())
