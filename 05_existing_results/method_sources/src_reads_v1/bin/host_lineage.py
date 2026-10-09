#!/usr/bin/env python3
"""Declared supported-clade lineage assignment; outcome-independent and never calibrated by P."""
import argparse,json,os,sys
from pathlib import Path
ROOT=Path(os.environ.get('CRC_SPARK_REPO','/srv/CRC-PHIRE/projects/crc-pipeline'));sys.path.insert(0,str(ROOT/'bin'))
from mainline_common import read_table,write_table,save,sha,new_output
import spark_runtime
def execute(a):
 spark_runtime.allocation()
 from Bio import Phylo
 c=json.loads(Path(a.config).read_text())
 if c.get('status') not in ('synthetic','frozen') or c.get('method')!='supported_clade_max_patristic_diameter':raise ValueError('Explicit supported-clade lineage rule required')
 for k in ('rule_version','max_patristic_distance','minimum_support','minimum_members','outgroup'):
  if c.get(k) is None:raise ValueError('Lineage choices pending')
 if not 0<=c['minimum_support']<=100 or c['max_patristic_distance']<=0 or c['minimum_members']<2:raise ValueError('Invalid lineage rule')
 if c['status']=='frozen' and not all(c.get('approval',{}).get(k) for k in ('reviewer','date','decision_ref')):raise ValueError('Frozen lineage rule requires approval provenance')
 tree=Phylo.read(a.tree,'newick');tips=[t.name for t in tree.get_terminals()]
 declared=[r['genome_id'] for r in read_table(a.tips)]
 if len(tips)!=len(set(tips)) or set(tips)!=set(declared):raise ValueError('Lineage tree/tip identity mismatch')
 if any(x not in tips for x in c['outgroup']) or not c['outgroup']:raise ValueError('Explicit existing outgroup required')
 if any(t.branch_length is not None and t.branch_length<0 for t in tree.find_clades()):raise ValueError('Negative branch length')
 import hashlib
 tree.root_with_outgroup(*c['outgroup'])
 rows=[];assigned=set()
 def visit(node):
  leaves=node.get_terminals();ids=sorted(t.name for t in leaves)
  support=node.confidence
  diameter=max((tree.distance(x,y) for i,x in enumerate(leaves) for y in leaves[i+1:]),default=0)
  if len(ids)>=c['minimum_members'] and support is not None and support>=c['minimum_support'] and diameter<=c['max_patristic_distance']:
   ident='lin_'+hashlib.sha256(json.dumps([c['rule_version'],sha(a.tree),ids]).encode()).hexdigest()[:20]
   for gid in ids:rows.append(dict(genome_id=gid,lineage=ident,lineage_rule_version=c['rule_version'],status='assigned',support=support,diameter=diameter,reason=''));assigned.add(gid)
  else:
   for child in node.clades:visit(child)
 visit(tree.root)
 for gid in sorted(set(tips)-assigned):rows.append(dict(genome_id=gid,lineage='unknown',lineage_rule_version=c['rule_version'],status='insufficient_evidence',support='',diameter='',reason='no_supported_clade_meets_predeclared_rule'))
 spark_runtime.validate_output_inputs(a.outdir,[a.config,a.tree,a.tips])
 out=new_output(a.outdir);write_table(out/'lineages.tsv',rows,['genome_id','lineage','lineage_rule_version','status','support','diameter','reason'])
 save(out/'lineage_manifest.json',dict(status='completed',config_status=c['status'],approval=c.get('approval',{}),tree_sha256=sha(a.tree),tips_sha256=sha(a.tips),config_sha256=sha(a.config),method=c['method'],disease_data_used=False,scientific_calibration=False,outputs={'lineages.tsv':sha(out/'lineages.tsv')}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('tree','tips','config','outdir'):p.add_argument('--'+k,required=True)
 execute(p.parse_args())
