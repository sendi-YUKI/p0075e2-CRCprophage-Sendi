#!/usr/bin/env python3
"""Reference-coordinate SKA2 -> Gubbins for explicitly declared close-relative groups."""
import argparse,json,os,sys,shutil,hashlib
from pathlib import Path
ROOT=Path(os.environ.get('CRC_SPARK_REPO','/srv/CRC-PHIRE/projects/crc-pipeline'));sys.path.insert(0,str(ROOT/'bin'))
from mainline_common import fasta,write_fasta,read_table,write_table,sha,save,new_output,run,unique
from viral_catalog import load_units
import spark_runtime
def config(path,execute=False):
 c=json.loads(Path(path).read_text())
 required={'schema_version','status','scope','applicability_review','assembly_manifest','genome_ids','reference_genome_id','reference_contig','mobile_intervals','ska_k','gubbins','memory_gib'}
 if not required<=set(c) or set(c)-required-{'quality_policy'} or c['schema_version']!=1 or c['status'] not in ('proposed','synthetic','frozen'):raise ValueError('Recombination config schema mismatch')
 if c['scope']!='declared_recent_common_ancestor_group':raise ValueError('Gubbins must not run automatically over species-wide diversity')
 if execute:
  if c['status']=='proposed' or not c['reference_genome_id'] or not c['reference_contig'] or not c['mobile_intervals']:raise ValueError('Recombination method/reference choices pending')
  if c['status']=='frozen' and not all(c['applicability_review'].get(k) for k in ('reviewer','date','decision_ref')):raise ValueError('Close-relative applicability review required')
  if not c.get('quality_policy') or c['quality_policy'].get('status')!=c['status']:raise ValueError('Explicit matching recombination quality policy required')
  if c['status']=='frozen' and not c['applicability_review'].get('evidence'):raise ValueError('Hash-bound applicability assessment required')
  if len(set(c['genome_ids']))!=len(c['genome_ids']) or len(c['genome_ids'])<4:raise ValueError('Declare >=4 unique close relatives')
  if c['reference_genome_id'] not in c['genome_ids']:raise ValueError('Reference must be a declared group member')
  if not isinstance(c['ska_k'],int) or not 5<=c['ska_k']<=63 or c['ska_k']%2!=1:raise ValueError('SKA k must be explicit odd5..63')
  g=c['gubbins'];keys={'iterations','min_snps','min_window_size','max_window_size','p_value','trimming_ratio','filter_percentage','seed','tree_builder','model'}
  if set(g)!=keys or any(v is None for v in g.values()):raise ValueError('Explicit versioned Gubbins parameters required')
  if not 2<=g['iterations']<=20 or not 1<=g['min_snps'] or not 1<=g['min_window_size']<=g['max_window_size'] or not 0<g['p_value']<1 or not 0<g['trimming_ratio']<=1 or not 0<=g['filter_percentage']<=100 or not isinstance(g['seed'],int):raise ValueError('Gubbins parameter bounds')
  if g['tree_builder'] not in ('fasttree','iqtree') or g['model'] not in ('GTR','GTRGAMMA') or (g['tree_builder']=='fasttree' and g['model']=='GTR'):raise ValueError('Unreviewed tree configuration')
  if not 1<=c['memory_gib']<=46:raise ValueError('Memory ceiling')
 return c

def quality_report(seqs,c,source_manifest,original):
 ids=list(seqs);L=len(next(iter(seqs.values())))
 if not L or any(len(s)!=L for s in seqs.values()):raise ValueError('Empty/nonrectangular reference alignment')
 per=[dict(genome_id=g,called_bases=sum(b in 'ACGT' for b in s),called_fraction=sum(b in 'ACGT' for b in s)/L,
           input_contigs=len(original[g]),input_bp=sum(map(len,original[g].values()))) for g,s in seqs.items()]
 shared=sum(all(b in 'ACGT' for b in col) for col in zip(*seqs.values()))
 variable=sum(len({b for b in col if b in 'ACGT'})>1 for col in zip(*seqs.values()))
 pairs=[]
 for i,g in enumerate(ids):
  for h in ids[i+1:]:
   n=d=0
   for a,b in zip(seqs[g],seqs[h]):
    if a in 'ACGT' and b in 'ACGT':n+=1;d+=a!=b
   pairs.append(dict(genome_a=g,genome_b=h,comparable_sites=n,mismatches=d,distance=d/n if n else None))
 p=c.get('quality_policy');reasons=[]
 if not p:reasons.append('quality_policy_pending')
 else:
  for k in ('min_called_fraction','min_shared_fraction','max_pair_distance'):
   if not isinstance(p.get(k),(int,float)) or not 0<=p[k]<=1:raise ValueError('Invalid recombination quality policy '+k)
  if not isinstance(p.get('min_variable_sites'),int) or p['min_variable_sites']<1:raise ValueError('Invalid minimum variable sites')
  if p.get('status') not in ('proposed','synthetic','frozen') or not p.get('rule_version'):raise ValueError('Quality policy version/status required')
  if any(r['called_fraction']<p['min_called_fraction'] for r in per):reasons.append('insufficient_called_sequence')
  if shared/L<p['min_shared_fraction']:reasons.append('insufficient_shared_reference_columns')
  if any(r['distance'] is None or r['distance']>p['max_pair_distance'] for r in pairs):reasons.append('excess_divergence_or_no_pair_overlap')
  if variable<p['min_variable_sites']:reasons.append('insufficient_variable_sites')
  if p['status']!=c['status']:reasons.append('quality_policy_not_approved_for_run_status')
 review=c['applicability_review']
 if c['status']=='frozen':
  from genome_callability import bound_json,signature
  if not review.get('evidence'):reasons.append('applicability_evidence_missing')
  else:
   ev=bound_json(review['evidence'])
   if (ev.get('source_manifest_sha256')!=sha(source_manifest) or ev.get('quality_policy_sha256')!=signature(p) or
       ev.get('genome_ids')!=c['genome_ids'] or ev.get('reference_contig')!=c['reference_contig'] or
       ev.get('reference_genome_id')!=c['reference_genome_id'] or ev.get('scope')!=c['scope'] or
       not ev.get('recent_common_ancestor_basis') or not ev.get('limitations')):
    reasons.append('applicability_evidence_scope_mismatch')
 return dict(status='passed' if not reasons else 'not_applicable_or_pending',reasons=reasons,policy=p,
     reference_columns=L,shared_called_columns=shared,variable_columns=variable,genomes=per,pairs=pairs,
     scope='selected_reference_contig_only; reference_unmapped_and_other_contigs_not_assessed',scientific_calibration=False)
def execute(a):
 a.config=str(Path(a.config).resolve())
 alloc=spark_runtime.allocation();c=config(a.config,True)
 if c['memory_gib']*1024>alloc['memory_mib'] or not 1<=a.cpus<=min(6,alloc['cpus']):raise ValueError('Resource declaration exceeds allocation')
 source_manifest=str(Path(getattr(a,'assembly_manifest',None) or c['assembly_manifest']).resolve())
 units=load_units(source_manifest,verify=not bool(a.source_fastas))
 if a.source_fastas:
  if len(units)!=len(a.source_fastas):raise ValueError('Staged assembly count mismatch')
  for u,p in zip(units,a.source_fastas):
   if sha(p)!=u['sha256']:raise ValueError('Staged assembly checksum mismatch')
   u['assembly_fasta']=str(Path(p).resolve())
 byid=unique(units,'assembly_id')
 if not set(c['genome_ids'])<=set(byid):raise ValueError('Missing declared genomes')
 dest=Path(a.outdir).resolve();spark_runtime.validate_output_inputs(dest,[a.config,c['assembly_manifest'],c['mobile_intervals'],*[byid[g]['assembly_fasta'] for g in c['genome_ids']]])
 out=new_output(dest);original={g:fasta(byid[g]['assembly_fasta']) for g in c['genome_ids']};reference=original[c['reference_genome_id']]
 if c['reference_contig'] not in reference:raise ValueError('Reference contig missing')
 ref=reference[c['reference_contig']];L=len(ref)
 if c['status']=='synthetic' and (L>20000 or len(c['genome_ids'])>12):raise ValueError('Artificial recombination fixture bound')
 estimate=len(c['genome_ids'])*L*48+len(c['genome_ids'])**2*128
 if estimate>c['memory_gib']*2**30*.5:raise ValueError('Conservative alignment/sketch/matrix resource estimate exceeds half allocation')
 write_fasta(out/'reference.fna',{c['reference_contig']:ref})
 # Controlled local names avoid shell metacharacters in third-party scripts.
 local=out/'inputs';local.mkdir()
 paths=[]
 for g in c['genome_ids']:
  from mainline_common import safe_id;safe_id(g)
  p=local/(g+'.fna');p.symlink_to(Path(byid[g]['assembly_fasta']).resolve());paths.append((g,p))
 (out/'assemblies.tsv').write_text(''.join(g+'\t'+str(p)+'\n' for g,p in paths))
 old=Path.cwd();os.chdir(out)
 try:
  run(['ska','build','-o',out/'assemblies','-k',str(c['ska_k']),'-f',out/'assemblies.tsv','--threads',str(a.cpus)],out,'ska_build',timeout=a.timeout)
  run(['ska','map','-o',out/'reference_ordered.fna','--ambig-mask','--repeat-mask','--threads',str(a.cpus),out/'reference.fna',out/'assemblies.skf'],out,'ska_map',timeout=a.timeout)
  seqs=fasta(out/'reference_ordered.fna')
  if set(seqs)!=set(c['genome_ids']) or any(len(s)!=L for s in seqs.values()):raise ValueError('SKA reference-coordinate/sample identity mismatch')
  intervals=read_table(c['mobile_intervals'],['genome_id','contig_id','start0','end0'])
  mobile=set()
  for r in intervals:
   if r['genome_id']==c['reference_genome_id'] and r['contig_id']==c['reference_contig']:
    lo,hi=int(r['start0']),int(r['end0'])
    if not 0<=lo<hi<=L:raise ValueError('Mobile mask out of reference bounds')
    mobile.update(range(lo,hi))
  masked={g:''.join('N' if i in mobile or b not in 'ACGT-' else b for i,b in enumerate(s)) for g,s in seqs.items()}
  write_fasta(out/'gubbins_input.fna',masked)
  quality=quality_report(masked,c,source_manifest,original)
  save(out/'applicability.json',quality)
  if quality['status']!='passed':raise ValueError('Recombination applicability: '+','.join(quality['reasons']))
  gg=c['gubbins'];cmd=['run_gubbins.py','--prefix',out/'gubbins','--threads',str(a.cpus),'--no-cleanup']
  for k,v in gg.items():cmd+=['--'+k.replace('_','-'),str(v)]
  cmd.append(out/'gubbins_input.fna');run(cmd,out,'gubbins',timeout=a.timeout)
  gff=out/'gubbins.recombination_predictions.gff';tree=out/'gubbins.final_tree.tre'
  if not gff.is_file() or not tree.is_file():raise ValueError('Incomplete Gubbins outputs')
  remove=set(mobile);events=[]
  for line in gff.read_text().splitlines():
   if not line or line.startswith('#'):continue
   cc=line.split('\t')
   if len(cc)!=9:raise ValueError('Malformed recombination GFF')
   lo,hi=int(cc[3])-1,int(cc[4])
   if not 0<=lo<hi<=L:raise ValueError('Recombination interval out of reference bounds')
   remove.update(range(lo,hi));events.append(dict(start0=lo,end0=hi,attributes=cc[8],evidence_type='gubbins_predicted_recombination'))
  keep=[i for i in range(L) if i not in remove];final={g:''.join(s[i] for i in keep) for g,s in masked.items()}
  for name in ('masked','tree'): (out/name).mkdir()
  write_fasta(out/'masked/alignment.fna',final)
  tips=[dict(tip_id=g,genome_id=g,source_host_genome=g) for g in c['genome_ids']]
  for name in ('masked','tree'):write_table(out/name/'tree_tips.tsv',tips,['tip_id','genome_id','source_host_genome'])
  # Reject taxa silently filtered by Gubbins.
  from Bio import Phylo
  observed={leaf.name for leaf in Phylo.read(tree,'newick').get_terminals()}
  if observed!=set(c['genome_ids']):raise ValueError('Gubbins removed declared taxa; revise applicability/missingness explicitly')
  shutil.copyfile(tree,out/'tree/final_tree.treefile')
  save(out/'masked/mask_manifest.json',dict(status='completed',strategy='reference_ordered_SKA_repeat_mobile_plus_global_union_Gubbins',
      input_alignment_sha256=sha(out/'reference_ordered.fna'),config_sha256=sha(a.config),reference_genome=c['reference_genome_id'],reference_contig=c['reference_contig'],original_columns=L,final_columns=len(keep),removed_columns=len(remove),scientific_calibration=False))
  save(out/'tree/tree_status.json',dict(status='completed',input_alignment_sha256=sha(out/'gubbins_input.fna'),
      distance_alignment_sha256=sha(out/'masked/alignment.fna'),alignment_relationship='Gubbins_branch_specific_inference_vs_global_union_distance_mask',
      applicability_sha256=sha(out/'applicability.json'),method='Gubbins inferred clonal tree; distance/kinship uses conservative global union of predicted intervals'))
  write_table(out/'recombination_intervals.tsv',events,['start0','end0','attributes','evidence_type'])
  # A run-length crosswalk preserves original physical positions after removals.
  ranges=[]
  for finalpos,original in enumerate(keep):
   if ranges and ranges[-1]['reference_end0']==original:ranges[-1]['reference_end0']+=1;ranges[-1]['final_end0']+=1
   else:ranges.append(dict(reference_start0=original,reference_end0=original+1,final_start0=finalpos,final_end0=finalpos+1))
  write_table(out/'masked/reference_column_crosswalk.tsv',ranges,['reference_start0','reference_end0','final_start0','final_end0'])
  save(out/'recombination_manifest.json',dict(status='completed',method='Gubbins3.4.3_SKA2_0.5.1_reference_coordinate',
    scope=c['scope'],applicability_review=c['applicability_review'],applicability=quality,config_status=c['status'],config_sha256=sha(a.config),
    source_manifest_sha256=sha(source_manifest),mobile_intervals_sha256=sha(c['mobile_intervals']),scientific_calibration=False,
    limitations='Restricted to declared close relatives and one reference contig; reference coordinate mobile mask; not whole-species recombination inference; reference-unmapped regions unknown',
    outputs={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink()}))
  # Gubbins --no-cleanup retains dangling internal temporary symlinks.
  # Preserve raw work for audit; publish only regular, hash-bound deliverables.
  delivery=out/'deliverables';delivery.mkdir()
  for name in ('masked','tree'):shutil.copytree(out/name,delivery/name)
  for name in ('applicability.json','gubbins_input.fna','recombination_intervals.tsv','gubbins.recombination_predictions.gff','gubbins.final_tree.tre',
               'ska_build.command.json','ska_map.command.json','gubbins.command.json','gubbins.stdout.log','gubbins.stderr.log'):
   shutil.copyfile(out/name,delivery/name)
  record=json.loads((out/'recombination_manifest.json').read_text())
  record['raw_manifest_sha256']=sha(out/'recombination_manifest.json')
  record['outputs']={str(p.relative_to(delivery)):sha(p) for p in delivery.rglob('*') if p.is_file()}
  save(delivery/'recombination_manifest.json',record)
 finally:os.chdir(old)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['check','run']);p.add_argument('--config',required=True);p.add_argument('--outdir');p.add_argument('--assembly-manifest');p.add_argument('--source-fastas',nargs='*');p.add_argument('--cpus',type=int,default=1);p.add_argument('--timeout',type=int,default=900)
 a=p.parse_args()
 if a.action=='check':print(json.dumps(dict(config_status=config(a.config)['status'],analysis_started=False)))
 else:
  if not a.outdir or not 1<=a.timeout<=3600:raise ValueError('Output and bounded timeout required')
  execute(a)
