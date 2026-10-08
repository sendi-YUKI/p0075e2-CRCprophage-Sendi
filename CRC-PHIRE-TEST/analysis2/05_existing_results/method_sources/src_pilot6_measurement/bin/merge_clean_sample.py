"""Join already-QCed PE and SE readsets from one biological sample for assembly.

Primary PE fragments stay the quantification denominator. Separately deposited SE
reads and fastp orphans are retained explicitly as singletons for assembly.
This adapter does not claim that SE reads are paired or create a second patient.
"""
import argparse,copy,json
from pathlib import Path
from reads_io import sha,dump,ident

def merge(manifest_path,sample_id,output):
    manifest_path=Path(manifest_path);m=json.loads(manifest_path.read_text())
    units=[u for u in m['units'] if u['biological_sample_id']==sample_id]
    if len(units)!=len(m['units']) or not units:raise ValueError('Expected one sample only; prevent cross-patient merging')
    for key in ['subject_id','material_type','molecule','platform','preprocessing_id','profile_sha256','quality_encoding']:
        if len({u.get(key) for u in units})!=1:raise ValueError('Incompatible clean readset attribute: '+key)
    if units[0]['material_type']!='bulk_metagenome':raise ValueError('Adapter is for bulk metagenomes')
    paired=[u for u in units if u['layout']=='PE'];single=[u for u in units if u['layout']=='SE']
    if len(paired)!=1 or len(paired)+len(single)!=len(units):raise ValueError('First adapter requires exactly one PE readset and explicit SE readsets')
    if any(u['status'] not in ['success','zero_clean_reads'] for u in units):raise ValueError('QC readset incomplete')
    if any(u['metadata'].get('human_remove')!='true' for u in units):raise ValueError('All input readsets must pass the configured human-removal branch')
    seen=set()
    for u in units:
        for path,h in u['checksums'].items():
            if path in seen:raise ValueError('Duplicate clean input file')
            if sha(path)!=h:raise ValueError('Clean input checksum changed: '+path)
            seen.add(path)
    p=copy.deepcopy(paired[0]);p['unit_id']=ident(sample_id+'_assembly','unit_id')
    p['source_readsets']=sorted(r for u in units for r in u['source_readsets'])
    if len(p['source_readsets'])!=len(set(p['source_readsets'])):raise ValueError('Duplicate source readset')
    p['source_libraries']=sorted({u['library_id'] for u in units})
    if len(p['source_libraries'])!=1:raise ValueError('This engineering adapter requires one documented physical library')
    for u in single:
        if u.get('reads_2'):raise ValueError('SE readset cannot have mate2')
        p['singletons'].extend([u['reads_1'],*u.get('singletons',[])])
        p['singleton_counts'].extend([u['counts'],*u.get('singleton_counts',[])])
        p['checksums'].update(u['checksums'])
    p['merge_policy']='one_sample_one_library_primary_PE_plus_SE_singletons_v1'
    p['merge_source_manifest_sha256']=sha(manifest_path)
    p['readset_evidence']=[{k:u.get(k) for k in ['unit_id','layout','library_id','source_readsets','counts','raw_counts','human_filter','checksums','status']} for u in units]
    p['measurement_policy']='counts_are_primary_PE_fragments; singleton_counts_separate; assembly_uses_both'
    expected={x for x in [p['reads_1'],p['reads_2'],*p['singletons']] if x}
    if expected!=set(p['checksums']):raise ValueError('Merged paths and checksum keys disagree')
    result={k:v for k,v in m.items() if k!='units'}
    result.update(units=[p],merge_policy=p['merge_policy'],source_clean_manifest=str(manifest_path.resolve()),source_clean_manifest_sha256=sha(manifest_path))
    dump(output,result);return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--manifest',required=True);a.add_argument('--sample',required=True);a.add_argument('--output',required=True);x=a.parse_args();merge(x.manifest,x.sample,x.output)
