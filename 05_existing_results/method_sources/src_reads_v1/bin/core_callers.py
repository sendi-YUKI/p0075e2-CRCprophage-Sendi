#!/usr/bin/env python3
"""A4 caller adapters; canonical catalog is read-only. Runtime needs pinned callers image."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import runpy
import shutil
import subprocess
import sys
import time

VERSIONS = {'PhiSpy': '5.0.10', 'PhageBoost': '0.1.7+46239e3994be9cb5b2b714f8742fedd350d432c7'}
PHAGEBOOST_MIN_CONTIG_BP = 10000
POLICIES = {'PhiSpy': 'phispy_5.0.10_reported_1inclusive_v1', 'PhageBoost': 'phageboost_46239e3_reported_1inclusive_v1'}
HIT_FIELDS = ['caller_hit_id','genome_id','contig_id','original_contig_id','caller','tool_version','caller_original_id','candidate_type','start0','end0','strand','length','sequence_sha256','boundary_version','raw_start','raw_end','raw_coordinate_system','raw_sequence_sha256','raw_sequence_check','normalization_rule','caller_score','caller_confidence','integration_evidence','boundary_confidence','annotation_source','canonical_catalog_member']
COMPARISON_FIELDS = ['genome_id','contig_id','left_hit_id','left_caller','right_hit_id','right_caller','relation','overlap_bp','left_overlap_fraction','right_overlap_fraction','left_contains_right','right_contains_left','adjacent','gap_bp','right_minus_left_start_bp','right_minus_left_end_bp']

def digest(data):
    return hashlib.sha256(data).hexdigest()

def filehash(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def read_fasta(path):
    path=Path(path)
    with path.open('rb') as f: compressed=f.read(2)==b'\x1f\x8b'
    seqs={}; current=None
    with (gzip.open(path,'rt') if compressed else path.open()) as f:
        for raw in f:
            line=raw.strip()
            if not line:continue
            if line.startswith('>'):
                current=line[1:].split()[0]
                if not current or current in seqs:raise ValueError('Empty or duplicate FASTA ID')
                seqs[current]=[]
            elif current is None:raise ValueError('FASTA sequence before header')
            else:seqs[current].append(line.upper())
    return {k:''.join(v) for k,v in seqs.items()}

def read_tsv(path):
    with open(path,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path,fields,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def checked_slice(sequence,start0,end0):
    if not 0<=start0<end0<=len(sequence):raise ValueError(f'Invalid linear interval [{start0},{end0}) / {len(sequence)}; origin wrapping unsupported')
    return sequence[start0:end0]

def normalize_sequence(caller,source,start1,end1,raw_sequence):
    start1,end1=int(start1),int(end1)
    canonical=checked_slice(source,start1-1,end1)
    if raw_sequence==canonical:return canonical,'exact_canonical_match','convert_1based_inclusive_to_0halfopen'
    # Only the two source-proven writer formulas are accepted. Other mismatches fail.
    expected=source[start1:end1+1] if caller=='PhageBoost' else source[start1:end1]
    if raw_sequence != expected:raise ValueError(f'{caller}: raw FASTA matches neither canonical nor documented writer formula')
    return canonical,('known_writer_shift_right_1bp' if caller=='PhageBoost' else 'known_writer_omits_first_base'),('reextract_original_source_start_minus_1_to_end;raw_preserved')

def mapping_for(genome_id,contig_map,seqs):
    rows=[r for r in read_tsv(contig_map) if r['genome_id']==genome_id]
    by_orig={r['original_contig_id']:r for r in rows}
    if len(rows)!=len(by_orig) or set(seqs)!=set(by_orig):raise ValueError('Contig map must exactly match this original FASTA')
    return by_orig

def original_sequences(genome_id,contig_map,seqs):
    rows=[r for r in read_tsv(contig_map) if r['genome_id']==genome_id]
    original={r['original_contig_id']:r for r in rows}; canonical={r['contig_id']:r for r in rows}
    if len(rows)!=len(original) or len(rows)!=len(canonical):raise ValueError('Duplicate mapping IDs')
    if set(seqs)==set(original):result=seqs
    elif set(seqs)==set(canonical):result={canonical[k]['original_contig_id']:v for k,v in seqs.items()}
    else:raise ValueError('FASTA IDs match neither original nor canonical contig map')
    for k,v in result.items():
        r=original[k]
        if int(r['length'])!=len(v) or r['sequence_sha256']!=digest(v.encode()):raise ValueError('FASTA sequence disagrees with contig map checksum')
    return result

def validate_genbank(path,seqs):
    from Bio import SeqIO
    with open(path,'rb') as f:compressed=f.read(2)==b'\x1f\x8b'
    with (gzip.open(path,'rt') if compressed else open(path)) as f:records=list(SeqIO.parse(f,'genbank'))
    seen={}; receipt=[]
    for r in records:
        if r.id in seen:raise ValueError('Duplicate GenBank contig')
        seen[r.id]=str(r.seq).upper()
        cds=[f for f in r.features if f.type=='CDS']
        products=sum(bool(f.qualifiers.get('product')) for f in cds)
        if cds and not products:raise ValueError('PhiSpy needs genuine CDS product annotations; input is not assessed')
        receipt.append({'contig':r.id,'length':len(r),'sequence_sha256':digest(seen[r.id].encode()),'record_date':r.annotations.get('date'),'cds_count':len(cds),'cds_with_product':products,'phispy_contig_eligibility':('no_CDS_not_assessed' if not cds else 'below_min_contig_size_not_assessed' if len(r)<=5000 else 'annotated_input_eligible'),'phispy_min_contig_size_exclusive':5000,'annotation_comments':r.annotations.get('comment','')})
    if seen!=seqs:raise ValueError('GenBank sequences differ from input FASTA; no coordinate transfer allowed')
    if not any(x['cds_with_product'] for x in receipt):raise ValueError('PhiSpy needs genuine CDS product annotations; input is not assessed')
    return receipt

def invoke(command,log):
    with open(log,'w') as f:
        result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'Caller returned {result.returncode}; see {log}')
    return Path(log).read_text(errors='replace')

def make_hit(caller,gid,mapping,raw_id,contig,start,end,score,rawseq,source,annotation):
    seq,check,rule=normalize_sequence(caller,source,start,end,rawseq)
    version=POLICIES[caller]; start0=int(start)-1;end0=int(end)
    stable='|'.join([gid,mapping[contig]['contig_id'],str(start0),str(end0),version,digest(seq.encode())])
    return {'caller_hit_id':'pilot_'+digest(stable.encode())[:32], 'genome_id':gid,'contig_id':mapping[contig]['contig_id'],'original_contig_id':contig,'caller':caller,'tool_version':VERSIONS[caller],'caller_original_id':raw_id,'candidate_type':'provirus_locus_prediction','start0':start0,'end0':end0,'strand':'+','length':len(seq),'sequence_sha256':digest(seq.encode()),'boundary_version':version,'raw_start':start,'raw_end':end,'raw_coordinate_system':'1based_inclusive_reported','raw_sequence_sha256':digest(rawseq.encode()),'raw_sequence_check':check,'normalization_rule':rule,'caller_score':score,'caller_confidence':'caller_model_output_not_calibrated_between_callers','integration_evidence':'caller_predicted_locus','boundary_confidence':'caller_only_coordinate_convention_source_audited','annotation_source':annotation,'canonical_catalog_member':'false'},seq

def parse_phageboost(raw,gid,mapping,seqs,log):
    files=list(raw.glob('phages_*.gff'))
    if not files:
        if 'no phages found' in log and 'time after predictions:' in log:return []
        raise ValueError('Missing PhageBoost GFF without explicit successful zero-result marker')
    if len(files)!=1:raise ValueError('Expected one PhageBoost GFF')
    hits=[]
    for line in files[0].read_text().splitlines():
        if not line or line.startswith('#'):continue
        cells=line.split('\t')
        if len(cells)!=9:raise ValueError('Malformed PhageBoost GFF')
        contig,_,_,start,end,score,_,_,attributes=cells
        attrs=dict(x.split('=',1) for x in attributes.split(';'))
        rawid=attrs['phage_id'];name=contig.split('.')[0]+'_'+rawid
        fasta=read_fasta(raw/(name+'.fasta'))
        if len(fasta)!=1:raise ValueError('Expected one raw phage sequence')
        hits.append(make_hit('PhageBoost',gid,mapping,rawid,contig,start,end,score,next(iter(fasta.values())),seqs[contig],'internal_Pyrodigal_features_not_A2_gene_IDs'))
    return hits

def parse_phispy(raw,gid,mapping,seqs):
    coord=raw/'prophage_coordinates.tsv'
    if not coord.exists():raise ValueError('Missing PhiSpy coordinates')
    fastas=read_fasta(raw/'phage.fasta'); hits=[]
    for line in coord.read_text().splitlines():
        if not line:continue
        fields=line.split('\t'); rawid,contig,start,end=fields[:4]
        # Numeric output can be float text after optional attachment-site refinement.
        a,b=float(start),float(end)
        if not a.is_integer() or not b.is_integer():raise ValueError('Nonintegral PhiSpy coordinates')
        a,b=int(a),int(b)
        key=f'{contig}_{start}_{end}'
        if key not in fastas:raise ValueError('Missing PhiSpy raw sequence for coordinate record')
        hit,seq=make_hit('PhiSpy',gid,mapping,rawid,contig,a,b,'NA',fastas[key],seqs[contig],'sequence_verified_external_GenBank_CDS_products')
        if len(fields)>4 and any(x for x in fields[4:]):
            hit['boundary_confidence']='caller_attachment_refinement_convention_unresolved'
        hits.append((hit,seq))
    return hits

def run_caller(args):
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    caller='PhiSpy' if args.command=='run-phispy' else 'PhageBoost'
    status={'caller':caller,'genome_id':args.genome_id,'status':'running','exit_code':None,'candidate_count':None,'started_unix':time.time(),'scope':'engineering_pilot_not_accuracy_benchmark'}
    statusfile=out/'caller_status.json';write_json(statusfile,status)
    try:
        seqs=original_sequences(args.genome_id,args.contig_map,read_fasta(args.fasta))
        if not seqs:raise ValueError('Empty FASTA')
        mapping=mapping_for(args.genome_id,args.contig_map,seqs)
        raw=out/'raw'
        if raw.exists():raise ValueError('Output raw directory already exists; use a new task directory or Nextflow resume')
        raw.mkdir()
        caller_input=out/'caller_input.fna'
        with caller_input.open('w') as f:
            for name,seq in sorted(seqs.items()):f.write('>'+name+'\n'+seq+'\n')
        manifest={'caller':caller,'version':VERSIONS[caller],'boundary_policy':POLICIES[caller],'input_fasta_sha256':filehash(args.fasta),'contig_map_sha256':filehash(args.contig_map),'adapter_sha256':filehash(__file__),'coordinate_system':'0based_halfopen','canonical_catalog_changed':False,'models':model_manifest(caller)}
        if caller=='PhageBoost':
            eligible=[name for name,seq in seqs.items() if len(seq)>=PHAGEBOOST_MIN_CONTIG_BP]
            manifest['input_domain']={'min_contig_bp_inclusive':PHAGEBOOST_MIN_CONTIG_BP,
                'n_input_contigs':len(seqs),'n_eligible_contigs':len(eligible),
                'input_bp':sum(map(len,seqs.values())),
                'eligible_bp':sum(len(seqs[name]) for name in eligible)}
            if not eligible:
                reason='no_contig_at_or_above_PhageBoost_min_contig_bp'
                manifest.update(executed=False,reason=reason)
                write_json(out/'caller_manifest.json',manifest)
                write_tsv(out/'caller_hits_pilot.tsv',HIT_FIELDS,[])
                (out/'caller_candidates.fna').write_text('')
                status.update(status='not_assessed',exit_code=None,candidate_count=None,
                              reason=reason,finished_unix=time.time())
                write_json(statusfile,status)
                return
        if caller=='PhiSpy':
            if not args.genbank:raise ValueError('PhiSpy GenBank not supplied: not assessed; do not synthesize annotations')
            manifest['genbank']={'sha256':filehash(args.genbank),'records':validate_genbank(args.genbank,seqs)}
            command=[sys.executable,str(Path(__file__).resolve()),'_phispy',str(Path(args.genbank).resolve()),'-o',str(raw.resolve()),'--threads',str(args.threads),'--phage_genes','1','--min_contig_size','5000','--output_choice','511','--keep']
            manifest['random_seed']=1729;manifest['annotation_mode']='include_existing_products'
        else:
            command=['PhageBoost','-f',str(caller_input.resolve()),'-o',str(raw.resolve()),'-j','1','-t','0.9','-l','10','-g','5','-cs',str(PHAGEBOOST_MIN_CONTIG_BP),'-meta','0']
            manifest['parameters']={'threshold':0.9,'min_genes':10,'gap_genes':5,'min_contig_bp':PHAGEBOOST_MIN_CONTIG_BP,'force_meta':False,'repeat_refinement':False,'prediction_threads':1}
        manifest['command']=command
        write_json(out/'caller_manifest.json',manifest)
        log=invoke(command,out/'caller.log')
        results=parse_phispy(raw,args.genome_id,mapping,seqs) if caller=='PhiSpy' else parse_phageboost(raw,args.genome_id,mapping,seqs,log)
        results.sort(key=lambda x:x[0]['caller_hit_id'])
        if len({r['caller_hit_id'] for r,s in results})!=len(results):raise ValueError('Duplicate canonical caller hit ID')
        write_tsv(out/'caller_hits_pilot.tsv',HIT_FIELDS,[r for r,s in results])
        with (out/'caller_candidates.fna').open('w') as f:
            for row,sequence in results:f.write('>'+row['caller_hit_id']+'\n'+sequence+'\n')
        manifest['raw_files']=[{'path':str(p.relative_to(out)),'sha256':filehash(p),'bytes':p.stat().st_size} for p in sorted(raw.rglob('*')) if p.is_file()]
        manifest['raw_sequence_checks']={x:sum(r['raw_sequence_check']==x for r,s in results) for x in sorted({r['raw_sequence_check'] for r,s in results})}
        write_json(out/'caller_manifest.json',manifest)
        status.update(status='completed_hits' if results else 'completed_zero',exit_code=0,candidate_count=len(results),finished_unix=time.time())
        write_json(statusfile,status)
    except Exception as e:
        status.update(status='failed',exit_code=1,error=str(e),finished_unix=time.time());write_json(statusfile,status);raise

def model_manifest(caller):
    if caller=='PhiSpy':
        import PhiSpyModules
        base=Path(PhiSpyModules.__file__).parent/'data'
    else:
        import PhageBoost
        base=Path(PhageBoost.__file__).parent/'models'
    return [{'relative_path':str(p.relative_to(base)),'sha256':filehash(p),'bytes':p.stat().st_size} for p in sorted(base.rglob('*')) if p.is_file() and '__pycache__' not in str(p)]

def mark_unassessed(args):
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'caller_hits_pilot.tsv',HIT_FIELDS,[])
    (out/'caller_candidates.fna').write_text('')
    write_json(out/'caller_status.json',{'caller':'PhiSpy','genome_id':args.genome_id,'status':'not_assessed','exit_code':None,'candidate_count':None,'reason':args.reason,'scope':'engineering_pilot_not_accuracy_benchmark'})
    write_json(out/'caller_manifest.json',{'caller':'PhiSpy','expected_version':VERSIONS['PhiSpy'],'executed':False,'reason':args.reason,'adapter_sha256':filehash(__file__),'canonical_catalog_changed':False})

def relation(a,b):
    a0,a1,b0,b1=int(a['start0']),int(a['end0']),int(b['start0']),int(b['end0'])
    if a0>=a1 or b0>=b1:raise ValueError('Empty/reversed comparison interval')
    overlap=max(0,min(a1,b1)-max(a0,b0));ac=a0<=b0 and a1>=b1;bc=b0<=a0 and b1>=a1
    adjacent=a1==b0 or b1==a0
    kind='identical' if a0==b0 and a1==b1 else 'nested' if ac or bc else 'overlap' if overlap else 'adjacent' if adjacent else 'disjoint'
    return {'relation':kind,'overlap_bp':overlap,'left_overlap_fraction':f'{overlap/(a1-a0):.8f}','right_overlap_fraction':f'{overlap/(b1-b0):.8f}','left_contains_right':str(ac).lower(),'right_contains_left':str(bc).lower(),'adjacent':str(adjacent).lower(),'gap_bp':max(0,max(a0,b0)-min(a1,b1)),'right_minus_left_start_bp':b0-a0,'right_minus_left_end_bp':b1-a1}

def compare(args):
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    rows=[];statuses=[]
    for d in args.pilot_dirs:
        d=Path(d);status=json.loads((d/'caller_status.json').read_text());statuses.append(status)
        if status['status'] not in ('completed_hits','completed_zero','not_assessed'):raise ValueError('Cannot aggregate failed caller as zero')
        hits=read_tsv(d/'caller_hits_pilot.tsv')
        if status['status']=='not_assessed':
            if hits or status['candidate_count'] is not None:raise ValueError('Unassessed caller must retain unknown count and no hits')
        elif len(hits)!=status['candidate_count']:raise ValueError('Status candidate count mismatch')
        seqs=read_fasta(d/'caller_candidates.fna')
        if set(seqs)!={r['caller_hit_id'] for r in hits}:raise ValueError('Candidate FASTA and table IDs differ')
        for r in hits:
            if digest(seqs[r['caller_hit_id']].encode())!=r['sequence_sha256']:raise ValueError('Candidate hash mismatch')
        rows.extend(hits)
    for r in read_tsv(args.baseline):
        rows.append(dict(r,caller_hit_id=r['candidate_id'],caller='geNomad',tool_version='1.12.0',caller_original_id=r['caller_seq_name'],raw_start=r['start0'],raw_end=r['end0'],raw_coordinate_system='0based_halfopen_A1_canonical',raw_sequence_sha256=r['caller_sequence_sha256'],raw_sequence_check='A1_acceptance_verified',normalization_rule='baseline_unchanged',caller_score=r['virus_score'],caller_confidence=r['viral_identity_confidence'],annotation_source='geNomad_internal',canonical_catalog_member='true'))
    rows.sort(key=lambda r:r['caller_hit_id'])
    if len({r['caller_hit_id'] for r in rows})!=len(rows):raise ValueError('Duplicate caller IDs')
    comparisons=[]
    for i,a in enumerate(rows):
        for b in rows[i+1:]:
            if (a['genome_id'],a['contig_id'])!=(b['genome_id'],b['contig_id']) or a['caller']==b['caller']:continue
            comparisons.append(dict(genome_id=a['genome_id'],contig_id=a['contig_id'],left_hit_id=a['caller_hit_id'],left_caller=a['caller'],right_hit_id=b['caller_hit_id'],right_caller=b['caller'],**relation(a,b)))
    write_tsv(out/'caller_hits_pilot.tsv',HIT_FIELDS,rows)
    write_tsv(out/'caller_comparisons.tsv',COMPARISON_FIELDS,comparisons)
    write_tsv(out/'caller_sample_status.tsv',['genome_id','caller','status','exit_code','candidate_count','reason'],sorted(statuses,key=lambda r:(r['genome_id'],r['caller'])))
    write_json(out/'caller_comparison_manifest.json',{'schema_version':1,'baseline_sha256':filehash(args.baseline),'baseline_unchanged':True,'caller_counts':{c:sum(r['caller']==c for r in rows) for c in sorted({r['caller'] for r in rows})},'relations':{c:sum(r['relation']==c for r in comparisons) for c in sorted({r['relation'] for r in comparisons})},'comparison_scope':'all_cross_caller_pairs_on_same_genome_and_contig','adjacency_definition':'zero_gap_only','union_or_voting_rule':None,'scientific_benchmark':'not_performed','caller_runs':[{'status':s,'manifest_sha256':filehash(Path(d)/'caller_manifest.json')} for s,d in zip(statuses,args.pilot_dirs)]})

def main():
    if len(sys.argv)>1 and sys.argv[1]=='_phispy':
        import numpy as np
        random.seed(1729);np.random.seed(1729)
        executable=shutil.which('PhiSpy.py')
        sys.argv=[executable]+sys.argv[2:]
        runpy.run_path(executable,run_name='__main__');return
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    for name in ('run-phispy','run-phageboost'):
        s=sub.add_parser(name);s.add_argument('--genome-id',required=True);s.add_argument('--fasta',required=True);s.add_argument('--contig-map',required=True);s.add_argument('--genbank');s.add_argument('--outdir',required=True);s.add_argument('--threads',type=int,default=6)
    s=sub.add_parser('mark-unassessed');s.add_argument('--genome-id',required=True);s.add_argument('--reason',default='required_GenBank_CDS_product_annotations_not_supplied');s.add_argument('--outdir',required=True)
    s=sub.add_parser('compare');s.add_argument('--baseline',required=True);s.add_argument('--pilot-dirs',nargs='+',required=True);s.add_argument('--outdir',required=True)
    args=p.parse_args()
    if args.command=='compare':compare(args)
    elif args.command=='mark-unassessed':mark_unassessed(args)
    else:
        if not 1<=args.threads<=6:p.error('threads must be 1..6')
        run_caller(args)

if __name__=='__main__':main()
