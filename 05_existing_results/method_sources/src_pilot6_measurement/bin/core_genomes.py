#!/usr/bin/env python3
"""A2 assembly QC, canonical bacterial genes and reference-limited taxonomy.

Workers use stdlib; analysis programs use the selected pinned container or native
locked runtime. No quality estimate is invented for missing/failed output.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from urllib.parse import quote

from phageflow import ContractError, read_fasta, read_tsv, write_tsv, write_fasta, write_json, file_hash, sequence_hash, reverse_complement, unique_index

VERSION = 'core_genomes_v1'
CHECKM2_IMAGE = 'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d@sha256:fcd923f5d60786deede345f52e8ad3ec73954ca8ea2eb914b56e7d38c8e1159b'
FASTANI_IMAGE = 'quay.io/biocontainers/fastani:1.34--hb66fcc3_7@sha256:0723121d8154f7e4917304142824d12bff93494d2b2ed4ab072d51eccf6c8494'
GENE_FIELDS = ['genome_id','contig_id','original_contig_id','gene_id','tool_gene_id','tool_fasta_id','start0','end0','strand','protein_sha256','cds_sha256','partial','gene_order','translation_table','gene_caller','caller_mode','protein_length','cds_length','translation_verified','slice_verified']
QC_FIELDS = ['genome_id','checkm2_status','completeness','contamination','completeness_model','genome_qc_status','eligibility','high_confidence','eligibility_reason','completeness_min','contamination_max','high_completeness_min','high_contamination_max','length','contigs','N50','fragmentation_status','fragmentation_n50_warning_below','replicon_status','mag_status','duplicate_status','duplicate_of','assembly_content_sha256','ncbi_declared_organism','ncbi_declared_taxid','ncbi_declared_assembly_level','ncbi_assembly_accession','metadata_source','qc_policy_version']
TAX_FIELDS = ['genome_id','ncbi_declared_organism','ncbi_declared_taxid','taxonomy_status','reference_supported_species','reference_supported_species_taxid','reference_id','ani_forward_percent','ani_reverse_percent','query_fragment_fraction','reference_fragment_fraction','reason','reference_catalog_version','ani_min_percent','fragment_fraction_min','taxonomy_policy_version']
PAIR_FIELDS = ['genome_id','reference_id','reference_species','reference_species_taxid','status','ani_forward_percent','ani_reverse_percent','query_fragment_fraction','reference_fragment_fraction','passes_thresholds','query_sha256','reference_sha256']
CODONS = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))

def normalizer_runtime_provenance(legacy_image, container_field='container'):
    """Identify this normalizer; never attribute an upstream tool to its environment."""
    if os.environ.get('CRC_PHAGE_RUNTIME')!='spark-native' and not os.environ.get('CRC_SPARK_ENVIRONMENT'):
        return {container_field:legacy_image}
    return {container_field:None, 'legacy_container_reference':legacy_image,
            'container_executed':False, 'runtime':'native-locked-prefix',
            'normalizer_environment':os.environ.get('CRC_SPARK_ENVIRONMENT'),
            'normalizer_environment_lock_sha256':os.environ.get('CRC_SPARK_ENV_LOCK_SHA256'),
            'normalizer_runtime_fingerprint':os.environ.get('CRC_SPARK_ENV_FINGERPRINT')}


def fasta_any(path):
    """Read protein/CDS FASTA without DNA-only validation, keeping headers."""
    result = {}; header = None; parts = []
    def finish():
        if header is None: return
        key = header.split()[0]
        if key in result: raise ContractError(f'Duplicate FASTA key {key}')
        sequence = ''.join(parts).upper()
        if not sequence: raise ContractError(f'Empty FASTA sequence {key}')
        result[key] = (header, sequence)
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            finish(); header = line[1:]; parts=[]
        elif line.strip():
            if header is None: raise ContractError('Sequence before FASTA header')
            parts.append(line.strip())
    finish(); return result

def attrs(text):
    result={}
    for field in text.strip().rstrip(';').split(';'):
        if not field: continue
        key, value=field.split('=',1)
        if key in result: raise ContractError('Duplicate attribute '+key)
        result[key]=value
    return result

def translate(cds, start_type):
    cds=cds.upper()
    if len(cds)%3: raise ContractError('CDS length is not divisible by three')
    peptide=''.join(CODONS.get(cds[p:p+3], 'X') for p in range(0,len(cds),3))
    if start_type != 'Edge' and peptide: peptide='M'+peptide[1:]
    return peptide[:-1] if peptide.endswith('*') else peptide

def canonical_gene_id(genome_id, contig, start0, end0, strand, cds):
    identity=json.dumps([VERSION,genome_id,contig,start0,end0,strand,sequence_hash(cds)],separators=(',',':'))
    return 'gene_'+hashlib.sha256(identity.encode()).hexdigest()[:32]

def prepared_contract(prepared, genome_id):
    prepared=Path(prepared); dna=read_fasta(prepared/'genome.fna')
    _, rows=read_tsv(prepared/'contig_map.tsv',['genome_id','contig_id','original_contig_id','length','sequence_sha256'])
    mapping=unique_index(rows,'contig_id','contig map')
    if set(dna)!=set(mapping): raise ContractError('FASTA and mapping contig sets differ')
    for c,(_,s) in dna.items():
        if mapping[c]['genome_id']!=genome_id or int(mapping[c]['length'])!=len(s) or mapping[c]['sequence_sha256']!=sequence_hash(s):
            raise ContractError('Prepared contig metadata/hash mismatch '+c)
    return dna,mapping

def genes_normalize(args):
    dna,mapping=prepared_contract(args.prepared,args.genome_id)
    raw=Path(args.raw); out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    proteins=fasta_any(raw/'prodigal.faa');cdss=fasta_any(raw/'prodigal.ffn')
    if set(proteins)!=set(cdss): raise ContractError('Prodigal protein/CDS key mismatch')
    by_tool={}
    for fasta_id,(header,seq) in proteins.items():
        parts=header.split(' # ')
        if len(parts)!=5: raise ContractError('Unexpected Prodigal FASTA header')
        meta=attrs(parts[4]); tool_id=meta['ID']
        if tool_id in by_tool: raise ContractError('Duplicate Prodigal internal ID')
        by_tool[tool_id]=(fasta_id,meta,int(parts[1])-1,int(parts[2]),'+' if parts[3]=='1' else '-',seq)
    rows=[]; sequences={}; seen=set()
    for line in (raw/'prodigal.gff').read_text().splitlines():
        if not line or line.startswith('#'):continue
        fields=line.split('\t')
        if len(fields)!=9 or fields[2]!='CDS' or fields[7]!='0':raise ContractError('Unsupported Prodigal GFF record')
        contig=fields[0];start0=int(fields[3])-1;end0=int(fields[4]);strand=fields[6];meta=attrs(fields[8]);tid=meta['ID']
        if tid in seen or tid not in by_tool or contig not in dna:raise ContractError('Duplicate/missing GFF identity or contig')
        seen.add(tid);fid,fmeta,fs,fe,fstrand,protein=by_tool[tid]
        if (fs,fe,fstrand)!=(start0,end0,strand) or fmeta['partial']!=meta['partial']:raise ContractError('GFF and protein coordinate mismatch')
        if not(0<=start0<end0<=len(dna[contig][1])) or strand not in ('+','-'):raise ContractError('Invalid gene coordinates')
        cds=dna[contig][1][start0:end0]
        if strand=='-':cds=reverse_complement(cds)
        if cds!=cdss[fid][1]:raise ContractError('Gene original-source slice mismatch '+tid)
        translated=translate(cds,meta['start_type']); peptide=protein[:-1] if protein.endswith('*') else protein
        if translated!=peptide or '*' in peptide:raise ContractError('Gene translation mismatch '+tid)
        gene_id=canonical_gene_id(args.genome_id,contig,start0,end0,strand,cds)
        sequences[gene_id]=(peptide,cds)
        rows.append(dict(genome_id=args.genome_id,contig_id=contig,original_contig_id=mapping[contig]['original_contig_id'],gene_id=gene_id,tool_gene_id=tid,tool_fasta_id=fid,start0=start0,end0=end0,strand=strand,protein_sha256=sequence_hash(peptide),cds_sha256=sequence_hash(cds),partial=meta['partial'],translation_table=11,gene_caller='Prodigal2.6.3',caller_mode='single',protein_length=len(peptide),cds_length=len(cds),translation_verified='true',slice_verified='true'))
    if seen!=set(by_tool):raise ContractError('Unmapped Prodigal FASTA records')
    rows.sort(key=lambda r:(r['contig_id'],r['start0'],r['end0'],r['strand'],r['gene_id']))
    counts={}
    for row in rows:
        counts[row['contig_id']]=counts.get(row['contig_id'],0)+1;row['gene_order']=counts[row['contig_id']]
    unique_index(rows,'gene_id','canonical genes')
    write_tsv(out/'gene_table.tsv',GENE_FIELDS,rows)
    write_tsv(out/'crosswalk.tsv',['genome_id','gene_id','tool_gene_id','tool_fasta_id','contig_id','original_contig_id'],rows)
    write_fasta(out/'genes.faa',[(r['gene_id'],sequences[r['gene_id']][0]) for r in rows])
    write_fasta(out/'genes.ffn',[(r['gene_id'],sequences[r['gene_id']][1]) for r in rows])
    with (out/'genes.gff').open('w') as w:
        w.write('##gff-version 3\n')
        for row in rows:
            attributes=';'.join(k+'='+quote(str(v),safe='._-') for k,v in {'ID':row['gene_id'],'tool_gene_id':row['tool_gene_id'],'partial':row['partial'],'transl_table':11}.items())
            w.write('\t'.join(map(str,[row['contig_id'],'Prodigal2.6.3','CDS',row['start0']+1,row['end0'],'.',row['strand'],0,attributes]))+'\n')
    write_json(out/'gene_manifest.json',dict(schema_version=VERSION,genome_id=args.genome_id,gene_count=len(rows),caller='Prodigal',version='2.6.3',mode='single',translation_table=11,**normalizer_runtime_provenance(CHECKM2_IMAGE),coordinates='0-based half-open; GFF retains standard 1-based inclusive',input_sha256=file_hash(Path(args.prepared)/'genome.fna'),raw_files={p.name:file_hash(p) for p in sorted(raw.glob('prodigal.*')) if p.is_file()},outputs={p.name:file_hash(p) for p in sorted(out.glob('gene*')) if p.is_file()}))

def numeric(value,name,upper=100):
    if value is None or str(value).strip().lower() in ('','na','nan','none','unknown'):return None
    result=float(value)
    if not math.isfinite(result) or result<0 or result>upper:raise ContractError('Invalid '+name)
    return result

def qc_decision(completeness,contamination,complete_min=90,contam_max=5,high_complete=95,high_contam=2):
    if completeness is None or contamination is None:return 'unknown','unknown','unknown','missing_quality_estimate'
    reasons=[]
    if completeness<complete_min:reasons.append('completeness_below_threshold')
    if contamination>contam_max:reasons.append('contamination_above_threshold')
    passed=not reasons
    return ('pass' if passed else 'fail','eligible' if passed else 'ineligible','true' if passed and completeness>=high_complete and contamination<=high_contam else 'false','thresholds_met' if passed else ';'.join(reasons))

def read_metadata(path,genome_id):
    if not path:return {}
    _,rows=read_tsv(Path(path),['genome_id']); index=unique_index(rows,'genome_id','genome metadata')
    return index.get(genome_id,{})

def qc_normalize(args):
    dna,_=prepared_contract(args.prepared,args.genome_id)
    _,stats=read_tsv(Path(args.prepared)/'genome_stats.tsv',['genome_id','length','contigs','N50'])
    if len(stats)!=1 or stats[0]['genome_id']!=args.genome_id:raise ContractError('Genome statistic identity mismatch')
    meta=read_metadata(args.metadata,args.genome_id); s=stats[0]
    if not(0<=args.completeness_min<=args.high_completeness_min<=100 and 0<=args.high_contamination_max<=args.contamination_max<=100):raise ContractError('Invalid quality policy thresholds')
    row=dict(genome_id=args.genome_id,checkm2_status='not_assessed',completeness='unknown',contamination='unknown',completeness_model='unknown',genome_qc_status='not_assessed',eligibility='unknown',high_confidence='unknown',eligibility_reason='checkm2_not_assessed',completeness_min=args.completeness_min,contamination_max=args.contamination_max,high_completeness_min=args.high_completeness_min,high_contamination_max=args.high_contamination_max,length=s['length'],contigs=s['contigs'],N50=s['N50'],fragmentation_status='warning_low_N50' if int(s['N50'])<args.fragmentation_n50 else 'no_N50_warning',fragmentation_n50_warning_below=args.fragmentation_n50,replicon_status=meta.get('replicon_status','unknown'),mag_status=meta.get('mag_status','unknown'),duplicate_status='not_assessed',duplicate_of='',assembly_content_sha256=hashlib.sha256('\n'.join(sorted(sequence_hash(v[1]) for v in dna.values())).encode()).hexdigest(),ncbi_declared_organism=meta.get('ncbi_declared_organism','unknown'),ncbi_declared_taxid=meta.get('ncbi_declared_taxid','unknown'),ncbi_declared_assembly_level=meta.get('ncbi_declared_assembly_level','unknown'),ncbi_assembly_accession=meta.get('ncbi_assembly_accession','unknown'),metadata_source=meta.get('metadata_source','not_provided'),qc_policy_version='checkm2_90_5_configurable_v1')
    if args.checkm2_report:
        _,reports=read_tsv(Path(args.checkm2_report),['Name','Completeness','Contamination','Completeness_Model_Used'])
        if len(reports)!=1 or reports[0]['Name']!=args.checkm2_name:raise ContractError('CheckM2 output genome identity/count mismatch')
        q=reports[0];complete=numeric(q['Completeness'],'completeness');contam=numeric(q['Contamination'],'contamination',100000)
        status,eligible,high,reason=qc_decision(complete,contam,args.completeness_min,args.contamination_max,args.high_completeness_min,args.high_contamination_max)
        row.update(checkm2_status='completed',completeness=complete if complete is not None else 'unknown',contamination=contam if contam is not None else 'unknown',completeness_model=q['Completeness_Model_Used'],genome_qc_status=status,eligibility=eligible,high_confidence=high,eligibility_reason=reason)
    out=Path(args.outdir);write_tsv(out/'genome_qc.tsv',QC_FIELDS,[row])
    write_json(out/'qc_manifest.json',dict(genome_id=args.genome_id,**normalizer_runtime_provenance(CHECKM2_IMAGE,'checkm2_container'),checkm2_report_sha256=file_hash(Path(args.checkm2_report)) if args.checkm2_report else None,metadata_sha256=file_hash(Path(args.metadata)) if args.metadata else None,policy=row))

def deduplicate(args):
    rows=[]
    for filename in args.inputs:
        _,part=read_tsv(Path(filename),QC_FIELDS);rows.extend(part)
    unique_index(rows,'genome_id','genome QC');groups={}
    for row in rows:groups.setdefault(row['assembly_content_sha256'],[]).append(row)
    for group in groups.values():
        def rank(r):return (r['eligibility']!='eligible',-(float(r['completeness']) if r['completeness']!='unknown' else -1),float(r['contamination']) if r['contamination']!='unknown' else 1e9,r['genome_id'])
        group.sort(key=rank)
        for i,row in enumerate(group):
            row['duplicate_status']='unique' if len(group)==1 else 'representative' if i==0 else 'exact_sequence_duplicate'
            if i:row['duplicate_of']=group[0]['genome_id'];row['eligibility']='ineligible';row['eligibility_reason']+=';exact_sequence_duplicate'
    write_tsv(Path(args.outdir)/'genome_qc.tsv',QC_FIELDS,sorted(rows,key=lambda r:r['genome_id']))

def docker(command,image,mounts,log,threads=6):
    argv=['docker','run','--rm','--cpus',str(threads),'--memory','12g','--user',f'{os.getuid()}:{os.getgid()}','-e','OMP_NUM_THREADS='+str(threads),'-e','TF_NUM_INTRAOP_THREADS='+str(threads),'-e','TF_NUM_INTEROP_THREADS=1']
    for host,container,mode in mounts:argv.extend(['-v',f'{Path(host).resolve()}:{container}:{mode}'])
    argv.append(image);argv.extend(map(str,command));Path(log).parent.mkdir(parents=True,exist_ok=True)
    with Path(log).open('w') as w:
        w.write(json.dumps(argv)+'\n');w.flush();result=subprocess.run(argv,stdout=w,stderr=subprocess.STDOUT)
    if result.returncode:raise ContractError(f'Container failed exit={result.returncode}; see {log}')

def run_genes(args):
    prepared_contract(args.prepared,args.genome_id)
    out=Path(args.outdir).resolve();raw=out/'raw';raw.mkdir(parents=True,exist_ok=True)
    # A0 prepare is deterministic; sorting again also protects callers of this interface.
    dna=read_fasta(Path(args.prepared)/'genome.fna');write_fasta(raw/'training_input.fna',[(k,v[1]) for k,v in sorted(dna.items())])
    docker(['prodigal','-i','/output/training_input.fna','-p','single','-g','11','-f','gff','-o','/output/prodigal.gff','-a','/output/prodigal.faa','-d','/output/prodigal.ffn','-q'],CHECKM2_IMAGE,[(raw,'/output','rw')],out/'prodigal.log',1)
    args.raw=str(raw);genes_normalize(args)

def run_checkm2(args):
    prepared_contract(args.prepared,args.genome_id)
    db=Path(args.database).resolve()
    if not db.is_file() or db.stat().st_size!=3082500605 or file_hash(db)!='1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f':raise ContractError('Missing/invalid pinned full CheckM2 database')
    out=Path(args.outdir).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'raw/quality_report.tsv').exists():raise ContractError('Refuse overwrite of previous CheckM2 result; use a new output directory or Nextflow cache')
    docker(['checkm2','predict','--threads',str(args.threads),'--lowmem','--input','/input/genome.fna','--output-directory','/output/raw','--database_path','/database/'+db.name],CHECKM2_IMAGE,[(args.prepared,'/input','ro'),(out,'/output','rw'),(db.parent,'/database','ro')],out/'checkm2.log',args.threads)
    args.checkm2_report=str(out/'raw/quality_report.tsv');args.checkm2_name='genome';qc_normalize(args)

def parse_fastani(path):
    lines=Path(path).read_text().splitlines()
    if not lines:return None
    if len(lines)!=1:raise ContractError('Expected exactly one FastANI pair')
    cols=lines[0].split('\t')
    if len(cols)!=5:raise ContractError('Invalid FastANI result')
    ani=numeric(cols[2],'ANI');matched=int(cols[3]);total=int(cols[4])
    if total<=0 or matched<0 or matched>total:raise ContractError('Invalid FastANI fragment counts')
    return ani,matched/total

def taxonomy_decision(pairs):
    passing=[r for r in pairs if r['passes_thresholds']=='true']
    if not passing:return 'unresolved',None,'no_trusted_independent_reference_passing_ANI_and_fragment_thresholds'
    taxa=set(r['reference_species_taxid'] for r in passing)
    if len(taxa)>1:return 'ambiguous',None,'multiple_species_references_pass; limited_panel_cannot_resolve'
    passing.sort(key=lambda r:(-min(float(r['ani_forward_percent']),float(r['ani_reverse_percent'])),r['reference_id']))
    return 'reference_supported',passing[0],'supported_within_supplied_reference_panel; not_exhaustive_species_assignment'

def validate_reference_manifest(manifest):
    """Validate the existing catalog identity contract without reading sequences."""
    if not isinstance(manifest, dict):
        raise ContractError('Reference manifest must be an object')
    version = manifest.get('catalog_version')
    if not isinstance(version, str) or not version.strip():
        raise ContractError('Missing/invalid reference catalog_version')
    refs = manifest.get('references')
    if not isinstance(refs, list) or any(not isinstance(r, dict) for r in refs):
        raise ContractError('Reference manifest references must be a list of objects')
    unique_index(refs, 'reference_id', 'reference manifest')
    return manifest


def taxonomy(args):
    prepared_contract(args.prepared,args.genome_id)
    out=Path(args.outdir).resolve();raw=Path(args.raw).resolve() if getattr(args,'raw',None) else out/'raw';raw.mkdir(parents=True,exist_ok=True)
    query=Path(args.prepared).resolve()/'genome.fna';meta=read_metadata(args.metadata,args.genome_id)
    manifest=json.loads(Path(args.reference_manifest).read_text()) if args.reference_manifest else {'catalog_version':'none','references':[]}
    refs=validate_reference_manifest(manifest)['references'];pairs=[]
    if not(0<=args.ani_min<=100 and 0<=args.fragment_fraction_min<=1):raise ContractError('ANI percent/fraction policy out of range')
    query_content=sorted(sequence_hash(v[1]) for v in read_fasta(query).values())
    for ref in sorted(refs,key=lambda r:r['reference_id']):
        path=Path(ref['fasta']);path=Path(args.reference_dir)/path.name if getattr(args,'reference_dir',None) else path if path.is_absolute() else Path(args.reference_manifest).resolve().parent/path
        if file_hash(path)!=ref['fasta_sha256']:raise ContractError('Reference checksum mismatch '+ref['reference_id'])
        pair=dict(genome_id=args.genome_id,reference_id=ref['reference_id'],reference_species=ref['species'],reference_species_taxid=str(ref['species_taxid']),status='not_run',ani_forward_percent='unknown',ani_reverse_percent='unknown',query_fragment_fraction='unknown',reference_fragment_fraction='unknown',passes_thresholds='false',query_sha256=file_hash(query),reference_sha256=ref['fasta_sha256'])
        if query_content==sorted(sequence_hash(v[1]) for v in read_fasta(path).values()):pair['status']='excluded_self_sequence';pairs.append(pair);continue
        if not ref.get('trust_source'):pair['status']='excluded_untrusted_reference';pairs.append(pair);continue
        forward=raw/(ref['reference_id']+'.forward.tsv');reverse=raw/(ref['reference_id']+'.reverse.tsv')
        mounts=[(query.parent,'/query','ro'),(path.parent,'/reference','ro'),(raw,'/output','rw')]
        for q,r,result in [('/query/'+query.name,'/reference/'+path.name,forward),('/reference/'+path.name,'/query/'+query.name,reverse)]:
            if args.command=='taxonomy-normalize':continue
            docker(['fastANI','-q',q,'-r',r,'--threads',str(args.threads),'--fragLen','3000','--minFraction','0.2','-o','/output/'+result.name],FASTANI_IMAGE,mounts,result.with_suffix('.log'),args.threads)
        f=parse_fastani(forward);r=parse_fastani(reverse)
        if f and r:
            pair.update(status='completed',ani_forward_percent=f[0],ani_reverse_percent=r[0],query_fragment_fraction=f[1],reference_fragment_fraction=r[1],passes_thresholds=str(min(f[0],r[0])>=args.ani_min and min(f[1],r[1])>=args.fragment_fraction_min).lower())
        else:pair['status']='no_reportable_match'
        pairs.append(pair)
    status,best,reason=taxonomy_decision(pairs)
    row=dict(genome_id=args.genome_id,ncbi_declared_organism=meta.get('ncbi_declared_organism','unknown'),ncbi_declared_taxid=meta.get('ncbi_declared_taxid','unknown'),taxonomy_status=status,reference_supported_species=best['reference_species'] if best else 'unknown',reference_supported_species_taxid=best['reference_species_taxid'] if best else 'unknown',reference_id=best['reference_id'] if best else 'unknown',reason=reason,reference_catalog_version=manifest['catalog_version'],ani_min_percent=args.ani_min,fragment_fraction_min=args.fragment_fraction_min,taxonomy_policy_version='fastani_bidirectional_reference_support_v1')
    for k in ('ani_forward_percent','ani_reverse_percent','query_fragment_fraction','reference_fragment_fraction'):row[k]=best[k] if best else 'unknown'
    write_tsv(out/'taxonomy.tsv',TAX_FIELDS,[row]);write_tsv(out/'reference_ani.tsv',PAIR_FIELDS,pairs)
    write_json(out/'taxonomy_manifest.json',dict(**normalizer_runtime_provenance(FASTANI_IMAGE),reference_manifest=manifest,reference_manifest_sha256=file_hash(Path(args.reference_manifest)) if args.reference_manifest else None,query_sha256=file_hash(query),policy=row))

def parser():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    def basic(name,fn):
        q=sub.add_parser(name);q.set_defaults(fn=fn);q.add_argument('--outdir',required=True);return q
    for name,fn in [('genes',run_genes),('genes-normalize',genes_normalize),('qc-normalize',qc_normalize),('checkm2',run_checkm2),('taxonomy',taxonomy),('taxonomy-normalize',taxonomy)]:
        q=basic(name,fn);q.add_argument('--genome-id',required=True);q.add_argument('--prepared',required=True)
        if name=='genes-normalize':q.add_argument('--raw',required=True)
        if name in ('qc-normalize','checkm2','taxonomy','taxonomy-normalize'):q.add_argument('--metadata')
        if name in ('qc-normalize','checkm2'):
            q.add_argument('--checkm2-report');q.add_argument('--checkm2-name',default='genome')
            q.add_argument('--completeness-min',type=float,default=90);q.add_argument('--contamination-max',type=float,default=5)
            q.add_argument('--high-completeness-min',type=float,default=95);q.add_argument('--high-contamination-max',type=float,default=2)
            q.add_argument('--fragmentation-n50',type=int,default=10000)
        if name in ('checkm2','taxonomy','taxonomy-normalize'):q.add_argument('--threads',type=int,choices=range(1,7),default=6)
        if name=='checkm2':q.add_argument('--database',required=True)
        if name in ('taxonomy','taxonomy-normalize'):
            q.add_argument('--reference-dir');q.add_argument('--raw');q.add_argument('--reference-manifest');q.add_argument('--ani-min',type=float,default=95);q.add_argument('--fragment-fraction-min',type=float,default=0.65)
    q=basic('deduplicate',deduplicate);q.add_argument('--inputs',nargs='+',required=True)
    return p

def main():
    args=parser().parse_args()
    try:args.fn(args)
    except (ContractError,OSError,ValueError,KeyError,csv.Error) as exc:
        write_json(Path(args.outdir)/'failure.json',dict(status='failed',stage=args.command,error=str(exc),version=VERSION));print(str(exc),file=sys.stderr);return 2
    (Path(args.outdir)/'failure.json').unlink(missing_ok=True)
    return 0

if __name__=='__main__':sys.exit(main())
