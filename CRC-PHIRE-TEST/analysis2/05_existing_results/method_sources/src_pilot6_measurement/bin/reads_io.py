#!/usr/bin/env python3
"""Strict identities, FASTQ validation and resumable explicit-accession acquisition."""
import argparse,csv,gzip,hashlib,io,itertools,json,os,re,shutil,sqlite3,subprocess,sys,tempfile,urllib.parse,urllib.request,zipfile
from pathlib import Path

def write_if_changed(path,data):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    content=data.encode() if isinstance(data,str) else data
    if p.exists() and p.read_bytes()==content:return False
    with tempfile.NamedTemporaryFile(dir=p.parent,prefix=p.name+'.',suffix='.tmp',delete=False) as handle:
        handle.write(content);temp=Path(handle.name)
    temp.replace(p);return True
def dump(path,data):
    return write_if_changed(path,json.dumps(data,indent=2,sort_keys=True)+'\n')
def sha(path,algo='sha256'):
    h=hashlib.new(algo)
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def table(path):
    with open(path,newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError('Missing/duplicate TSV headers')
        rows=list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise ValueError('Malformed TSV row')
        return rows
def write_tsv(path,rows,fields):
    stream=io.StringIO(newline='')
    writer=csv.DictWriter(stream,fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return write_if_changed(path,stream.getvalue())
def ident(s,label):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,100}',s or ''):raise ValueError(f'Invalid {label}: {s!r}')
    return s
def fastq(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        while True:
            a=f.readline()
            if not a:break
            b,c,d=f.readline().rstrip('\r\n'),f.readline(),f.readline().rstrip('\r\n')
            if not a.startswith('@') or not c.startswith('+') or not b or len(b)!=len(d):raise ValueError(f'Malformed FASTQ: {path}')
            if re.search('[^ACGTNacgtn]',b) or any(not 33<=ord(x)<=126 for x in d):raise ValueError(f'Unsupported FASTQ alphabet/quality: {path}')
            yield a.rstrip('\r\n'),b,c.rstrip('\r\n'),d
def qname(record):
    name=record[0][1:].split()
    if not name or not name[0] or len(name[0])>254 or any(ord(c)<33 for c in name[0]):
        raise ValueError('Missing/invalid FASTQ query name')
    return re.sub(r'/[12]$','',name[0])
def fq_counts(r1,r2=None):
    n=b=0
    # Disk-backed uniqueness keeps large FASTQs out of RAM; duplicate names would
    # otherwise make any-mate filtering or fragment counting silently ambiguous.
    with tempfile.TemporaryDirectory(prefix='crc-fastq-names-') as temp:
        db=sqlite3.connect(str(Path(temp)/'names.sqlite'))
        db.execute('PRAGMA synchronous=OFF');db.execute('PRAGMA journal_mode=OFF')
        db.execute('CREATE TABLE names (name TEXT PRIMARY KEY)')
        pairs=itertools.zip_longest(fastq(r1),fastq(r2)) if r2 else ((a,None) for a in fastq(r1))
        try:
            for a,c in pairs:
                if a is None or (r2 and (c is None or qname(a)!=qname(c))):
                    raise ValueError('PE files are not synchronized')
                name=qname(a)
                try:db.execute('INSERT INTO names VALUES (?)',(name,))
                except sqlite3.IntegrityError as error:
                    raise ValueError('Duplicate FASTQ QNAME within a readset: '+name) from error
                n+=1;b+=len(a[1])+(len(c[1]) if c is not None else 0)
        finally:db.close()
    return {'fragments':n,'reads':n*(2 if r2 else 1),'bases':b}
def run(cmd,log):
    Path(log).parent.mkdir(parents=True,exist_ok=True)
    with open(log,'a') as f:
        f.write(json.dumps([str(x) for x in cmd])+'\n');f.flush()
        p=subprocess.run([str(x) for x in cmd],stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'Command exit {p.returncode}; see {log}')
def reserve(path,needed):
    if shutil.disk_usage(path).free < needed+5*1024**3:raise RuntimeError(f'Insufficient space: need {needed+5*1024**3} bytes including reserve at {path}')
def download(url,path,md5=None,size=None):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists() and (size is None or p.stat().st_size==int(size)) and (not md5 or sha(p,'md5')==md5):return {'url':url,'path':str(p.resolve()),'sha256':sha(p),'bytes':p.stat().st_size,'source_md5':md5}
    part=p.with_suffix(p.suffix+'.part');reserve(p.parent,max(0,(size or 10*1024**3)-(part.stat().st_size if part.exists() else 0)))
    run(['curl','--fail','--location','--retry','3','--continue-at','-','--output',part,url],p.with_suffix(p.suffix+'.download.log'))
    if size and part.stat().st_size!=int(size):raise ValueError(f'Download size mismatch: {url}')
    if md5 and sha(part,'md5')!=md5:raise ValueError(f'Download MD5 mismatch: {url}')
    part.replace(p);return {'url':url,'path':str(p.resolve()),'sha256':sha(p),'bytes':p.stat().st_size,'source_md5':md5}
def docker(tool,args,dirs):
    import spark_runtime
    if spark_runtime.enabled():
        spark_runtime.allocation()
        return spark_runtime.command(tool,args,writable=dirs)
    assets=json.loads((Path(__file__).resolve().parents[1]/'assets/reads_tools.json').read_text())
    cmd=['docker','run','--rm','--cpus','2','--memory','4g','--user',f'{os.getuid()}:{os.getgid()}']
    for d in sorted(set(str(Path(x).resolve()) for x in dirs)):cmd+=['-v',f'{d}:{d}']
    return cmd+[assets[tool]['container']]+list(map(str,args))
def ena_files(meta,layout):
    fields=[meta.get(k,'').split(';') if meta.get(k) else [] for k in ('fastq_ftp','fastq_md5','fastq_bytes')]
    urls,hashes,sizes=fields
    if not urls:return None
    if len({len(x) for x in fields})!=1:raise ValueError('ENA FASTQ/hash/size list lengths differ')
    entries=[]
    for u,m,size in zip(urls,hashes,sizes):
        if not re.fullmatch('[a-fA-F0-9]{32}',m) or not size.isdigit() or int(size)<=0:
            raise ValueError('Invalid ENA MD5/size metadata')
        url=u if u.startswith('https://') else 'https://'+u
        if urllib.parse.urlparse(url).scheme!='https':raise ValueError('ENA download must use HTTPS')
        name=Path(urllib.parse.urlparse(url).path).name
        if not name:raise ValueError('ENA file URL has no basename')
        entries.append(dict(url=url,md5=m.lower(),size=int(size),name=name))
    if len({e['name'] for e in entries})!=len(entries):raise ValueError('Duplicate ENA filename')
    if layout=='SE':return (entries,[]) if len(entries)==1 else None
    mates={}
    for e in entries:
        match=re.search(r'_([12])\.(?:fastq|fq)(?:\.gz)?$',e['name'])
        if match:
            if match[1] in mates:raise ValueError('Ambiguous ENA mate filenames')
            mates[match[1]]=e
    if set(mates)=={'1','2'} and len(entries) in (2,3):
        orphan=[e for e in entries if e not in mates.values()]
        return [mates['1'],mates['2']],orphan
    return None
def validate_ena_source(meta):
    source=(meta.get('library_source') or 'unknown').upper()
    if 'RNA' in source or 'TRANSCRIPT' in source:raise ValueError('ENA library_source indicates RNA/transcriptomic material outside DNA scope')
    return source

def acquire_reads(row,cache):
    acc=row['run_accession'];ident(acc,'run_accession')
    if not re.fullmatch(r'[SED]RR\d+',acc):raise ValueError('Require explicit ENA/SRA run accession')
    dest=Path(cache)/acc;dest.mkdir(parents=True,exist_ok=True)
    saved=dest/'acquisition_manifest.json'
    if saved.exists():
        old=json.loads(saved.read_text())
        validate_ena_source(old)
        if old['accession']!=acc or old['layout']!=row['layout']:raise ValueError('Cached accession/layout mismatch')
        if all(Path(f['path']).is_file() and sha(f['path'])==f['sha256'] for f in old['files']):
            files=old['primary'];return files[0],files[1] if len(files)==2 else None
    fields='run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,library_strategy,library_source,instrument_platform'
    url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode({'accession':acc,'result':'read_run','fields':fields,'format':'json'})
    try:
        data=json.loads(urllib.request.urlopen(url,timeout=60).read());meta=data[0] if len(data)==1 else {}
    except Exception as error:meta={'ena_error':str(error)}
    dump(dest/'ena_metadata.json',{'url':url,'response':meta})
    if meta.get('run_accession') and meta['run_accession']!=acc:raise ValueError('ENA returned a different run')
    if meta.get('instrument_platform') and meta['instrument_platform']!='ILLUMINA':raise ValueError('ENA platform contradicts Illumina contract')
    if meta.get('library_strategy') and meta['library_strategy'] not in ('WGS',):raise ValueError('ENA strategy outside DNA WGS scope')
    if meta.get('library_layout') and meta['library_layout']!=('PAIRED' if row['layout']=='PE' else 'SINGLE'):raise ValueError('ENA layout mismatch')
    library_source=validate_ena_source(meta)
    selection=ena_files(meta,row['layout'])
    receipts=[];singletons=[]
    if selection:
        primary,orphans=selection
        files=[]
        for item in primary+orphans:
            path=dest/item['name']
            receipts.append(download(item['url'],path,item['md5'],item['size']))
            (files if item in primary else singletons).append(str(path.resolve()))
    else:
        max_gib=int(row.get('sra_max_download_gib') or 30)
        if max_gib<=0:raise ValueError('SRA download cap must be positive')
        reserve(dest,max_gib*1024**3)
        dump(dest/'sra_space_preflight.json',{'phase':'before_prefetch','unknown_size_download_cap_gib':max_gib,'available_bytes':shutil.disk_usage(dest).free})
        run(docker('reads-fetch',['prefetch',acc,'--max-size',str(max_gib)+'G','-O',dest/acc],[dest]),dest/'sra.log')
        archive_files=list((dest/acc).rglob('*.sra'))+list((dest/acc).rglob('*.sralite'))
        archive_bytes=sum(f.stat().st_size for f in archive_files)
        if not archive_bytes:raise ValueError('Cannot determine downloaded SRA archive size for conversion preflight')
        available=shutil.disk_usage(dest).free
        required=17*archive_bytes+5*1024**3
        dump(dest/'sra_space_preflight.json',{'phase':'before_fasterq','archive_bytes':archive_bytes,
            'conversion_factor':17,'required_free_bytes_including_reserve':required,'available_bytes':available,
            'deficit_bytes':max(0,required-available),
            'basis':'https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump',
            'policy':'17x actual downloaded archive plus 5 GiB; fasterq internal disk checks remain enabled'})
        reserve(dest,17*archive_bytes)
        run(docker('reads-fetch',['fasterq-dump',dest/acc,'--split-3','--force','--threads','2','--mem','2GB','--temp',dest/'tmp','-O',dest],[dest]),dest/'sra.log')
        files=[str((dest/f'{acc}_{i}.fastq').resolve()) for i in (1,2)] if row['layout']=='PE' else [str((dest/f'{acc}.fastq').resolve())]
        orphan=dest/f'{acc}.fastq'
        if row['layout']=='PE' and orphan.is_file():singletons=[str(orphan.resolve())]
        if any(not Path(f).is_file() for f in files):raise ValueError('SRA output layout differs from declared layout; inspect download directory')
        receipts=[{'accession':acc,'source':'NCBI SRA fallback','path':f,'sha256':sha(f),'bytes':Path(f).stat().st_size} for f in files+singletons]
    fq_counts(files[0],files[1] if len(files)==2 else None)
    for f in singletons:fq_counts(f)
    dump(dest/'download_manifest.json',receipts)
    dump(saved,{'accession':acc,'layout':row['layout'],'primary':files,'singletons':singletons,'files':receipts,'library_source':library_source,'molecule_evidence':'ENA_genomic_source_and_WGS' if library_source in ('GENOMIC','METAGENOMIC') else 'declared_DNA_source_unknown','quality_encoding':'phred33'})
    return files[0],files[1] if len(files)==2 else None
def prepare(sheet,out,merge_map=None):
    sheet=Path(sheet).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    rows=table(sheet);seen=set();libs={};library_layouts={};samples={};units=[];accessions=set();used_paths=set();run_keys={}
    for row in rows:
        for k in ('biological_sample_id','library_id','readset_id'):ident(row.get(k),k)
        declared_encoding=bool(row.get('quality_encoding'))
        encoding=(row.get('quality_encoding') or 'phred33').lower()
        if encoding!='phred33':raise ValueError('Only declared/contract phred33 quality encoding is supported; no automatic inference')
        row['quality_encoding']=encoding
        row['quality_encoding_basis']='explicit_declaration' if declared_encoding else 'contract_phred33_not_inferred'
        row['human_remove']=row.get('human_remove','true').lower()
        if row['human_remove'] not in ('true','false'):raise ValueError('human_remove must be true/false')
        if row['human_remove']=='false' and not row.get('human_removal_reason'):
            raise ValueError('Disabled human removal requires an explicit reason')
        rid=row['readset_id']
        if rid in seen:raise ValueError(f'Duplicate readset_id {rid}')
        seen.add(rid)
        if row.get('material_type') not in ('bacterial_isolate','bulk_metagenome','vlp'):raise ValueError('Unsupported material_type')
        if row.get('molecule')!='DNA' or row.get('platform')!='ILLUMINA' or row.get('layout') not in ('SE','PE'):raise ValueError('Only DNA/ILLUMINA SE/PE supported')
        if row.get('genbank'):
            gb=Path(row['genbank']);gb=gb if gb.is_absolute() else sheet.parent/gb
            if not gb.is_file():raise ValueError('Provided GenBank input is missing: '+str(gb))
            row['genbank']=str(gb.resolve())
            row['genbank_sha256']=sha(gb)
        subject=row.get('subject_id') or None
        if subject is not None and subject.strip().lower() in ('unknown','na','null','none',''):subject=None
        key=(row['biological_sample_id'],row['material_type'],subject,row.get('mda','unknown'),row.get('human_remove','true'))
        library_layouts.setdefault(row['library_id'],set()).add(row['layout'])
        if row['library_id'] in libs and libs[row['library_id']]!=key:raise ValueError('Library foreign key/attributes conflict')
        libs[row['library_id']]=key
        run_id=row.get('run_id') or row.get('run_accession') or ('local:'+rid)
        if run_id.strip().lower() in ('unknown','na','none','null'):run_id='local:'+rid
        run_key=(row['library_id'],row['biological_sample_id'],row['platform'],row['molecule'],row['layout'])
        if run_id in run_keys and run_keys[run_id]!=run_key:raise ValueError('Run maps to conflicting libraries or samples')
        run_keys[run_id]=run_key;row['run_id']=run_id
        if row['biological_sample_id'] in samples and samples[row['biological_sample_id']]!=subject:raise ValueError('Sample maps to conflicting subjects')
        samples[row['biological_sample_id']]=subject
        source_singletons=[]
        if row.get('run_accession'):
            if row['run_accession'] in accessions:raise ValueError('Repeated run accession must not become independent readsets')
            accessions.add(row['run_accession'])
            if row.get('reads_1') or row.get('reads_2'):raise ValueError('Specify local files OR run accession')
            r1,r2=acquire_reads(row,out/'downloads')
            acquired=json.loads((out/'downloads'/row['run_accession']/'acquisition_manifest.json').read_text())
            source_singletons=acquired.get('singletons',[])
            row['public_library_source']=acquired.get('library_source','unknown')
            row['molecule_evidence']=acquired.get('molecule_evidence','declared_DNA_source_unknown')
        else:
            def local(k):
                if not row.get(k):return None
                p=Path(row[k]);p=p if p.is_absolute() else sheet.parent/p
                if not p.is_file():raise ValueError(f'Input missing {p}')
                return str(p.resolve())
            r1,r2=local('reads_1'),local('reads_2')
        if not r1 or (row['layout']=='PE')!=bool(r2):raise ValueError('Input layout/file mismatch')
        for path in [r1,r2]+source_singletons:
            if path and path in used_paths:raise ValueError('Same FASTQ path is reused by multiple readsets/mates')
            if path:used_paths.add(path)
        counts=fq_counts(r1,r2)
        if not counts['fragments']:raise ValueError('Raw input contains zero reads')
        u={**row,'unit_id':rid,'subject_id':subject,'reads_1':r1,'reads_2':r2,'source_readsets':[rid],'source_singletons':source_singletons,'source_singleton_checksums':{f:sha(f) for f in source_singletons},'counts':counts,'checksums':{f:sha(f) for f in [r1,r2] if f},'metadata':{k:v for k,v in row.items() if k not in ('reads_1','reads_2')}}
        units.append(u)
    if not units:raise ValueError('Empty samplesheet')
    readset_units=list(units)
    if merge_map:
        groups={};assigned=set();byid={u['unit_id']:u for u in units}
        for m in table(merge_map):
            ident(m.get('merge_id'),'merge_id');rid=m.get('readset_id')
            if rid not in byid or rid in assigned:raise ValueError('Invalid/duplicate merge-map readset')
            assigned.add(rid);groups.setdefault(m['merge_id'],[]).append(byid[rid])
        result=[u for u in units if u['unit_id'] not in assigned]
        for gid,parts in sorted(groups.items()):
            if gid in byid:raise ValueError('merge_id must not collide with any original readset ID')
            keys=[(p['library_id'],p['layout'],p['material_type'],p['molecule'],p['platform'],p['metadata'].get('adapter_r1'),p['metadata'].get('adapter_r2'),p['metadata'].get('genbank_sha256')) for p in parts]
            if len(set(keys))!=1:raise ValueError('Only compatible technical runs within one library can merge')
            u=dict(parts[0]);u['unit_id']=gid;u['source_singletons']=[f for p in parts for f in p.get('source_singletons',[])];u['source_singleton_checksums']={f:h for p in parts for f,h in p.get('source_singleton_checksums',{}).items()};u['source_readsets']=sorted(p['unit_id'] for p in parts);u['merge_source_checksums']={p['unit_id']:p['checksums'] for p in parts}
            for mate in ('reads_1','reads_2'):
                if not u[mate]:continue
                target=out/'merged'/f'{gid}.{mate}.fastq.gz';target.parent.mkdir(exist_ok=True)
                inputs=[{'readset_id':part['unit_id'],'source_sha256':part['checksums'][part[mate]]} for part in sorted(parts,key=lambda part:part['unit_id'])]
                receipt=target.with_suffix('.merge.json')
                previous=json.loads(receipt.read_text()) if receipt.exists() else {}
                source_manifest={'inputs':inputs,'namespace_policy':'readset_prefix_v1','gzip_mtime':0}
                if not (target.exists() and previous.get('source_manifest')==source_manifest and previous.get('output_sha256')==sha(target)):
                    with tempfile.NamedTemporaryFile(dir=target.parent,prefix=target.name+'.',delete=False) as binary:
                        temp=Path(binary.name)
                        with gzip.GzipFile(filename='',fileobj=binary,mode='wb',mtime=0) as compressed:
                            with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as output:
                                for part in sorted(parts,key=lambda part:part['unit_id']):
                                    for rec in fastq(part[mate]):
                                        output.write('@'+part['unit_id']+':'+rec[0][1:]+'\n'+'\n'.join(rec[1:])+'\n')
                    if target.exists() and sha(target)==sha(temp):temp.unlink()
                    else:temp.replace(target)
                    dump(receipt,{'source_manifest':source_manifest,'output_sha256':sha(target)})
                u[mate]=str(target)
            u['counts']=fq_counts(u['reads_1'],u['reads_2']);u['checksums']={u[k]:sha(u[k]) for k in ('reads_1','reads_2') if u[k]};result.append(u)
        units=result
    if len({u['unit_id'] for u in units})!=len(units):raise ValueError('Duplicate final unit IDs')
    for u in units:dump(out/'units'/f'{u["unit_id"]}.json',u)
    write_tsv(out/'readsets.tsv',[{**u,'run_id':u.get('run_id') or u.get('run_accession') or 'local:'+u['readset_id']} for u in readset_units],['readset_id','run_id','biological_sample_id','library_id','subject_id','material_type','layout','molecule','platform','reads_1','reads_2'])
    write_tsv(out/'libraries.tsv',[{'library_id':lid,'biological_sample_id':key[0],'material_type':key[1],'subject_id':key[2],'layout':';'.join(sorted(library_layouts[lid]))} for lid,key in sorted(libs.items())],['library_id','biological_sample_id','material_type','subject_id','layout'])
    write_tsv(out/'samples.tsv',[{'biological_sample_id':sid,'subject_id':subject,'subject_status':'known' if subject else 'unknown'} for sid,subject in sorted(samples.items())],['biological_sample_id','subject_id','subject_status'])
    write_tsv(out/'merge_map.tsv',[{'unit_id':u['unit_id'],'readset_id':rid} for u in units for rid in u['source_readsets']],['unit_id','readset_id'])
    manifest={'schema_version':1,'kind':'raw_reads_manifest','samplesheet_sha256':sha(sheet),'merge_map_sha256':sha(merge_map) if merge_map else None,'units':units}
    dump(out/'raw_reads_manifest.json',manifest);return manifest
def genomes(sheet,out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True);rows=[]
    for r in table(sheet):
        acc=r.get('assembly_accession','')
        if not re.fullmatch(r'GC[AF]_\d+\.\d+',acc):raise ValueError('Require explicit assembly accession.version')
        gid=ident(r.get('genome_id'),'genome_id');d=out/acc;d.mkdir(exist_ok=True);z=d/'ncbi_dataset.zip'
        if not z.exists():run(docker('reads-fetch',['datasets','download','genome','accession',acc,'--include','genome,gbff','--filename',z],[out]),d/'download.log')
        reserve(d,4*z.stat().st_size)
        with zipfile.ZipFile(z) as f:
            for n in f.namelist():
                if not (d/n).resolve().is_relative_to(d):raise ValueError('Unsafe archive member')
            f.extractall(d)
        fas=list((d/'ncbi_dataset/data'/acc).glob('*genomic.fna'));gb=list((d/'ncbi_dataset/data'/acc).glob('*.gbff'))
        if len(fas)!=1:raise ValueError('Downloaded accession/version not found exactly once')
        rows.append({'genome_id':gid,'fasta':str(fas[0]),'genbank':str(gb[0]) if len(gb)==1 else '', 'assembly_accession':acc,'sha256':sha(fas[0]),'source_archive_sha256':sha(z)})
    write_tsv(out/'genomes.tsv',rows,['genome_id','fasta','genbank','assembly_accession','sha256','source_archive_sha256']);dump(out/'source_manifest.json',rows)
def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('prepare');a.add_argument('--samplesheet',required=True);a.add_argument('--outdir',required=True);a.add_argument('--merge-map')
    a=sub.add_parser('genomes');a.add_argument('--samplesheet',required=True);a.add_argument('--outdir',required=True)
    a=p.parse_args()
    if a.cmd=='prepare':prepare(a.samplesheet,a.outdir,a.merge_map)
    else:genomes(a.samplesheet,a.outdir)
if __name__=='__main__':
    try:main()
    except Exception as e:print(f'ERROR: {e}',file=sys.stderr);sys.exit(1)
