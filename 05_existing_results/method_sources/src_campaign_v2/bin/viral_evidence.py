#!/usr/bin/env python3
"""Optional V1 evidence; never changes catalog, vOTUs, B1 counts or primary taxonomy."""
import argparse, csv, hashlib, json, math, re, subprocess
from pathlib import Path

def sha(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def table(path, required=(), delimiter='\t'):
    with Path(path).open(newline='') as f:
        reader=csv.DictReader(f,delimiter=delimiter)
        fields=reader.fieldnames
        if not fields or len(fields)!=len(set(fields)) or not set(required)<=set(fields):
            raise ValueError('Invalid table header: '+str(path))
        rows=list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError('Invalid row width: '+str(path))
        return rows

def write(path, rows, fields):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def fasta(path):
    records={}; name=None; parts=[]
    def add():
        if name is not None:
            if name in records: raise ValueError('Duplicate FASTA identifier')
            records[name]=''.join(parts).upper()
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                add(); name=line[1:].split()[0]; parts=[]
            elif line.strip():
                if name is None: raise ValueError('Sequence before header')
                parts.append(line.strip())
    add(); return records

def load_catalog(root):
    rows=table(root/'viral_sequence_master.tsv',('viral_sequence_id','sequence_sha256','length','virus_scope','checkv_quality'))
    seqs=fasta(root/'viral_sequences.fna')
    ids=[r['viral_sequence_id'] for r in rows]
    if len(ids)!=len(set(ids)) or set(ids)!=set(seqs): raise ValueError('Catalog/FASTA IDs disagree')
    for row in rows:
        uid=row['viral_sequence_id']; seq=seqs[uid]
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',uid): raise ValueError('Unsafe catalog ID')
        if len(seq)!=int(row['length']) or hashlib.sha256(seq.encode()).hexdigest()!=row['sequence_sha256']:
            raise ValueError('Catalog sequence identity mismatch')
    return rows,seqs

def eligibility(row,tool):
    if tool=='bacphlip':
        if row['virus_scope']!='bacteriophage': return 'not_assessed_scope'
        if row['checkv_quality']!='Complete': return 'not_assessed_incomplete'
    elif row['virus_scope'] not in ('bacteriophage','archaeal_virus'): return 'not_assessed_scope'
    return 'eligible'

def parse_bacphlip(path,expected):
    rows=table(path,('','Virulent','Temperate'))
    out=[]; seen=set()
    for r in rows:
        uid=r['']
        if uid not in expected or uid in seen: raise ValueError('BACPHLIP unexpected/duplicate ID')
        seen.add(uid); v=float(r['Virulent']); t=float(r['Temperate'])
        if not all(math.isfinite(x) and 0<=x<=1 for x in (v,t)) or abs(v+t-1)>1e-6: raise ValueError('Invalid probabilities')
        out.append(dict(viral_sequence_id=uid,virulent_probability=v,temperate_probability=t,
                        lifestyle='unknown' if v==t else ('temperate' if t>v else 'virulent')))
    if seen!=set(expected): raise ValueError('Missing BACPHLIP predictions')
    return out

def parse_vcontact(path,expected):
    rows=table(path,('Genome','Reference'),','); out=[]; seen=set()
    for row in rows:
        uid=row['Genome']
        if uid not in expected:
            if row['Reference'].lower() not in ('true','1'): raise ValueError('Unexpected non-reference query classification')
            continue  # verified database records are retained only in raw output
        if uid in seen: raise ValueError('Duplicate query classification')
        seen.add(uid)
        if row['Reference'].lower() not in ('false','0',''): raise ValueError('Query classified as database reference')
        tax={k:v for k,v in row.items() if k.endswith(('_prediction','_reference'))}
        if not tax: raise ValueError('Missing version-specific taxonomy columns')
        out.append(dict(viral_sequence_id=uid,taxonomy_json=json.dumps(tax,sort_keys=True),
                        genus_prediction=row.get('genus_prediction',''),status='completed'))
    if seen!=set(expected): raise ValueError('Missing query classifications')
    return out

def unique_file(root,pattern):
    hits=list(root.rglob(pattern))
    if len(hits)!=1: raise ValueError('Expected one '+pattern+', found '+str(len(hits)))
    return hits[0]

def parse_vibrant(root,expected,seqs,tool_genes=None):
    quality=table(unique_file(root,'VIBRANT_genome_quality_*.tsv'),('scaffold','type','Quality'))
    outputseq=fasta(unique_file(root,'*.phages_combined.fna'))
    if {r['scaffold'] for r in quality}!=set(outputseq): raise ValueError('VIBRANT quality/sequence IDs disagree')
    mapped={}; out=[]
    for r in quality:
        local=r['scaffold']
        if local in mapped: raise ValueError('Duplicate VIBRANT scaffold')
        parents=[uid for uid in expected if local==uid or re.fullmatch(re.escape(uid)+r'_fragment_\d+',local)]
        if len(parents)!=1: raise ValueError('Cannot map VIBRANT scaffold to catalog')
        uid=parents[0]; seq=outputseq[local]; start=seqs[uid].find(seq)
        if not seq or start<0: raise ValueError('VIBRANT slice missing')
        if seqs[uid].find(seq,start+1)>=0:
            hints=[g for g in (tool_genes or []) if g['tool_scaffold_id']==local]
            # 1.2.1 fragment export is query[start-1:stop-1]; FAA/FFN retain query coordinates.
            left=min((int(g['start0']) for g in hints),default=-1)
            right=max((int(g['end0']) for g in hints),default=-1)-1
            if local==uid or left<0 or seqs[uid][left:right]!=seq: raise ValueError('VIBRANT slice ambiguous without coordinate evidence')
            start=left
        mapped[local]=uid
        out.append(dict(viral_sequence_id=uid,tool_scaffold_id=local,start0=start,end0=start+len(seq),
                        lifestyle=r['type'],tool_quality=r['Quality'],sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
    amgfile=unique_file(root,'VIBRANT_AMG_individuals_*.tsv')
    # VIBRANT 1.2.1 exports six fixed columns; keep source descriptions verbatim.
    with amgfile.open(newline='') as f:
        reader=csv.reader(f,delimiter='\t'); header=next(reader,None); raw=list(reader)
    if header!=['protein','scaffold','AMG KO','AMG KO name','Pfam','Pfam name']: raise ValueError('Unexpected VIBRANT AMG header')
    annotation=unique_file(root,'VIBRANT_annotations_*.tsv')
    mother_rows=table(annotation,('protein','scaffold'))
    mother={row['protein'] for row in mother_rows}
    if len(mother)!=len(mother_rows): raise ValueError('Duplicate mother gene')
    parent_by_gene={row['protein']:row['scaffold'] for row in mother_rows}
    entries=[]; seen=set()
    for row in raw:
        if len(row)!=6: raise ValueError('Invalid VIBRANT AMG row')
        # Actual header ordering is verified below rather than assumed silently.
        values=dict(zip(header,row))
        gene=values.get('protein',values.get('gene',''))
        scaffold=values.get('scaffold','')
        if not gene or gene not in mother or scaffold not in mapped or parent_by_gene[gene]!=scaffold: raise ValueError('AMG/mother/scaffold join failed')
        signature=tuple(row)
        if signature in seen: raise ValueError('Duplicate AMG entry')
        seen.add(signature)
        entries.append(dict(viral_sequence_id=mapped[scaffold],tool_scaffold_id=scaffold,tool_gene_id=gene,
                            annotation_json=json.dumps(values,sort_keys=True)))
    return out,entries

def verify_assets(config):
    for key in ('model_manifest','database_manifest'):
        if not config.get(key): continue
        path=Path(config[key])
        if sha(path)!=config[key+'_sha256']: raise ValueError('Asset manifest changed: '+key)
        doc=json.loads(path.read_text()); root=Path(doc.get('path',str(path.parent)))
        for item in doc['files']:
            p=(root/item['path']).resolve()
            if not p.is_relative_to(root.resolve()) or not p.is_file() or p.stat().st_size!=item['bytes']:
                raise ValueError('Missing/changed resource '+str(p))
            if key=='model_manifest' and sha(p)!=item['sha256']: raise ValueError('Changed model: '+str(p))
        # File hashes are validated at install time; manifest identity is bound to each run.

def run(args):
    catalog=Path(args.catalog).resolve(); out=Path(args.outdir).resolve()
    if out==catalog or out.is_relative_to(catalog) or catalog.is_relative_to(out): raise ValueError('Output overlaps input')
    if out.exists(): raise ValueError('Require new output directory')
    for extra in [getattr(args,'genes_dir',None),*getattr(args,'canonical_genes',[])]:
        if extra:
            p=Path(extra).resolve()
            if out==p or out.is_relative_to(p) or p.is_relative_to(out): raise ValueError('Output overlaps gene input')
    contract=json.loads(Path(args.contract).read_text()); config=contract['tools'][args.tool]
    if config['state']!='prepared': raise ValueError('Tool is not prepared: '+args.tool)
    verify_assets(config)
    rows,seqs=load_catalog(catalog); selected={r['viral_sequence_id'] for r in rows if eligibility(r,args.tool)=='eligible'}
    out.mkdir(parents=True); raw=out/'raw'; raw.mkdir(); inp=raw/'query.fna'
    inp.write_text(''.join('>'+uid+'\n'+seqs[uid]+'\n' for uid in sorted(selected)))
    statuses=[dict(viral_sequence_id=r['viral_sequence_id'],status=eligibility(r,args.tool),lifestyle='unknown') for r in rows]
    receipt=dict(tool=args.tool,tool_version=config['version'],contract_sha256=sha(args.contract),
        catalog_manifest_sha256=sha(catalog/'catalog_manifest.json'),selected_sequences=len(selected),
        status='empty_input' if not rows else ('not_assessed' if not selected else 'running'),scientifically_calibrated=False)
    def save(): (out/'status.json').write_text(json.dumps(receipt,indent=2)+'\n')
    save()
    try:
        if args.tool=='vibrant':
            import evidence_crosswalk as ec
            if not getattr(args,'genes_dir',None): raise ValueError('VIBRANT requires --genes-dir')
        if selected:
            if args.tool=='bacphlip': cmd=['bacphlip','-i','query.fna','--multi_fasta']
            elif args.tool=='vibrant': cmd=['VIBRANT_run.py','-i','query.fna','-folder','vibrant','-t',str(args.cpus),'-virome','-no_plot','-d',config['database_path'],'-m',config['model_path']]
            else: cmd=['vcontact3','run','--nucleotide','query.fna','--output','vcontact','--db-path',config['database_path'],'--db-version',str(config['database_version']),'--db-domain','prokaryotes','--threads',str(args.cpus),'--no-progress']
            receipt['command']=cmd; save()
            with (out/'tool.log').open('w') as log:
                rc=subprocess.run(cmd,cwd=raw,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout).returncode
            receipt['exit_code']=rc
            if rc: raise RuntimeError('Tool returned '+str(rc))
            if args.tool=='bacphlip':
                result=parse_bacphlip(raw/'query.fna.bacphlip',selected)
                write(out/'lifestyle.tsv',result,['viral_sequence_id','virulent_probability','temperate_probability','lifestyle'])
            elif args.tool=='vcontact3':
                result=parse_vcontact(unique_file(raw,'final_assignments.csv'),selected)
                write(out/'taxonomy.tsv',result,['viral_sequence_id','taxonomy_json','genus_prediction','status'])
            else:
                tool_genes=ec.parse_tool_genes(raw,selected,seqs)
                result,entries=parse_vibrant(raw,selected,seqs,tool_genes)
                write(out/'viral_regions.tsv',result,['viral_sequence_id','tool_scaffold_id','start0','end0','lifestyle','tool_quality','sequence_sha256'])
                write(out/'amg_entries.tsv',entries,['viral_sequence_id','tool_scaffold_id','tool_gene_id','annotation_json'])
                receipt['canonical_gene_crosswalk']=ec.export(raw,catalog,args.genes_dir,args.canonical_genes,tool_genes,result,entries,seqs,out)
            hits={r['viral_sequence_id'] for r in result}
            for row in statuses:
                if row['viral_sequence_id'] in selected: row['status']='completed' if row['viral_sequence_id'] in hits else 'completed_no_hit'
            receipt['status']='completed'
        elif args.tool=='vibrant':
            write(out/'viral_regions.tsv',[],['viral_sequence_id','tool_scaffold_id','start0','end0','lifestyle','tool_quality','sequence_sha256'])
            write(out/'amg_entries.tsv',[],['viral_sequence_id','tool_scaffold_id','tool_gene_id','annotation_json'])
            receipt['canonical_gene_crosswalk']=ec.export(raw,catalog,args.genes_dir,args.canonical_genes,[],[],[],seqs,out)
        write(out/'sequence_status.tsv',statuses,['viral_sequence_id','status','lifestyle'])
        receipt['raw_files']=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=sha(p)) for p in raw.rglob('*') if p.is_file()]
        (out/'OUTPUTS.md').write_text('Optional evidence only. Join by viral_sequence_id; do not sum abundance after one-to-many joins.\n'
            'sequence_status.tsv reports eligibility/execution; its lifestyle is unknown unless read from the tool-specific table.\n'
            'BACPHLIP: complete bacteriophages only; probabilities are model outputs, not calibrated accuracy.\n'
            'VIBRANT: all coordinates are 0-based half-open; tool ORFs remain independent IDs. gene_crosswalk.tsv maps to viral and host_canonical namespaces.\n'
            'Only unambiguous identical coordinates/CDS/protein yield exact; partial/ambiguous/unmapped are retained, without annotation transfer.\n'
            'VIBRANT 1.2.1 fragment FASTA may omit the terminal base; ORF coordinates retain the full original-query CDS extent.\n'
            'AMG entry counts and unique tool genes are separate. Candidate AMG is not experimental metabolic evidence.\n'
            'vConTACT3: query records only in taxonomy.tsv; references remain in raw. Predictions do not overwrite geNomad; vOTU is not VC.\n'
            'Unknown, not assessed, no hit and failure are distinct. Missing files cause failure, never an assumed negative.\n')
        receipt['output_files']={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p!=out/'status.json'}
        save()
    except Exception as exc:
        receipt['status']='failed'; receipt['error']=str(exc); save(); raise

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tool',choices=['bacphlip','vibrant','vcontact3'],required=True)
    p.add_argument('--catalog',required=True);p.add_argument('--contract',required=True);p.add_argument('--outdir',required=True)
    p.add_argument('--genes-dir');p.add_argument('--canonical-genes',nargs='*',default=[])
    p.add_argument('--cpus',type=int,default=1);p.add_argument('--timeout',type=int,default=21600)
    a=p.parse_args()
    if not 1<=a.cpus<=6: p.error('cpus must be 1..6')
    run(a)
