"""Bounded sequence-context and read-junction review; no induction claims.

Uses existing coordinate-sorted assembly readback BAMs. New indexes live only
inside the extension output, never beside archived input symlinks.
"""
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime, timezone
import argparse, csv, hashlib, json, math, os, re, subprocess

def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def table(p):
    with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
    with Path(p).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def checked_manifest(p,synthetic_test=False):
    d=read(p);assert d['status']=='completed',str(p)
    assert bool(d.get('synthetic_preview',False))==synthetic_test,'Synthetic/production evidence mismatch'
    assert d['outputs'],'Missing output inventory'
    for rel,h in d['outputs'].items():
        output=(p.parent/rel).resolve();assert output.is_relative_to(p.parent.resolve()),rel
        assert sha(output)==h,rel
    for path,h in d.get('input_sha256',{}).items():assert sha(Path(path))==h,path
    return d
def fasta_selected(path,ids):
    found={};name=None;parts=[]
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name in ids:found[name]=''.join(parts).upper()
                name=line[1:].split()[0];parts=[]
            elif name in ids:parts.append(line.strip())
    if name in ids:found[name]=''.join(parts).upper()
    assert set(found)==set(ids),(str(path),set(ids)-set(found))
    return found
def entropy(seq):
    return -sum((seq.count(b)/len(seq))*math.log2(seq.count(b)/len(seq)) for b in 'ACGT' if seq.count(b)) if seq else 0
def boundary_repeat(seq,start,end):
    """Exploratory same-orientation repeat within +/-250 bp, NOT att calling."""
    left0=max(0,start-250);left1=min(len(seq),start+250)
    right0=max(0,end-250);right1=min(len(seq),end+250)
    if left1>right0:return dict(repeat_status='not_assessed_overlapping_windows')
    left=seq[left0:left1];right=seq[right0:right1];best=None
    for n in range(60,14,-1):
        seen={}
        for i in range(len(left)-n+1):
            x=left[i:i+n]
            if set(x)<=set('ACGT') and entropy(x)>=1.5:seen.setdefault(x,[]).append(i)
        hits=[]
        for j in range(len(right)-n+1):
            x=right[j:j+n]
            if x in seen:
                for i in seen[x]:
                    if left0+i+n<=right0+j:hits.append((abs(left0+i-start)+abs(right0+j-end),left0+i,right0+j,x))
        if hits:
            best=min(hits);break
    if best:
        _,a,b,x=best
        return dict(repeat_status='direct_repeat_candidate_not_confirmed_att',repeat_length=len(x),repeat_sequence=x,repeat_left_start0=a,repeat_right_start0=b)
    return dict(repeat_status='no_repeat_under_this_heuristic_not_prophage_absence')
def spanning_read(fields,boundary,anchor=20,minq=20):
    """One uninterrupted aligned block with high-quality anchors on both sides."""
    if len(fields)<11:return False
    flag=int(fields[1]);mapq=int(fields[4])
    if flag & (4|256|512|1024|2048) or mapq<20:return False
    cigar=fields[5];qual=fields[10]
    if cigar=='*' or qual=='*':return False
    ops=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in ops)!=cigar:return False
    ref=int(fields[3])-1;query=0
    for n,op in ops:
        n=int(n)
        if op in 'M=X':
            if ref<=boundary-anchor and ref+n>=boundary+anchor:
                a=query+boundary-anchor-ref;b=query+boundary+anchor-ref
                if 0<=a<b<=len(qual) and min(ord(c)-33 for c in qual[a:b])>=minq:return True
            ref+=n;query+=n
        elif op in 'DN':ref+=n
        elif op in 'IS':query+=n
    return False
def choose(occ,diff,links,limit=20):
    broad={r['m_sequence_id'] for r in links if r['match']!='partial_or_fragment_match'}
    eligible=[r for r in occ if r['source_kind']!='legacy_bacterial_genome' and r['candidate_type']=='provirus_locus' and r.get('votu_id') in diff]
    def key(r):
        d=diff[r['votu_id']];pooled=int(d.get('crc_detected') or 0)+int(d.get('control_detected') or 0)
        return (r['viral_sequence_id'] not in broad,r.get('checkv_quality') not in {'Complete','High-quality'},-pooled,-int(r['length']),r['source_occurrence_id'])
    return sorted(eligible,key=key)[:limit]
def review(root,limit=20,samtools=None,*,_synthetic_test=False):
    assert os.environ.get('SLURM_JOB_ID'),'Production context review runs through Slurm'
    root=Path(root).resolve();assert 0<limit<=20
    if _synthetic_test:assert (root/'SYNTHETIC_TEST_DATA_ONLY.txt').is_file()
    report=root/'12_reports/pilot10';evidence=report/'prophage_evidence';stats=report/'exploratory'
    er=checked_manifest(evidence/'evidence_manifest.json',_synthetic_test)
    sr=checked_manifest(stats/'exploratory_manifest.json',_synthetic_test)
    assert 'occurrence_evidence.tsv' in er['outputs'] and 'votu_differential.tsv' in sr['outputs']
    assert sr['input_sha256'][str(evidence/'evidence_manifest.json')]==sha(evidence/'evidence_manifest.json')
    out=report/'candidate_review';out.mkdir(exist_ok=True)
    receipt=out/'candidate_review_manifest.json'
    if receipt.exists():
        previous=checked_manifest(receipt,_synthetic_test)
        assert previous['script_sha256']==sha(Path(__file__)) and previous['maximum_occurrences']==limit
        for bam in previous['bams']:assert sha(bam['path'])==bam['sha256'],bam['path']
        print(json.dumps(dict(status='reused_verified_completed',selected_occurrences=previous['selected_occurrences'])))
        return
    assert not any(out.iterdir()),'Inspect preserved partial output before recovery; do not overwrite'
    indexdir=out/'bam_indexes';indexdir.mkdir()
    occ=table(evidence/'occurrence_evidence.tsv');diff={r['votu_id']:r for r in table(stats/'votu_differential.tsv')}
    links=table(report/'G_M_sequence_matches.tsv');selected=choose(occ,diff,links,limit)
    frozen_path=root/'01_manifests/pilot10_frozen_assemblies.json';frozen=read(frozen_path)
    assert sr['input_sha256'][str(report/'G_M_sequence_matches.tsv')]==sha(report/'G_M_sequence_matches.tsv')
    assert sr['input_sha256'][str(frozen_path)]==sha(frozen_path)
    units={u['biological_sample_id']:u for u in frozen['units']}
    assert len(units)==len(frozen['units'])==10 and set(units)==set(frozen['groups'])
    assert Counter(frozen['groups'].values())=={'CRC':5,'control':5}
    assert len({r['source_occurrence_id'] for r in occ})==len(occ)
    assert len(diff)==len(table(stats/'votu_differential.tsv')),'Duplicate vOTU statistics'
    inputs={str(p):sha(p) for p in [evidence/'evidence_manifest.json',stats/'exploratory_manifest.json',report/'G_M_sequence_matches.tsv',frozen_path]}
    cat=root/'06_viral_catalog/pilot10_combined'
    gene_path=cat/'genes/gene_table.tsv';hit_path=cat/'functions/function_hits.tsv'
    terminal=read(cat/'pipeline_info/terminal_output_manifest.json')
    inputs[str(cat/'pipeline_info/terminal_output_manifest.json')]=sha(cat/'pipeline_info/terminal_output_manifest.json')
    for rel,key in [('genes/gene_calling_status.json','gene_table.tsv'),('functions/viral_function_status.json','function_hits.tsv')]:
        receipt=cat/rel;assert sha(receipt)==terminal[rel]
        obj=read(receipt);assert obj['status']=='completed';path=receipt.parent/key
        assert sha(path)==obj['output_files'][key];inputs[str(path)]=sha(path);inputs[str(receipt)]=sha(receipt)
    genes=defaultdict(list)
    for row in table(gene_path):genes[row['viral_sequence_id']].append(row)
    hits={r['gene_id']:r for r in table(hit_path) if r.get('accepted')=='true' and r.get('is_best')=='true'}
    function_receipt=read(cat/'functions/viral_function_status.json')
    system_path=cat/'functions/systems.tsv';system_hits_path=cat/'functions/system_hits.tsv'
    for path in [system_path,system_hits_path]:
        assert sha(path)==function_receipt['output_files'][path.name];inputs[str(path)]=sha(path)
    systems={r['system_id']:r for r in table(system_path)};gene_systems=defaultdict(set)
    for h in table(system_hits_path):
        s=systems[h['system_id']];gene_systems[h['gene_id']].add(s['activity']+':'+s['subtype'])
    ko_path=cat/'viral_kofam/kofam/viral_member_ko_counts.tsv';ko_receipt=ko_path.parent/'viral_ko_manifest.json'
    ko=read(ko_receipt);assert ko['status']=='completed' and sha(ko_path)==ko['outputs'][ko_path.name]
    assert ko['genes_sha256']==sha(gene_path);inputs[str(ko_path)]=sha(ko_path);inputs[str(ko_receipt)]=sha(ko_receipt)
    membership=cat/'catalog/vOTU_members.tsv'
    assert sha(membership)==terminal['catalog/vOTU_members.tsv'] and ko['membership_sha256']==sha(membership)
    inputs[str(membership)]=sha(membership)
    for rel,h in ko['outputs'].items():assert sha(ko_receipt.parent/rel)==h,rel
    gene_kos=defaultdict(set)
    for k in table(ko_path):
        for gid in k['gene_ids'].split(';'):
            if gid:gene_kos[gid].add(k['ko'])
    defaults=root/'02_pipeline/src_campaign_v2/assets/spark_dependencies.json'
    if samtools is None:
        samtools=Path(read(defaults)['environments']['reads-assembly']['prefix'])/'bin/samtools'
        inputs[str(defaults)]=sha(defaults)
    version=subprocess.run([str(samtools),'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    bydonor=defaultdict(list)
    for row in selected:bydonor[row['biological_sample_id']].append(row)
    junction=[];repeats=[];gene_rows=[];bam_audit=[];commands=[]
    for sid,rows in sorted(bydonor.items()):
        u=units[sid];fasta=Path(u['assembly_fasta']);assert sha(fasta)==u['sha256']
        inputs[str(fasta)]=u['sha256'];seqs=fasta_selected(fasta,{r['original_contig_id'] for r in rows})
        bam=fasta.with_name('readback.bam');available=bam.is_file();bam_reason='missing_existing_readback_bam'
        index=indexdir/(sid+'.bai')
        if available:
            receipt_path=root/f'13_handoff/archive_{sid}.json';archive=read(receipt_path)
            assert archive['status']=='completed' and archive['sample_id']==sid
            archived=[r for r in archive['files'] if r['path']==str(bam)]
            assert len(archived)==1 and archived[0]['status']=='archived_verified_link',str(bam)
            assert bam.stat().st_size==archived[0]['bytes'];digest=sha(bam)
            assert digest==archived[0]['sha256'],'Changed archived BAM'
            inputs[str(receipt_path)]=sha(receipt_path)
            cmd=[str(samtools),'index','-@','3','-o',str(index),str(bam)]
            commands.append(cmd);subprocess.run(cmd,check=True)
            bam_audit.append(dict(sample_id=sid,path=str(bam),resolved_path=str(bam.resolve()),sha256=digest,index=str(index),status='verified_existing_bam_new_external_index'))
        for row in rows:
            oid=row['source_occurrence_id'];contig=row['original_contig_id'];seq=seqs[contig]
            assert row['strand']=='+','Current primary catalog normalizes to source orientation; unsupported orientation'
            a,b=int(row['start0']),int(row['end0']);assert 0<=a<b<=len(seq)
            assert hashlib.sha256(seq[a:b].encode()).hexdigest()==row['sequence_sha256']
            repeats.append(dict(source_occurrence_id=oid,viral_sequence_id=row['viral_sequence_id'],**boundary_repeat(seq,a,b),interpretation='sequence_repeat_screen_only_no_boundary_change_no_att_confirmation'))
            for side,boundary in [('left',a),('right',b)]:
                base=dict(source_occurrence_id=oid,biological_sample_id=sid,original_contig_id=contig,side=side,boundary0=boundary,anchor_bp=20,min_mapq=20,min_baseq=20)
                if boundary<20 or len(seq)-boundary<20:
                    junction.append(dict(base,status='not_assessable_contig_end',spanning_fragment_count=''));continue
                if not available:
                    junction.append(dict(base,status=bam_reason,spanning_fragment_count=''));continue
                region=f'{contig}:{boundary-19}-{boundary+20}'
                cmd=[str(samtools),'view','-X','-q','20','-F','3844',str(bam),str(index),region];commands.append(cmd)
                result=subprocess.run(cmd,capture_output=True,text=True,check=True)
                names=set();records=0
                for line in result.stdout.splitlines():
                    fields=line.split('\t');records+=1
                    if spanning_read(fields,boundary):names.add(fields[0])
                junction.append(dict(base,status='assessed_assembly_junction_not_excision',spanning_fragment_count=len(names),overlapping_filtered_alignment_records=records))
            for g in genes[row['viral_sequence_id']]:
                h=hits.get(g['gene_id'],{});ga,gb=int(g['start0']),int(g['end0']);assert 0<=ga<gb<=b-a
                desc=h.get('annotation','');category=h.get('category','');distance=min(ga,b-a-gb)
                gene_rows.append(dict(source_occurrence_id=oid,viral_sequence_id=row['viral_sequence_id'],gene_id=g['gene_id'],viral_start0=ga,viral_end0=gb,source_start0=a+ga,source_end0=a+gb,strand=g['strand'],partial=g['partial'],phrog=h.get('phrog',''),annotation=desc,category=category,defense_antidefense_systems=';'.join(sorted(gene_systems[g['gene_id']])),ko_homology=';'.join(sorted(gene_kos[g['gene_id']])),boundary_distance_bp=distance,integration_keyword_hint=bool(re.search(r'\b(integrase|excisionase|recombinase)\b',desc,re.I)),interpretation='own_sequence_homology_and_coordinates_only_not_activity_or_AMG_confirmation'))
    base_fields=list(selected[0]) if selected else ['source_occurrence_id','viral_sequence_id','votu_id','biological_sample_id']
    write(out/'selected_occurrences.tsv',selected,base_fields)
    write(out/'junction_read_support.tsv',junction,['source_occurrence_id','biological_sample_id','original_contig_id','side','boundary0','anchor_bp','min_mapq','min_baseq','status','spanning_fragment_count','overlapping_filtered_alignment_records'])
    write(out/'boundary_repeat_candidates.tsv',repeats,['source_occurrence_id','viral_sequence_id','repeat_status','repeat_length','repeat_sequence','repeat_left_start0','repeat_right_start0','interpretation'])
    write(out/'candidate_gene_context.tsv',gene_rows,['source_occurrence_id','viral_sequence_id','gene_id','viral_start0','viral_end0','source_start0','source_end0','strand','partial','phrog','annotation','category','defense_antidefense_systems','ko_homology','boundary_distance_bp','integration_keyword_hint','interpretation'])
    (out/'commands.json').write_text(json.dumps(commands,indent=2),encoding='utf-8')
    d=dict(status='completed',synthetic_preview=_synthetic_test,utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ['SLURM_JOB_ID'],selected_occurrences=len(selected),maximum_occurrences=limit,
        selection='M predicted provirus with mapped vOTU; broad G-M match first, then CheckV Complete/High, pooled detected donor count, length, occurrence ID; no P/effect selection',
        boundary_repeat_rule='Exact same-orientation 15-60bp repeat in +/-250bp windows, Shannon entropy>=1.5; heuristic candidate only, not confirmed att; does not alter boundaries',
        junction_rule='Existing same-sample assembly BAM; MAPQ>=20, exclude unmapped/secondary/supplementary/QCfail/duplicate; a continuous CIGAR block anchors >=20bp each side with baseQ>=20; fragment names deduplicated per boundary',
        interpretation='Supports local assembly continuity at proposed boundary, not exact biological attachment site, patient host species, excision or inducibility. No reads is not evidence of biological absence. Same reads were used to assemble.',
        samtools_version=version,bams=bam_audit,input_sha256=inputs,script_sha256=sha(Path(__file__)),outputs={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})
    (out/'candidate_review_manifest.json').write_text(json.dumps(d,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='completed',selected_occurrences=len(selected),bams=len(bam_audit))))
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--limit',type=int,default=20);p.add_argument('--samtools',type=Path);a=p.parse_args()
    assert 0<a.limit<=20;review(a.root.resolve(),a.limit,a.samtools)
if __name__=='__main__':main()
