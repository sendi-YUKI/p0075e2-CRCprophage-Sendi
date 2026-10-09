#!/usr/bin/env python3
"""Nonmobile core-marker producer from canonical genes and assessed coordinate masks."""
import argparse,json,os,sys,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(os.environ.get('CRC_SPARK_REPO','/srv/CRC-PHIRE/projects/crc-pipeline'))/'bin'))
from collections import defaultdict
from mainline_common import read_table,write_table,fasta,write_fasta,save,sha,new_output,run,unique
import spark_runtime
def execute(a):
    allocation=spark_runtime.allocation(require_lock=True)
    c=json.loads(Path(a.config).read_text())
    if c.get('schema_version')!=1 or c.get('status') not in ('synthetic','frozen'):raise ValueError('Explicit synthetic/frozen host-core rule required')
    for k in ('homology_identity_percent','homology_min_alignment_bp','blast_evalue','word_size'):
        if c.get(k) is None:raise ValueError('Host marker selection parameter pending: '+k)
    if not 0<float(c['homology_identity_percent'])<=100 or not 1<=int(c['homology_min_alignment_bp']) or not 4<=int(c['word_size'])<=64:raise ValueError('Invalid marker homology parameters')
    if not 1<=a.cpus<=min(6,allocation['cpus']):raise ValueError('CPU bounds')
    spark_runtime.validate_output_inputs(a.outdir,[a.config,c['phage_fasta'],*[u[k] for u in c['hosts'] for k in ('cluster_directory','alignment_directory','mask_directory')]])
    out=new_output(a.outdir);genes={};sequences={};owners={};source=[];excluded=defaultdict(set);picked=[];background={};background_genes={}
    seenhost=set()
    for unit in c['hosts']:
        host=unit['host_id'];gid=unit['representative_genome_id']
        if host in seenhost:raise ValueError('Duplicate host reference')
        seenhost.add(host)
        cluster=Path(unit['cluster_directory']);align=Path(unit['alignment_directory']);mask=Path(unit['mask_directory'])
        mm=json.loads((mask/'mask_manifest.json').read_text())
        if mm.get('status') not in ('completed','masked'):raise ValueError('Host mask is not assessed')
        gm=unique(read_table(cluster/'canonical_genes.tsv'),'gene_id');ss=fasta(cluster/'canonical.ffn')
        for gid_key,g in gm.items():
            if g['genome_id']!=gid:continue
            bid='background_'+hashlib.sha256((host+'::'+gid_key).encode()).hexdigest()[:24]
            if gid_key not in ss or hashlib.sha256(ss[gid_key].encode()).hexdigest()!=g['cds_sha256']:raise ValueError('Canonical background CDS changed')
            background[bid]=ss[gid_key];background_genes[bid]=(host,gid_key)
        intervals=read_table(mask/'mask_intervals.tsv')
        members=read_table(align/'core_members.tsv')
        bygroup=defaultdict(list)
        for r in members:
            if r['selected_core'].lower()=='true' and int(r['copy_count'])==1:
                bygroup[r['group_id']].append(r)
        for group,rr in sorted(bygroup.items()):
            selected=[r for r in rr if r['genome_id']==gid]
            if len(selected)!=1:continue
            key=host+'::'+group;item=selected[0];gene=gm[item['gene_id']]
            if gene['genome_id']!=gid or gene['gene_id'] not in ss:raise ValueError('Marker gene identity mismatch')
            seq=ss[gene['gene_id']]
            if hashlib.sha256(seq.encode()).hexdigest()!=gene['cds_sha256']:raise ValueError('Canonical CDS changed')
            rid='core_'+hashlib.sha256(key.encode()).hexdigest()[:24]
            genes[rid]=dict(reference_id=rid,host_id=host,group_id=group,genome_id=gid,gene_id=gene['gene_id'],contig_id=gene['contig_id'],start0=gene['start0'],end0=gene['end0'],strand=gene['strand'],length=len(seq),sequence_sha256=gene['cds_sha256'])
            sequences[rid]=seq;owners[rid]=(host,group)
            # Any observed mobile overlap in any member excludes the whole marker.
            for member in rr:
                g=gm[member['gene_id']]
                if any(m['genome_id']==g['genome_id'] and m['contig_id']==g['contig_id'] and max(int(g['start0']),int(m['start0']))<min(int(g['end0']),int(m['end0'])) for m in intervals):
                    excluded[rid].add('mobile_or_prophage_overlap_in_core_member')
            if set(seq)-set('ACGT'):excluded[rid].add('non_acgt_marker')
        source.extend(dict(path=str(p.resolve()),sha256=sha(p)) for p in [cluster/'canonical_genes.tsv',cluster/'canonical.ffn',align/'core_members.tsv',mask/'mask_intervals.tsv',mask/'mask_manifest.json'])
    # Search target phages, core markers and all canonical genes of declared reference genomes.
    phages=fasta(c['phage_fasta']);subjects=dict(sequences);subjects.update(background)
    for i,seq in enumerate(phages.values()):subjects['phage_'+str(i)]=seq
    write_fasta(out/'search_markers.fna',sequences);write_fasta(out/'competition_sequences.fna',subjects)
    raw=out/'homology.tsv'
    if sequences and subjects:
        run(['blastn','-task','blastn','-query',out/'search_markers.fna','-subject',out/'competition_sequences.fna','-dust','no',
            '-word_size',str(c['word_size']),'-evalue',str(c['blast_evalue']),'-num_threads',str(a.cpus),
            '-max_target_seqs',str(max(5,len(subjects))),'-outfmt','6 qseqid sseqid pident length qstart qend sstart send','-out',raw],out,'marker_homology',timeout=a.timeout)
    else:raw.write_text('')
    for line in raw.read_text().splitlines():
        q,t,identity,length,qa,qb,ta,tb=line.split('\t')
        if q not in sequences or t not in subjects:raise ValueError('Marker BLAST foreign key')
        if float(identity)<c['homology_identity_percent'] or int(length)<c['homology_min_alignment_bp']:continue
        if t.startswith('phage_'):excluded[q].add('phage_homology')
        elif t in background_genes:
            if background_genes[t]!=(genes[q]['host_id'],genes[q]['gene_id']):excluded[q].add('nonself_canonical_gene_homology_or_repeat')
        elif q!=t:excluded[q].add('cross_marker_homology_or_repeat')
        elif (int(qa),int(qb))!=(int(ta),int(tb)):excluded[q].add('within_marker_repeat')
    retained={rid:s for rid,s in sequences.items() if not excluded[rid]}
    write_fasta(out/'host_core.fna',retained)
    # Phage-homologous copies are excluded instead of duplicated as decoys.
    # Other failed core markers remain explicit competitors.
    write_fasta(out/'excluded_decoys.fna',{rid:s for rid,s in sequences.items() if excluded[rid] and 'phage_homology' not in excluded[rid]})
    fields=['reference_id','host_id','group_id','genome_id','gene_id','contig_id','start0','end0','strand','length','sequence_sha256']
    write_table(out/'host_core_markers.tsv',[genes[rid] for rid in retained],fields)
    write_table(out/'marker_exclusions.tsv',[dict(genes[rid],reason=';'.join(sorted(excluded[rid]))) for rid in genes if excluded[rid]],fields+['reason'])
    save(out/'host_core_reference_manifest.json',dict(status='completed',mask_evidence_status='assessed',config_sha256=sha(a.config),
        phage_fasta_sha256=sha(c['phage_fasta']),sources=source,retained_markers=len(retained),host_marker_counts={h:sum(genes[g]['host_id']==h for g in retained) for h in seenhost},
        interpretation='Evidence-filtered core markers under explicit rule, not guaranteed mobile-free or species-exclusive outside supplied panel',
        scientific_calibration=False,outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--outdir',required=True);p.add_argument('--cpus',type=int,default=1);p.add_argument('--timeout',type=int,default=900)
    a=p.parse_args()
    if not 1<=a.timeout<=3600:raise ValueError('Timeout cap')
    execute(a)
