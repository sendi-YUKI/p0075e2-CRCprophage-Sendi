"""Research community reference: exact-identical sequence groups, own host-core markers.
Close nonidentical targets remain separate and shared reads remain unresolved;
this is not a calibrated short-read species/strain classifier.
"""
import json,hashlib
from pathlib import Path
from collections import defaultdict
from cohort_common import fasta,table,read_json,write_json,write_table,write_fasta,sha256,digest,require,guard_output
from cohort_reference import _ranges
def verify(root,name):
    root=Path(root);m=read_json(root/name)
    require(m.get('status')=='completed','Artifact incomplete: '+name)
    inv=m.get('outputs',{})
    require(bool(inv),'Artifact lacks inventory')
    for n,h in inv.items():
        p=(root/n).resolve()
        require(p.is_relative_to(root.resolve()) and p.is_file() and sha256(p)==h,'Artifact changed: '+n)
    return m
def prepare(catalog_dir,outdir,config_file):
    c=read_json(config_file)
    require(c.get('schema_version')==1 and c.get('status') in ('synthetic','frozen'),'Research community reference requires explicit synthetic/frozen config')
    require(c.get('grouping')=='exact_sequence_and_reverse_complement_equivalence','Unimplemented grouping strategy')
    catalog=Path(catalog_dir);cm=read_json(catalog/'catalog_manifest.json')
    require(cm.get('catalog_version','').startswith('research_'),'Research catalogue adapter required')
    # Bridge exports are hash-bound, including occurrences.
    verify(catalog,'catalog_manifest.json')
    core=Path(c['host_core_directory']);hm=verify(core,'host_core_reference_manifest.json')
    require(hm.get('status')=='completed' and hm.get('mask_evidence_status')=='assessed','Host core mask evidence missing')
    out=guard_output(outdir,[catalog,core,config_file]);reps=table(catalog/'vOTU_representatives.tsv')
    seqs=fasta(catalog/'vOTU_representatives.fna');occ=table(catalog/'source_occurrences.tsv')
    hostmap=c['source_genome_to_reference_host']
    groups=defaultdict(list);canonical={}
    rc=lambda s:s.translate(str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN'))[::-1]
    for r in reps:
        sid=r.get('representative_id') or r['candidate_id']
        require(sid in seqs,'Missing representative sequence')
        ss=seqs[sid];strand='+' if ss<=rc(ss) else '-';ss=min(ss,rc(ss))
        key=hashlib.sha256(ss.encode()).hexdigest();canonical[key]=ss
        groups[key].append((r,sid,strand))
    references=[];records={};cross=[];excluded=[]
    for key,rr in sorted(groups.items()):
        qg='qg_'+digest([cm['catalog_version'],key])[:24];rid='target_'+qg;hosts=set()
        for r,sid,strand in rr:
            sh={hostmap[o['source_host_genome']] for o in occ if o.get('viral_sequence_id')==sid and o.get('source_host_genome') in hostmap and o.get('slice_verification') in ('verified_this_run','verified')}
            hosts.update(sh);cross.append(dict(quantification_group=qg,votu_id=r['votu_id'],viral_sequence_id=sid,orientation_to_reference=strand,reference_host_candidates=';'.join(sorted(sh)),catalog_tier=r.get('catalog_tier','unknown'),member_role='representative'))
        seq=canonical[key];records[rid]=seq
        refs=dict(reference_id=rid,reference_kind='target',catalog_version=cm['catalog_version'],votu_id='',
            quantification_group=qg,biological_votu_ids=';'.join(sorted({r['votu_id'] for r,_,_ in rr})),
            reference_host_id=next(iter(hosts)) if len(hosts)==1 else '',reference_host_candidates=';'.join(sorted(hosts)),
            host_assignment_status='single_reference_host_not_sample_phasing' if len(hosts)==1 else 'multiple_or_unknown_reference_hosts',
            candidate_id='',viral_sequence_id='',source_host_genome='',candidate_type='research_sequence_group',virus_scope='unknown',
            length=len(seq),sequence_sha256=key,assessable_intervals=_ranges(seq),masked_bases=0)
        references.append(refs)
    markerrows=table(core/'host_core_markers.tsv');markers=fasta(core/'host_core.fna')
    require({r['reference_id'] for r in markerrows}==set(markers),'Host marker inventory mismatch')
    for r in markerrows:
        rid=r['reference_id'];ss=markers[rid]
        require(rid not in records,'Reference collision')
        records[rid]=ss;references.append(dict(reference_id=rid,reference_kind='host',host_id=r['host_id'],length=len(ss),
            sequence_sha256=hashlib.sha256(ss.encode()).hexdigest(),assessable_intervals=_ranges(ss),masked_bases=ss.count('N'),host_reference_kind='nonmobile_core_marker'))
    for name,path in [('excluded_host_core',core/'excluded_decoys.fna'),('declared_near_relative',Path(c['decoy_fasta']))]:
        require(path.is_file(),'Decoys missing: '+str(path))
        for ident,ss in fasta(path).items():
            rid='decoy_'+digest([name,ident,ss])[:24]
            records[rid]=ss;references.append(dict(reference_id=rid,reference_kind='decoy',length=len(ss),sequence_sha256=hashlib.sha256(ss.encode()).hexdigest(),assessable_intervals=_ranges(ss),masked_bases=0,decoy_source=name))
    require(bool(records),'Empty unified reference')
    write_fasta(out/'reference.fna',records)
    write_table(out/'quantification_group_members.tsv',['quantification_group','votu_id','viral_sequence_id','orientation_to_reference','reference_host_candidates','catalog_tier','member_role'],cross)
    write_table(out/'reference_map.tsv',['reference_id','reference_kind','quantification_group','biological_votu_ids','reference_host_id','host_assignment_status','host_id','length','sequence_sha256'],references)
    write_table(out/'host_mask.bed',['reference_id','start0','end0','reason','host_id','original_contig_id'],[])
    write_table(out/'scope_exclusions.tsv',['candidate_id','votu_id','reason'],excluded)
    m=dict(schema_version=1,status='completed',measurement_contract='research_quantification_group_v1',catalog_version=cm['catalog_version'],
        catalog_scope='research_community',catalog_manifest_sha256=sha256(catalog/'catalog_manifest.json'),host_core_manifest_sha256=sha256(core/'host_core_reference_manifest.json'),
        config_sha256=sha256(config_file),grouping=c['grouping'],grouping_limit='Exact-equivalence collapse only; nonidentical near targets remain unresolved where reads are shared; no abundance spreading',
        discovery_samples=cm.get('discovery_samples',[]),discovery_provenance_status=cm.get('discovery_provenance_status','unknown'),
        competition_evidence='unified_nonmobile_host_core_phage_declared_decoys',specificity_scope='supplied_reference_only',references=references,
        reference_fasta_sha256=sha256(out/'reference.fna'),quantification_group_members_sha256=sha256(out/'quantification_group_members.tsv'),
        sample_host_phasing=False,scientific_calibration=False)
    m['mapping_reference_version']='mapref_'+digest(m)[:24];write_json(out/'reference_manifest.json',m);return m
