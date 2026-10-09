#!/usr/bin/env python3
"""Combine independently validated A2/A3/A4 outputs; never promote failed stages."""
import argparse
import csv
import hashlib
import html
import json
import shutil
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''):h.update(b)
    return h.hexdigest()


def table(path):
    with Path(path).open() as f:
        r=csv.DictReader(f,delimiter='\t'); rows=list(r)
        if not r.fieldnames or any(None in x or None in x.values() for x in rows):
            raise ValueError('Malformed table '+str(path))
        return r.fieldnames,rows


def write(path,fields,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def merge(dirs,name,out,required=True):
    fields=None; rows=[]
    for directory in dirs:
        f=Path(directory)/name
        if not f.exists():
            if required:raise ValueError('Missing required output '+str(f))
            continue
        flds,part=table(f)
        if fields is not None and flds != fields:raise ValueError('Incompatible schemas '+name)
        fields=flds;rows.extend(part)
    if fields:
        rows.sort(key=lambda r:tuple(r[k] for k in fields))
        write(out/name,fields,rows)
    return rows


def unique(rows,key):
    result={r[key]:r for r in rows}
    if len(result)!=len(rows):raise ValueError('Duplicate '+key)
    return result


def validate_caller_status(row, hits):
    n=sum(h['genome_id']==row['genome_id'] and h['caller']==row['caller'] for h in hits)
    if row['status']=='not_assessed':
        if row['caller'] not in {'PhiSpy','PhageBoost'} or row['exit_code'] or row['candidate_count'] or not row.get('reason') or n:
            raise ValueError('Inconsistent unassessed caller evidence')
        return
    if row['status'] not in {'completed_hits','completed_zero'} or str(row['exit_code'])!='0':
        raise ValueError('Caller failed or unresolved execution')
    if int(row['candidate_count'])!=n or (row['status']=='completed_zero')!=(n==0):
        raise ValueError('Caller count/hit mismatch')



def neighborhood_browser(out, genes, candidates, limit=100):
    """Display-only bounded preview; full coordinates remain in source TSVs."""
    import bisect
    groups={}
    for g in genes:
        key=(g.get('genome_id'),g.get('contig_id'))
        try:start,end=int(g['start0']),int(g['end0'])
        except (KeyError,ValueError):continue
        groups.setdefault(key,[]).append((start,end,g))
    for key in groups:groups[key].sort(key=lambda x:(x[0],x[1],x[2]['gene_id']))
    h=lambda x:html.escape(str(x),quote=True)
    panels=[]
    for c in sorted(candidates.values(),key=lambda r:r['candidate_id'])[:limit]:
        try:start,end=int(c['start0']),int(c['end0'])
        except (KeyError,ValueError):continue
        left=max(0,start-5000);right=end+5000;width=max(1,right-left)
        items=groups.get((c.get('genome_id'),c.get('contig_id')),[])
        stop=bisect.bisect_left([r[0] for r in items],right)
        selected=[r for r in items[:stop] if r[1]>left];shown=selected[:200]
        scale=lambda x:20+960*(x-left)/width
        svg='<svg viewBox="0 0 1000 110" role="img"><rect x="'+str(scale(start))+'" y="10" width="'+str(960*(end-start)/width)+'" height="90" fill="#def0f5"/>'
        for i,(a,b,g) in enumerate(shown):
            x,y,z=scale(max(left,a)),35+(i%2)*35,scale(min(right,b))
            direction=g.get('strand','unknown');tip=5
            points=[(x,y-8),(z-tip,y-8),(z,y),(z-tip,y+8),(x,y+8)] if direction in ('+','1') else [(x,y),(x+tip,y-8),(z,y-8),(z,y+8),(x+tip,y+8)]
            svg+='<polygon fill="#285a73" points="'+' '.join(str(u)+','+str(v) for u,v in points)+'"><title>'+h(g['gene_id'])+' '+h(a)+'..'+h(b)+' strand '+h(direction)+'</title></polygon>'
        svg+='</svg>'
        panels.append('<details><summary>'+h(c['candidate_id'])+' / '+h(c.get('genome_id'))+' / '+h(c.get('contig_id'))+'</summary><p>0-based half-open; shaded candidate '+h(start)+'..'+h(end)+'. Window '+h(left)+'..'+h(right)+'. Showing '+str(len(shown))+'/'+str(len(selected))+' genes. Hover a gene for ID/strand.</p>'+svg+'</details>')
    page='<!doctype html><meta charset="utf-8"><title>Candidate coordinate browser</title><style>body{font:16px system-ui;max-width:1150px;margin:2em auto}svg{width:100%}details{margin:1em}</style><h1>Candidate coordinate browser</h1><p>Display preview of first '+str(limit)+' sorted candidates; at most200 genes/window. No new biological inference or synteny grouping. Full results: <a href="prophage_master.tsv">candidates</a>, <a href="gene_table.tsv">genes</a>, <a href="function_hits.tsv">own-gene functions</a>, <a href="system_prophage_membership.tsv">boundary systems</a>. Preview limits do not filter analysis.</p>'+''.join(panels)
    (out/'candidate_neighborhoods.html').write_text(page)
    return dict(displayed_candidates=len(panels),candidate_preview_limit=limit,genes_per_window_limit=200,coordinate_system='zero_based_half_open',scientific_filter=False)


def collect_research_outputs(a, out, catalog, genemap, candidates):
    extensions = {}
    if bool(a.locus)!=bool(a.research_summary):raise ValueError('Incomplete research extension')
    if a.research_summary:
        receipt=json.loads((Path(a.research_summary)/'research_summary.json').read_text())
        if receipt['catalog_version']!=catalog['catalog_version']:raise ValueError('Research catalog version differs')
        for row in table(Path(a.research_summary)/'research_member_gene_links.tsv')[1]:
            if row['gene_id'] not in genemap or row['candidate_id'] not in candidates or row['source_host_genome']!=genemap[row['gene_id']]['genome_id']:raise ValueError('Research member/gene/source-host mismatch')
        for directory, manifest_name in [(a.locus,'locus_manifest.json'),(a.research_summary,'research_summary.json')]:
            md=json.loads((Path(directory)/manifest_name).read_text())
            if md['catalog_version']!=catalog['catalog_version']:raise ValueError('Mixed research output versions')
            for name, expected in md['outputs'].items():
                if Path(name).name!=name or sha(Path(directory)/name)!=expected:raise ValueError('Changed research output '+name)
            for f in sorted(Path(directory).iterdir()):
                if f.is_file() and f.suffix in ('.tsv','.json','.fna'):
                    if (out/f.name).exists():raise ValueError('Research output name collision '+f.name)
                    shutil.copyfile(f,out/f.name)
        extensions['catalog_research']=receipt
    return extensions.get('catalog_research', {'status':'not_enabled'})


def collect_host_source(directory,out,catalog_version,research_summary):
    if not directory:return {'status':'not_enabled'}
    folder=Path(directory);m=json.loads((folder/'host_source_manifest.json').read_text())
    if m['catalog_version']!=catalog_version or not research_summary or m['research_manifest_sha256']!=sha(Path(research_summary)/'research_summary.json'):
        raise ValueError('Host-source research/catalog provenance mismatch')
    target=out/'host_source';target.mkdir(exist_ok=False)
    for name,expected in m['outputs'].items():
        if Path(name).name!=name or sha(folder/name)!=expected:raise ValueError('Changed host-source output '+name)
        shutil.copyfile(folder/name,target/name)
    if (folder/'sensitivity_inventory.json').is_file():
        inventory=json.loads((folder/'sensitivity_inventory.json').read_text())
        for name,expected in inventory['outputs'].items():
            source=(folder/name).resolve()
            if not source.is_relative_to(folder.resolve()) or sha(source)!=expected:raise ValueError('Changed nested sensitivity output')
            destination=target/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,destination)
    shutil.copyfile(folder/'host_source_manifest.json',target/'host_source_manifest.json')
    return m


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['a1','catalog','qc','function-summary','caller-comparison','outdir']:p.add_argument('--'+name,required=True)
    for name in ['genes','taxonomy','functions']:p.add_argument('--'+name,nargs='+',required=True)
    p.add_argument('--kofam-summary');p.add_argument('--host-dirs',nargs='+')
    p.add_argument('--host-source')
    p.add_argument('--locus');p.add_argument('--research-summary')
    a=p.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    # Directory contracts below are strict: a missing module output fails.
    for folder in [a.catalog,a.qc,a.function_summary,a.caller_comparison]:
        for f in sorted(Path(folder).glob('*')):
            if f.is_file() and f.suffix in ['.tsv','.json','.fna']:
                if (out/f.name).exists():raise ValueError('Unexpected output name collision '+f.name)
                shutil.copyfile(f,out/f.name)
    for name in ['caller_hits.tsv','checkv_quality.tsv','coordinate_mapping.tsv','contig_map.tsv','genome_stats.tsv','candidates.fna','data_manifests.json']:
        if name=='candidates.fna' and (out/name).exists():continue  # Preserve hash-bound research catalogue FASTA.
        shutil.copyfile(Path(a.a1)/name,out/name)
    genes=merge(a.genes,'gene_table.tsv',out)
    # Crosswalk filename is part of A2 public interface; preserve both input IDs.
    for name in ['crosswalk.tsv']:
        merge(a.genes,name,out,required=False)
    for name in ['genes.faa','genes.ffn']:
        with (out/name).open('w') as f:
            for d in a.genes:f.write((Path(d)/name).read_text())
    taxonomy=merge(a.taxonomy,'taxonomy.tsv',out)
    merge(a.taxonomy,'reference_ani.tsv',out)
    names=['function_hits.tsv','systems.tsv','system_hits.tsv','member_functions.tsv','cargo_candidates.tsv','system_prophage_membership.tsv','candidate_function_status.tsv','gene_prophage_membership.tsv']
    for name in names:merge(a.functions,name,out)
    receipts=[]
    for directory in a.functions:
        r=json.loads((Path(directory)/'function_status.json').read_text())
        if r.get('status')!='completed':raise ValueError('Functions not completed: '+str(directory))
        receipts.append(r)
    _,master=table(out/'prophage_master.tsv'); candidates=unique(master,'candidate_id')
    genemap=unique(genes,'gene_id'); _,qc=table(out/'genome_qc.tsv');qcmap=unique(qc,'genome_id')
    _,base_status=table(Path(a.a1)/'sample_status.tsv');samples=unique(base_status,'genome_id')
    if set(qcmap)!=set(samples) or {g['genome_id'] for g in genes}!=set(samples):raise ValueError('Sample/QC/gene foreign keys differ')
    for name,key in [('gene_prophage_membership.tsv','gene_id'),('function_hits.tsv','gene_id'),('system_hits.tsv','gene_id')]:
        for row in table(out/name)[1]:
            if row[key] not in genemap:raise ValueError('Unknown gene FK '+name)
    for name in ['candidate_function_status.tsv','member_functions.tsv','gene_prophage_membership.tsv','cargo_candidates.tsv','system_prophage_membership.tsv']:
        for row in table(out/name)[1]:
            if row['candidate_id'] not in candidates:raise ValueError('Unknown candidate FK '+name)
    votus=table(out/'vOTU_members.tsv')[1]
    if set(unique(votus,'candidate_id'))!=set(candidates):raise ValueError('vOTU rows must preserve all candidates')
    status=[];taxmap=unique(taxonomy,'genome_id')
    if set(taxmap)!=set(samples):raise ValueError('Taxonomy/sample foreign keys differ')
    if set(unique(receipts,'genome_id'))!=set(samples):raise ValueError('Function receipts/sample foreign keys differ')
    function_candidates=unique(table(out/'candidate_function_status.tsv')[1],'candidate_id')
    if set(function_candidates)!=set(candidates):raise ValueError('Function candidate assessability rows missing or extra')
    caller_status=table(out/'caller_sample_status.tsv')[1]
    expected={(gid,caller) for gid in samples for caller in ['PhiSpy','PhageBoost']}
    if len(caller_status)!=len(expected) or {(r['genome_id'],r['caller']) for r in caller_status}!=expected:
        raise ValueError('Incomplete auxiliary caller sample status')
    pilot_hits=table(out/'caller_hits_pilot.tsv')[1]
    for r in caller_status:
        validate_caller_status(r,pilot_hits)
    for gid in sorted(samples):
        status.append(dict(genome_id=gid,A0_A1_status='completed',checkm2_status=qcmap[gid]['checkm2_status'],genome_qc_status=qcmap[gid]['genome_qc_status'],eligibility=qcmap[gid]['eligibility'],taxonomy_status=taxmap[gid].get('status',taxmap[gid].get('taxonomy_status','unknown')),genes_status='completed',functions_status='completed',n_candidates=sum(r['genome_id']==gid for r in master),n_genes=sum(r['genome_id']==gid for r in genes)))
    write(out/'sample_status_core.tsv',list(status[0]),status)
    counts=[]
    for name in ['genome_qc.tsv','gene_table.tsv','prophage_master.tsv','vOTU_representatives.tsv','vOTU_members.tsv','function_hits.tsv','systems.tsv','system_hits.tsv','cargo_candidates.tsv','caller_hits_pilot.tsv','caller_comparisons.tsv']:
        f=out/name
        if f.exists():counts.append(dict(stage=name.removesuffix('.tsv'),rows=len(table(f)[1])))
    write(out/'stage_counts_core.tsv',['stage','rows'],counts)
    catalog=json.loads((out/'catalog_manifest.json').read_text())
    extensions=dict(kofam={'status':'not_assessed','reason':'not_enabled'},host_phylogeny={'status':'not_assessed','reason':'not_enabled'})
    if a.kofam_summary:
        ko=Path(a.kofam_summary)
        receipt=json.loads((ko/'kofam_summary.json').read_text())
        if receipt['status']!='completed':raise ValueError('Incomplete KOfam summary')
        for row in table(ko/'gene_ko_status.tsv')[1]:
            if row['gene_id'] not in genemap or row['genome_id']!=genemap[row['gene_id']]['genome_id']:raise ValueError('KOfam gene/source-host FK differs')
        if set(unique(table(ko/'gene_ko_status.tsv')[1],'gene_id'))!=set(genemap):raise ValueError('KOfam did not assess every canonical gene')
        for f in ko.iterdir():
            if f.is_file():
                if (out/f.name).exists():raise ValueError('KOfam output name collision')
                shutil.copyfile(f,out/f.name)
        extensions['kofam']=receipt
    if a.host_dirs:
        expected={'host_alignment','preliminary_tree','host_mask','final_tree'}
        directories={Path(d).name:Path(d) for d in a.host_dirs}
        if set(directories)!=expected:raise ValueError('Incomplete host extension outputs')
        final=json.loads((directories['final_tree']/'tree_status.json').read_text())
        for row in table(directories['final_tree']/'tree_tips.tsv')[1]:
            if row['genome_id'] not in samples or row['source_host_genome']!=row['genome_id']:raise ValueError('Host tree/source-host FK differs')
        outputs={}
        for name,d in directories.items():
            for f in d.iterdir():
                if f.is_file() and f.suffix in ('.tsv','.json','.fna','.nex','.treefile'):
                    outputs['host_phylogeny/'+name+'/'+f.name]=sha(f)
        extensions['host_phylogeny']=dict(status=final['status'],outputs=outputs,kinship='method_pending')
    extensions['catalog_research']={'status':'not_enabled'}
    extensions['catalog_research']=collect_research_outputs(a,out,catalog,genemap,candidates)
    extensions['host_source']=collect_host_source(a.host_source,out,catalog['catalog_version'],a.research_summary)
    extensions['neighborhood_browser']=neighborhood_browser(out,genes,candidates)
    (out/'mainline_extensions.json').write_text(json.dumps(extensions,indent=2,sort_keys=True)+'\n')
    manifest=dict(current_target='A2_A3_A4' ,status='completed',catalog_version=catalog['catalog_version'],
                  baseline_master_sha256=sha(Path(a.a1)/'prophage_master.tsv'),
                  method_status='engineering_implemented_and_testable_scientific_rules_provisional',
                  retained_boundaries='geNomad raw v1 unchanged; PhiSpy/PhageBoost comparison only',
                  mainline_extensions=extensions,functions=receipts,outputs={f.name:sha(f) for f in sorted(out.iterdir()) if f.is_file()})
    if extensions['catalog_research'].get('status')=='completed_with_failures':manifest['status']='completed_with_failures'
    (out/'core_data_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    def h(x):return html.escape(str(x))
    rows=''.join('<tr><td>'+h(r['stage'])+'</td><td>'+h(r['rows'])+'</td></tr>' for r in counts)
    links=''.join('<li><a href="'+h(f.name)+'">'+h(f.name)+'</a></li>' for f in sorted(out.glob('*.tsv')))
    content='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CRC phage core report</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;line-height:1.6;color:#142c38}table{border-collapse:collapse}td,th{padding:8px 20px;border:1px solid #ccc}code{background:#eef3f6}a{color:#076080}</style><h1>CRC phage pipeline: A2 + A3 + A4</h1><p>Catalog: <code>'+h(catalog['catalog_version'])+'</code></p><p>Primary boundaries remain geNomad predictions. Host QC, viral quality, integration evidence and boundary confidence are separate. Functional absence means absence under these models and thresholds; AMG labels are candidates requiring review. Auxiliary callers are compared, not unioned into the catalog.</p><table><tr><th>Output</th><th>Rows</th></tr>'+rows+'</table><h2>Tables</h2><ul>'+links+'</ul><p>Runtime provenance, parameters, software and database manifests: <a href="pipeline_info/core_run_manifest.json">pipeline_info/core_run_manifest.json</a>. Nextflow trace/report/timeline are retained in pipeline_info.</p></html>'
    content=content.replace('</html>','<h2>Optional mainline evidence</h2><p>KOfam: '+h(extensions['kofam']['status'])+'. Host phylogeny: '+h(extensions['host_phylogeny']['status'])+'. KO means homology evidence, not AMG/activity; source host is not a predicted host. A distance matrix is not kinship.</p><p><a href="mainline_extensions.json">Extension identities and status</a></p></html>')
    if a.research_summary:content=content.replace('</html>','<h2>Research catalog and locus evidence</h2><p>Status: '+h(extensions['catalog_research']['status'])+'. Unknown is not absence. Locus empty-site support is not whole-genome vOTU absence; fragment functions apply only to observed member sequence.</p></html>')
    if a.host_source:content=content.replace('</html>','<h2>CRC-source model</h2><p>Conditional association; unknown is excluded, not zero. Scientific calibration and real workflow validation pending.</p><p><a href="host_source/model_results.tsv">Effects and non-estimable reasons</a> | <a href="host_source/F1_results.tsv">Global F1 status</a></p></html>')
    content=content.replace('</html>','<p><a href="candidate_neighborhoods.html">Candidate coordinate / gene neighborhood browser</a></p></html>')
    (out/'report_core.html').write_text(content)

if __name__=='__main__':main()
