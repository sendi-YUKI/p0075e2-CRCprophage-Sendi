"""Metadata-only selection. No sequences downloaded, no sequence QC claimed."""
import csv,json,re,hashlib,collections,urllib.request,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCOPE=json.loads((ROOT/'selection_scope.json').read_text(encoding='utf-8'))
def table(path,rows,fields=None):
    rows=list(rows); fields=fields or list(rows[0])
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
old=read(ROOT.parent/'reference258_manifest.tsv')
reports={}; queries=collections.defaultdict(set)
for unit in SCOPE['taxa']:
    for p in sorted((ROOT/'raw'/unit['key']).glob('page_*.json')):
        if '.receipt.' in p.name:continue
        for r in json.loads(p.read_text(encoding='utf-8'))['reports']:
            reports[r['accession']]=r;queries[r['accession']].add(unit['key'])
# Obtain authoritative current reports for historical inputs not returned by taxon query.
missing=[r['accession'] for r in old if r['accession'] not in reports]
cache=ROOT/'raw'/'old258_lookup';cache.mkdir(exist_ok=True)
for i in range(0,len(missing),25):
    subset=missing[i:i+25];p=cache/f'batch_{i//25+1:03}.json'
    if not p.exists():
        url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+','.join(subset)+'/dataset_report?page_size=1000'
        for retry in range(4):
            try:
                data=urllib.request.urlopen(url,timeout=90).read();p.write_bytes(data)
                p.with_suffix('.receipt.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(data).hexdigest()},indent=2));break
            except Exception:
                if retry==3:raise
                time.sleep(2*(retry+1))
    for r in json.loads(p.read_text()) .get('reports',[]):
        reports[r['accession']]=r;queries[r['accession']].add('old258_lookup')
    time.sleep(.4)
units={x['key']:x for x in SCOPE['taxa']}
oldmap={x['accession']:x for x in old}
missingwords={'','missing','not provided','not collected','not applicable','unknown','na','n/a','not available','not determined'}
def clean(x):
    v=str(x or '').replace('\t',' ').replace('\r',' ').replace('\n',' ').strip()
    return '' if v.lower() in missingwords else v
def num(v,default=0):
    try:return float(v)
    except:return default
def normalize(r):
    a=r.get('assembly_info',{});b=a.get('biosample',{});s=r.get('assembly_stats',{});o=r.get('organism',{});c=r.get('checkm_info',{})
    attr={x['name'].lower():clean(x.get('value')) for x in b.get('attributes',[])}
    def val(*ks):
        return next((clean(b.get(k) or attr.get(k)) for k in ks if clean(b.get(k) or attr.get(k))),'')
    accession=r['accession'];name=o.get('organism_name','');scientific=r.get('average_nucleotide_identity',{}).get('submitted_species') or name
    keys=[k for k in queries[accession] if k in units]
    exact=[k for k in keys if scientific.lower()==units[k]['query'].lower()]
    if exact:key=exact[0]
    elif keys:key=sorted(keys,key=lambda k:(units[k]['scope_rank']=='genus',len(units[k]['query'])),reverse=True)[0]
    else:
        oldtax=oldmap.get(accession,{}).get('adopted_taxon','')
        names=scientific+' '+oldtax
        key=next((k for k,v in units.items() if v['query'].lower() in names.lower()),'')
    unit=units.get(key,{})
    host=val('host');source=val('isolation_source','env_medium');site=val('body_site','host_body_site','env_local_scale')
    combined=' '.join([host,source,site]).lower()
    human=bool(re.search(r'\bhuman\b|homo[ _]sapiens|h\. sapiens|\bpatient\b',combined)) or host.lower()=='homo'
    gut=bool(re.search(r'fecal|faecal|feces|faeces|stool|colon|colorectal|intestin|\bgut\b|rectal',combined))
    nonhuman=bool(re.search(r'mus musculus|mouse|mice|chicken|gallus|porcine|swine|\bpig\b|bovine|cattle|\bcow\b|canine|\bdog\b|rat\b|poultry|sewage|wastewater|soil|river|water sample|food',combined))
    # Explicit non-human host labels cannot enter as generic unknown human hosts.
    if host and not human and not re.search(r'unknown|not |mammal|animal|metagenome|strain|free-living',host.lower()):nonhuman=True
    if human:sourcegroup='human_gut' if gut else 'human_other'
    elif nonhuman:sourcegroup='nonhuman_or_environmental'
    elif host and not re.search('mammal|animal|metagenome',host.lower()):sourcegroup='unresolved_host'
    else:sourcegroup='gut_host_unknown' if gut else 'unknown'
    if re.search(r'oral|mouth|dental|saliva|gingiva|plaque',combined):body='oral'
    elif gut:body='gut'
    elif re.search('blood|bacteremia',combined):body='blood'
    elif site:body=site
    else:body='unknown'
    explicit=val('host_disease','disease','host_health_state','health_state','disease_status')
    disease='unknown'
    if re.search(r'colorectal cancer|colon cancer|colorectal carcinoma|\bcrc\b',explicit.lower()):disease='CRC_reported_unverified'
    elif re.search(r'healthy|normal|no disease',explicit.lower()):disease='healthy_reported_unverified'
    elif explicit:disease='other_or_unresolved_reported'
    strain=clean(o.get('infraspecific_names',{}).get('strain') or val('strain') or o.get('infraspecific_names',{}).get('isolate'))
    typeflag=bool(r.get('type_material')) or 'type strain' in val('type-material').lower()
    return dict(stable_candidate_id='REF_'+accession.replace('.','_'),accession=accession,current_accession=r.get('current_accession',''),paired_accession=r.get('paired_accession') or a.get('paired_assembly',{}).get('accession',''),query_key=key,query_memberships=';'.join(sorted(queries[accession])),tax_id=o.get('tax_id',''),taxon=scientific,organism_name=name,strain_original=strain,strain_normalized=re.sub(r'\s+',' ',strain),culture_collection_id=val('culture_collection'),biosample=b.get('accession',''),bioproject=a.get('bioproject_accession',''),source_database=r.get('source_database',''),isolation_source=source,body_site=body,host=host,source_group=sourcegroup,country=val('geo_loc_name','geographic location','country'),collection_date=val('collection_date'),disease_label=disease,disease_label_original=explicit,disease_reliability='public_metadata_only_not_adjudicated',eligible_for_disease_comparison='false',assembly_level=a.get('assembly_level',''),assembly_status=a.get('assembly_status',''),release_date=a.get('release_date',''),genome_size=s.get('total_sequence_length',''),contig_count=s.get('number_of_contigs',''),contig_n50=s.get('contig_n50',''),checkm_completeness=c.get('completeness',''),checkm_contamination=c.get('contamination',''),checkm_method=c.get('checkm_version',''),checkm_missing_policy='absent_API_values_not_certified_zero',taxonomy_check_status=r.get('average_nucleotide_identity',{}).get('taxonomy_check_status',''),taxonomy_match_status=r.get('average_nucleotide_identity',{}).get('match_status',''),taxonomy_snapshot='NCBI_Datasets_2026-10-08',type_material=str(typeflag).lower(),public_lineage=val('mlst','sequence_type','serotype','phylogroup'),lineage_evidence='BioSample_attribute_only' if val('mlst','sequence_type','serotype','phylogroup') else 'not_available',candidate_role=unit.get('role','unresolved'),disease_direction=unit.get('disease_direction','unresolved'),source_ids=';'.join(unit.get('source_ids',[])),old258_exact=str(accession in oldmap).lower(),selection_status='',selection_reason='',preferred_accession='',identity_group='',sequence_QC='not_performed_here',pending_checks='sequence_QC;ANI_taxonomy;isolate_not_MAG;strain_aliases;phylogeny;FASTA_SHA256',metadata_url='https://www.ncbi.nlm.nih.gov/datasets/genome/'+accession+'/',biosample_url='https://www.ncbi.nlm.nih.gov/biosample/'+b.get('accession',''),fasta_sha256='not_downloaded_here',official_size='not_retrieved',official_md5='not_retrieved')
rows=[normalize(r) for r in reports.values()]
# Union only explicit paired assemblies and same BioSample; never merge by species ANI.
parent={r['accession']:r['accession'] for r in rows}
def find(x):
    while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
    return x
def union(a,b):
    if b in parent:parent[find(a)]=find(b)
bios={}
for r in rows:
    union(r['accession'],r['paired_accession'])
    if r['biosample']:
        if r['biosample'] in bios:union(r['accession'],bios[r['biosample']])
        else:bios[r['biosample']]=r['accession']
def quality(r):
    fail=int((r['checkm_completeness']!='' and num(r['checkm_completeness'])<90) or (r['checkm_contamination']!='' and num(r['checkm_contamination'])>5))
    return (r['assembly_status'] not in ('','current'),fail,{'Complete Genome':0,'Chromosome':1,'Scaffold':2,'Contig':3}.get(r['assembly_level'],4),num(r['contig_count'],1e9),-num(r['contig_n50']),not r['accession'].startswith('GCF'),r['accession'])
groups=collections.defaultdict(list)
for r in rows:groups[find(r['accession'])].append(r)
representatives=[]
for members in groups.values():
    best=min(members,key=quality);gid='IDENTITY_'+hashlib.sha256('|'.join(sorted(x['accession'] for x in members)).encode()).hexdigest()[:16]
    for r in members:
        r['identity_group']=gid;r['preferred_accession']=best['accession']
        if r is not best:r['selection_status']='duplicate_assembly';r['selection_reason']='paired_GCA_GCF_or_same_BioSample;see_preferred_accession'
    if any(r['assembly_status']=='suppressed' for r in members):best['pending_checks']+=';paired_or_same_BioSample_record_suppressed_check_reason'
    representatives.append(best)
# Named strain collisions across BioSamples are flagged, not automatically collapsed.
aliases=collections.defaultdict(list)
for r in representatives:
    v=re.sub('[^a-z0-9]','',r['strain_original'].lower())
    if len(v)>=3:aliases[(r['taxon'],v)].append(r)
for members in aliases.values():
    if len(members)>1:
        for r in members:r['pending_checks']+=';same_strain_name_multiple_BioSamples'
for r in representatives:
    reason=''
    quality_warning=(r['checkm_completeness']!='' and num(r['checkm_completeness'])<90) or (r['checkm_contamination']!='' and num(r['checkm_contamination'])>5)
    # A complete/type structural anchor with a legacy marker-QC warning is retained
    # conditionally, never certified as passing; avoids deleting an entire lineage.
    anchor_exception=quality_warning and (r['assembly_level']=='Complete Genome' or r['type_material']=='true') and num(r['contig_count'],1e9)<=200 and num(r['checkm_completeness'],100)>=70 and num(r['checkm_contamination'])<=15
    if anchor_exception:r['pending_checks']+=';MANDATORY_legacy_CheckM_warning_structural_anchor_recheck'
    if not r['query_key']:reason='taxonomy_scope_unresolved'
    elif r['assembly_status'] not in ('','current'):reason='noncurrent_assembly'
    elif r['checkm_completeness']!='' and num(r['checkm_completeness'])<90 and not anchor_exception:reason='reported_completeness_below_90_recheck'
    elif r['checkm_contamination']!='' and num(r['checkm_contamination'])>5 and not anchor_exception:reason='reported_contamination_above_5_recheck'
    elif num(r['contig_count'])>500:reason='over_500_contigs'
    elif r['source_group']=='nonhuman_or_environmental' and r['type_material']!='true' and r['query_key']!='butyricicoccus_pullicaecorum':reason='nonhuman_context_outside_primary_reference_scope'
    if reason:
        r['selection_status']='review' if reason=='taxonomy_scope_unresolved' or reason.endswith('_recheck') else 'excluded';r['selection_reason']=reason
    elif r['source_group']=='nonhuman_or_environmental':
        r['selection_status']='context_only';r['selection_reason']='nonhuman_type_or_experimental_context_not_patient_host'
    elif num(r['contig_count'])>200:
        r['selection_status']='reserve';r['selection_reason']='201_to_500_contigs_pending_structural_quality'
    elif r['taxonomy_check_status'] not in ('','OK','Inconclusive'):
        r['selection_status']='review';r['selection_reason']='NCBI_taxonomy_warning'
    else:r['selection_status']='eligible'
priority={'human_gut':0,'human_other':1,'gut_host_unknown':2,'unknown':3,'unresolved_host':4}
for key,unit in units.items():
    eligible=[r for r in representatives if r['query_key']==key and r['selection_status']=='eligible']
    buckets=collections.defaultdict(list)
    for r in eligible:buckets[(r['bioproject'] or r['accession'],r['country'].split(':')[0],r['body_site'],r['taxon'])].append(r)
    for b in buckets.values():b.sort(key=lambda r:(priority.get(r['source_group'],5),quality(r)))
    chosen=[]
    # Retain type representatives first, then round-robin across study/site/taxon strata.
    for r in sorted(eligible,key=quality):
        if r['type_material']=='true' and len(chosen)<unit['candidate_budget_cap']:chosen.append(r)
    while len(chosen)<unit['candidate_budget_cap']:
        progressed=False
        for k in sorted(buckets):
            b=buckets[k]
            while b and b[0] in chosen:b.pop(0)
            if b and len(chosen)<unit['candidate_budget_cap']:chosen.append(b.pop(0));progressed=True
        if not progressed:break
    ids={r['accession'] for r in chosen}
    for r in eligible:
        r['selection_status']='preselected' if r['accession'] in ids else 'reserve'
        r['selection_reason']='metadata_eligible_stratified_budget_selection_pending_sequence_QC' if r['accession'] in ids else 'eligible_beyond_taxon_budget'
rows.sort(key=lambda r:(r['query_key'],r['accession']))
table(ROOT/'all_assembly_candidates.tsv',rows)
for status in ['preselected','reserve','excluded','review','context_only','duplicate_assembly']:
    table(ROOT/(status+'.tsv'),[r for r in rows if r['selection_status']==status],list(rows[0]))
selected=[r for r in rows if r['selection_status']=='preselected']
table(ROOT/'download_manifest.tsv',[dict(accession=r['accession'],stable_candidate_id=r['stable_candidate_id'],taxon=r['taxon'],strain=r['strain_original'],role=r['candidate_role'],expected_stage='download_then_QC',metadata_url=r['metadata_url'],download_command='datasets download genome accession '+r['accession']+' --include genome --filename '+r['accession']+'.zip',official_md5='retrieve_at_download',local_sha256='calculate_after_download') for r in selected])
(ROOT/'download_accessions.txt').write_text('\n'.join(r['accession'] for r in selected)+'\n')
byacc={r['accession']:r for r in rows};reconciliation=[]
for x in old:
    r=byacc.get(x['accession']);best=byacc.get(r['preferred_accession']) if r else None
    reconciliation.append(dict(old_accession=x['accession'],old_taxon=x['adopted_taxon'],old_strain=x['strain_name'],preferred_accession=best['accession'] if best else '',new_taxon=best['taxon'] if best else '',new_status=best['selection_status'] if best else 'review_missing_public_report',retained_identity=str(bool(best and best['selection_status']=='preselected')).lower(),exact_assembly_retained=str(bool(best and best['selection_status']=='preselected' and best['accession']==x['accession'])).lower(),reason=best['selection_reason'] if best else 'public_lookup_not_resolved',local_sha256=x['sha256'],sequence_equivalence='not_assumed_between_GCA_GCF_or_versions'))
table(ROOT/'old258_reconciliation.tsv',reconciliation)
counts=collections.Counter(r['selection_status'] for r in rows)
retained=sum(x['retained_identity']=='true' for x in reconciliation)
oldgroups={byacc[x['accession']]['identity_group'] for x in old if x['accession'] in byacc}
summary=dict(query_records=sum(x['records'] for x in json.loads((ROOT/'query_receipts.json').read_text())),unique_assembly_records=len(rows),metadata_identity_groups=len(groups),status_counts=dict(counts),preselected_genomes=len(selected),preselected_taxon_labels=len({r['taxon'] for r in selected}),preselected_query_units=len({r['query_key'] for r in selected}),old258_retained_identities=retained,old258_exact_assemblies_retained=sum(x['exact_assembly_retained']=='true' for x in reconciliation),new_identity_candidates=sum(r['identity_group'] not in oldgroups for r in selected),named_strain_count=sum(bool(r['strain_original']) for r in selected),final_sequence_QC='pending_Jinlong',sequence_level_independent_strains='not_yet_established')
(ROOT/'selection_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
taxrows=[]
for key,u in units.items():
    rr=[r for r in rows if r['query_key']==key];ct=collections.Counter(r['selection_status'] for r in rr)
    taxrows.append(dict(query_key=key,taxon_scope=u['query'],role=u['role'],direction=u['disease_direction'],budget_cap=u['candidate_budget_cap'],**{s:ct[s] for s in ['preselected','reserve','review','excluded','context_only','duplicate_assembly']},source_ids=';'.join(u['source_ids']),interpretation=u['interpretation']))
table(ROOT/'taxon_selection_summary.tsv',taxrows)
print(json.dumps(summary,ensure_ascii=False))
