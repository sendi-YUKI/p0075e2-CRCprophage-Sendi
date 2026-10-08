"""Conservative, scoped assembled-genome non-detection; no biological absence claim."""
import hashlib
import json
import math
from pathlib import Path
from mainline_common import read_table, write_table, sha, save, unique

METHOD = 'closed_gapless_catalog_target_nondetection_v1'
FIELDS = ['catalog_version','genome_id','votu_id','target_sha256','prepared_genome_sha256',
          'genome_sha256','carriage','value','callability','run_status','reason','rule_version',
          'policy_status','evidence_id','reviewer','review_date','assembly_evidence_sha256',
          'calibration_sha256','scope','limitations']

def signature(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def search_signature(rule):
    return signature({k:rule[k] for k in ('blast_evalue','word_size','edge_bp','rule_version')})

def bound_json(spec):
    p=Path(spec['path'])
    if not p.is_file() or sha(p)!=spec['sha256']:
        raise ValueError('Missing/changed scoped method evidence: '+str(p))
    return json.loads(p.read_text())

def validate_policy(p):
    if p.get('method')!=METHOD or p.get('status') not in ('proposed','synthetic','frozen'):
        raise ValueError('Unsupported whole-genome callability policy')
    if not p.get('rule_version'):raise ValueError('Callability rule version required')
    if not 0<=float(p['completeness_min'])<=100 or not 0<=float(p['contamination_max'])<=100:
        raise ValueError('Invalid callability quality bounds')
    return p

def quality_in_scope(qc, policy):
    if not policy or qc.get('status') != 'completed':
        return False
    try:
        values = [qc['completeness'], qc['contamination']]
        if any(isinstance(v, bool) for v in values):
            return False
        completeness, contamination = map(float, values)
        return (all(math.isfinite(v) and 0 <= v <= 100 for v in (completeness, contamination))
                and completeness >= float(policy['completeness_min'])
                and contamination <= float(policy['contamination_max']))
    except (KeyError, TypeError, ValueError):
        return False

def decide(records, catalog, preps, rule, out):
    """Positive calls use qualified catalogue members; zeros need additional bound evidence."""
    cm=json.loads((catalog/'catalog_manifest.json').read_text())
    members=read_table(catalog/'vOTU_members.tsv')
    p=rule.get('absence_rule')
    assemblies={};cal=None
    if p:
        validate_policy(p)
        if p.get('assembly_assessments'):
            package=bound_json(p['assembly_assessments'])
            assemblies=unique(package['genomes'],'genome_id')
        if p.get('calibration'):
            cal=bound_json(p['calibration'])
            results=bound_json(cal['results'])
            target_set=signature(sorted({(r['votu_id'],r['target_sha256']) for r in records}))
            if cal.get('search_signature')!=search_signature(rule) or cal.get('target_set_signature')!=target_set:
                raise ValueError('Calibration is for another search/target scope')
            if cal.get('assembly_class')!='closed_gapless_all_replicons' or cal.get('policy_rule_version')!=p['rule_version']:
                raise ValueError('Calibration applicability differs')
            if results.get('input_kind')!=('synthetic' if p['status']=='synthetic' else 'independent_biological_truth'):
                raise ValueError('Calibration truth type cannot authorize this policy')
            if not results.get('truth_manifest_sha256') or not results.get('strata'):
                raise ValueError('Calibration lacks truth identity/strata')
            if not all(isinstance(results.get(k),int) and results[k]>0 for k in ('positive_controls','negative_controls')):
                raise ValueError('Calibration requires both positive and negative controls')
            if results.get('acceptance_criteria_met') is not True or not results.get('criteria'):
                raise ValueError('Calibration criteria pending')
    decisions=[]
    for r in records:
        gid,v=r['genome_id'],r['votu_id'];observed=[m for m in members if m['genome_id']==gid and m['votu_id']==v]
        qualified=[m for m in observed if m.get('complete_element_eligible')=='true']
        source=preps[gid]/'input_manifest.json'
        identity=json.loads(source.read_text()) if source.is_file() else {}
        if identity and (identity.get('genome_id')!=gid or identity.get('normalized_fasta_sha256')!=r['prepared_genome_sha256']):
            raise ValueError('Normalization/genome identity mismatch')
        state='unknown';why=[];a=assemblies.get(gid,{})
        if qualified and cm['policy'].get('status') in ('synthetic','frozen'):
            state='present';why=['qualified_complete_catalog_member']
        else:
            if observed:why.append('fragment_or_unqualified_member_observed')
            if int(r['reported_hsps']):why.append('reported_homology_requires_review')
            if int(r['edge_hits']):why.append('edge_homology')
            if int(r['non_acgt_bp']):why.append('assembly_has_ambiguous_bases')
            if not int(r['assembled_bp']):why.append('empty_assembly')
            if r['search_status']!='no_reported_homology':why.append('search_not_complete_negative')
            if not identity.get('input_sha256'):why.append('original_genome_identity_missing')
            if not p or p['status']=='proposed':why.append('policy_pending')
            if not a:why.append('closed_assembly_assessment_missing')
            else:
                if a.get('prepared_genome_sha256')!=r['prepared_genome_sha256'] or a.get('genome_sha256')!=identity.get('input_sha256'):
                    raise ValueError('Assembly assessment identity mismatch')
                closure=bound_json(a['closure_evidence'])
                if (closure.get('genome_sha256')!=identity.get('input_sha256') or
                    closure.get('contig_count')!=int(r['contig_count']) or
                    closure.get('assembly_class')!='closed_gapless_all_replicons' or
                    not closure.get('evidence_sources') or not closure.get('reviewer') or not closure.get('review_date')):
                    why.append('all_replicons_closure_not_supported')
                qc=bound_json(a['quality_evidence'])
                if qc.get('genome_sha256')!=identity.get('input_sha256'):
                    raise ValueError('Quality evidence genome identity mismatch')
                if not quality_in_scope(qc, p):
                    why.append('quality_outside_policy_scope')
            if not cal or cal.get('status')!=('synthetic' if p and p['status']=='synthetic' else 'accepted') or not all(cal.get(k) for k in ('reviewer','review_date','decision_ref')):
                why.append('scoped_detection_calibration_pending')
            if p and p['status']=='frozen' and rule['status']!='frozen':why.append('search_not_frozen')
            if not why:state='absent_assessable';why=['no_reported_homology_with_scoped_closed_assembly_detection_evidence']
        row={k:r[k] for k in ('catalog_version','genome_id','votu_id','target_sha256','prepared_genome_sha256')}
        row.update(genome_sha256=identity.get('input_sha256',''),carriage=state,value=1 if state=='present' else 0 if state=='absent_assessable' else 'NA',
            callability='assessable' if state!='unknown' else 'unassessable',run_status='success',reason=';'.join(why),
            rule_version=p['rule_version'] if p else rule['rule_version'],policy_status=p['status'] if p else 'proposed',
            reviewer=cal.get('reviewer','') if cal else '',review_date=cal.get('review_date','') if cal else '',
            assembly_evidence_sha256=signature(a) if a else '',calibration_sha256=p['calibration']['sha256'] if p and p.get('calibration') else '',
            scope='whole_genome_target',limitations='Catalog-target non-detection in declared closed assembly/search domain; not biological absence or absence of all prophages')
        row['evidence_id']='wg_'+signature(row)[:32];decisions.append(row)
    write_table(out/'genome_callability.tsv',decisions,FIELDS)
    save(out/'callability_policy_receipt.json',dict(method=METHOD,policy=p,search_signature=search_signature(rule),
        decision_table_sha256=sha(out/'genome_callability.tsv'),original_identity_manifests={g:sha(d/'input_manifest.json') for g,d in preps.items() if (d/'input_manifest.json').is_file()},
        engineering_consumer_implemented=True,scientific_calibration_performed_this_run=False,
        evidence_truth_independently_verified=False,formal_research_authorized=False))
    return decisions

def load_decisions(directory, catalog):
    root=Path(directory).resolve();m=json.loads((root/'whole_genome_evidence_manifest.json').read_text())
    if m.get('status')!='completed' or m['catalog_manifest_sha256']!=sha(Path(catalog)/'catalog_manifest.json'):
        raise ValueError('Whole-genome evidence incomplete or stale catalogue')
    if not m.get('outputs'):raise ValueError('Missing whole-genome evidence inventory')
    for name,h in m['outputs'].items():
        p=(root/name).resolve()
        if not p.is_relative_to(root) or not p.is_file() or sha(p)!=h:raise ValueError('Whole-genome evidence changed')
    rows=read_table(root/'genome_callability.tsv',FIELDS)
    if len({(r['genome_id'],r['votu_id']) for r in rows})!=len(rows):raise ValueError('Duplicate genome target decisions')
    return rows,m
