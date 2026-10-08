#!/usr/bin/env python3
"""Opt-in CRC-source consumer. Metadata check never runs tools or reads sequences."""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(os.environ.get('CRC_SPARK_REPO', '/srv/CRC-PHIRE/projects/crc-pipeline')).resolve()
sys.path.insert(0, str(ROOT / 'bin'))
from mainline_common import sha, save, read_table, write_table, unique, fasta, new_output, safe_id
import study_contract as study
import spark_runtime

MODEL = 'crc_source_logistic_mds_v1'
KEY = ['catalog_version', 'species', 'votu_id', 'comparison_tier']
RESULT = KEY + ['model_id', 'status', 'reason', 'n', 'cases', 'controls', 'carriers', 'noncarriers',
    'logOR', 'OR', 'SE', 'CI_lower', 'CI_upper', 'ci_method', 'p_raw', 'formula', 'reference_levels',
    'dropped_constant_terms', 'config_sha256', 'structure_sha256', 'R_version']
SUBJECT = ['genome_id', 'participant_id', 'isolate_id', 'study_id', 'species', 'body_site',
    'CRC_source', 'selection_rank', 'lineage', 'lineage_rule_version', 'covariates_json']


def small_json(path):
    return study.strict_json(study.read_small(path))


def config(path, execute=False, production=False):
    from jsonschema import Draft202012Validator
    c = small_json(path)
    schema = small_json(ROOT / 'assets/schemas/host_source.schema.json')
    errors = list(Draft202012Validator(schema).iter_errors(c))
    if errors:
        raise ValueError('Host-source config: ' + '; '.join(e.message for e in errors[:5]))
    if c['status'] == 'frozen' and not all(c[k] for k in ('approved_by', 'approved_on', 'decision_ref')):
        raise ValueError('Frozen host-source choices require approval/date/decision provenance')
    pending = [k for k in ('species', 'catalog_version', 'min_pair_sites', 'min_pair_fraction', 'mds_k',
        'negative_eigen_policy', 'min_carriers', 'min_noncarriers', 'body_site', 'selection_rule_version') if c[k] is None]
    if execute and (pending or c['status'] == 'proposed'):
        raise ValueError('Host-source scientific configuration pending: ' + ','.join(pending))
    if production and c['status'] != 'frozen':
        raise ValueError('Production host-source config must be frozen; synthetic is fixture-only')
    names = [r['name'] for r in c['covariates']]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate covariate')
    for r in c['covariates']:
        if r['kind'] == 'categorical' and (len(r['levels']) < 2 or len(set(r['levels'])) != len(r['levels'])):
            raise ValueError('Categorical levels must be explicitly ordered and unique')
    return c, pending


def bundle(path):
    p = Path(path).resolve(); m = small_json(p / 'manifest.json')
    if m.get('schema_version') != 1 or m.get('purpose') != 'host_source_metadata':
        raise ValueError('Unknown host-source bundle')
    files = {}
    for row in m['files']:
        name = row['path']; f = (p / name).resolve()
        if Path(name).name != name or not f.is_relative_to(p) or name in files:
            raise ValueError('Unsafe/duplicate metadata member')
        study.read_small(f)
        if sha(f) != row['sha256']:
            raise ValueError('Changed metadata ' + name)
        files[name] = f
    if not {'subjects.tsv', 'absence_reviews.tsv', 'family.tsv', 'recombination.json'} <= set(files):
        raise ValueError('Incomplete host-source metadata bundle')
    contract = study.metadata_path(m['study_contract'], p) if m.get('study_contract') else None
    if contract and sha(contract) != m.get('study_contract_sha256'):
        raise ValueError('Changed clinical contract')
    return p, m, files, contract


def recombination(c, files, alignment_hash=None):
    r = small_json(files['recombination.json'])
    if r.get('status') not in ('not_assessed', 'reviewed'):
        raise ValueError('Unknown recombination evidence state')
    if r['status'] == 'reviewed':
        required = ('method', 'method_version', 'applicability', 'reviewer', 'review_date', 'alignment_sha256', 'evidence_file', 'evidence_sha256')
        if not all(r.get(k) for k in required):
            raise ValueError('Recombination review lacks method/applicability/source')
        if r['evidence_file'] not in files or sha(files[r['evidence_file']]) != r['evidence_sha256']:
            raise ValueError('Recombination source not bound to bundle')
        if alignment_hash and r['alignment_sha256'] != alignment_hash:
            raise ValueError('Recombination evidence applies to another alignment')
    return r


def clinical_pending(cfg):
    reasons=[]
    if cfg['status']!='frozen':reasons.append('clinical_contract_not_frozen')
    if cfg['stage']!='within_species':reasons.append('clinical_stage_not_within_species')
    if not cfg['master_universe_version']:reasons.append('master_universe_not_versioned')
    if not cfg['source_documents'] or not cfg['identity_references']:reasons.append('clinical_source_or_identity_provenance_missing')
    for key,value in dict(outcome='CRC',exposure='feature',unit='participant',missing_policy='exclude_unknown_no_zero',correction='BH').items():
        if cfg['statistics'][key]!=value:reasons.append('clinical_model_declaration_mismatch:'+key)
    return reasons


def inspect(path, metadata, production=False, planned_recombination=False):
    c, pending = config(path, production=production)
    result = dict(status='configuration_' + c['status'], pending=pending, analysis_started=False,
        whole_genome_search_evidence_producer='implemented_opt_in',whole_genome_absence_decision='not_calibrated', scientific_calibration=False)
    if not metadata:
        result['pending'].append('host_source_bundle'); return result
    p, m, files, contract = bundle(metadata)
    rec = recombination(c, files)
    result['recombination_status'] = rec['status']
    if rec['status'] != 'reviewed' and not planned_recombination: result['pending'].append('recombination_control')
    if planned_recombination:result['recombination_status']='planned_frozen_producer_runtime_alignment_check_required'
    if c['status']=='frozen' and m.get('status')!='frozen': result['pending'].append('metadata_bundle_not_frozen')
    if not contract: result['pending'].append('clinical_contract')
    else:
        cfg, tables, evidence, protected = study.load_contract(contract)
        result['pending'].extend(clinical_pending(cfg))
        errs, warns = study.relationships(cfg, tables)
        # Research-generated observation/feature FKs are resolved at runtime, not fabricated here.
        result['clinical_contract_schema_valid'] = True
        result['clinical_relationships_before_research_merge'] = errs
    result['declared_subjects'] = len(read_table(files['subjects.tsv']))
    result['family_hypotheses'] = len(family_rows(files['family.tsv']))
    if not result['declared_subjects']: result['pending'].append('subjects')
    if not result['family_hypotheses']: result['pending'].append('family_manifest')
    if production and result['pending']:
        raise ValueError('Host-source preparation incomplete: ' + ','.join(result['pending']))
    for key in ('python', 'r-sensitivity' if any(c.get('sensitivity',{}).get(k,False) for k in ('lmm','firth_on_separation')) else 'association'):
        _, item, prefix, lock, _ = spark_runtime.environment_definition(key)
        result.setdefault('environments', {})[key] = dict(prefix=str(prefix), lock_sha256=sha(lock), manifest_sha256=item['manifest_sha256'])
    return result


def verify_artifact(directory, manifest):
    p = Path(directory).resolve(); m = small_json(p / manifest)
    for name, expected in m['outputs'].items():
        f = p / name
        if Path(name).name != name or not f.is_file() or sha(f) != expected:
            raise ValueError('Changed/unsafe artifact: ' + name)
    return m


def verify_artifact_recursive(directory, manifest):
    root=Path(directory).resolve();m=small_json(root/manifest)
    if not m.get('outputs'):raise ValueError('Missing artifact inventory')
    for name,h in m['outputs'].items():
        p=(root/name).resolve()
        if not p.is_relative_to(root) or not p.is_file() or sha(p)!=h:raise ValueError('Changed recursive artifact '+name)
    return m

def proof(e, kind, ident, field, value):
    return (e.get('entity_type'), e.get('entity_id'), e.get('field'), e.get('value')) == (kind, ident, field, value) and e.get('level') in ('A', 'B') and e.get('review_state') == 'accepted' and bool(e.get('reviewer')) and bool(e.get('review_date'))


def merged_contract(contract, research, c, files):
    cfg, tables, evidence, protected = study.load_contract(contract)
    if c['status']=='frozen' and clinical_pending(cfg):
        raise ValueError('Clinical research configuration pending: '+','.join(clinical_pending(cfg)))
    if c['status'] == 'synthetic' and cfg['stage'] != 'preparation':
        raise ValueError('Synthetic fixtures require preparation contract')
    schema = small_json(study.SCHEMA_PATH)['$defs']['rows']
    for name in ('features', 'feature_links', 'observations', 'evidence'):
        old = {r[study.KEYS[name]]: r for r in tables[name]}
        for raw in read_table(Path(research) / ('research_' + name + '.tsv'), schema[name]['properties']):
            row = {k: study.row_value(raw[k], spec) for k, spec in schema[name]['properties'].items()}
            from jsonschema import Draft202012Validator, FormatChecker
            study.validate_schema(Draft202012Validator(schema[name], format_checker=FormatChecker()), row, name)
            key = row[study.KEYS[name]]
            if key in old and old[key] != row: raise ValueError('Conflicting clinical/research row ' + key)
            old[key] = row
        tables[name] = list(old.values())
    ev = unique(tables['evidence'], 'evidence_id'); gs = unique(tables['genomes'], 'genome_id')
    obs = unique(tables['observations'], 'observation_id'); seen = set()
    for r in read_table(files['absence_reviews.tsv']):
        oid = r['observation_id']
        if oid in seen: raise ValueError('Duplicate absence review')
        seen.add(oid); o = obs.get(oid, {}); g = gs.get(r['genome_id'], {})
        expected = (r['genome_id'], r['votu_id'], 'votu', 'genome')
        if (o.get('entity_id'), o.get('feature_id'), o.get('feature_unit'), o.get('entity_type')) != expected:
            raise ValueError('Absence endpoint/identity mismatch')
        if r['scope'] != 'whole_genome_target' or r['catalog_version'] != c['catalog_version'] or not r['rule_version']:
            raise ValueError('Local/unversioned absence cannot become genome 0')
        if not g.get('fasta_sha256') or g['fasta_sha256'] != r['genome_sha256']:
            raise ValueError('Absence genome identity not bound')
        if not proof(ev.get(r['evidence_id'], {}), 'observation', oid, 'callability', 'assessable'):
            raise ValueError('Absence lacks matching reviewed A/B evidence')
        if r['source_file'] not in files or sha(files[r['source_file']]) != r['source_sha256']:
            raise ValueError('Absence source hash not in bundle')
        raw = small_json(files[r['source_file']])
        claims = raw.get('claims', [])
        expected_claim = {k: r[k] for k in ('observation_id','genome_id','votu_id','catalog_version','scope','rule_version','genome_sha256')}
        matches = [x for x in claims if all(x.get(k) == v for k, v in expected_claim.items())]
        if len(matches) != 1 or not all(matches[0].get(k) for k in ('method', 'coverage_assessment', 'limitations')):
            raise ValueError('Absence source lacks endpoint-specific reviewed assessment')
        if o['run_status'] != 'success' or o['carriage'] == 'present':
            raise ValueError('Cannot override positive/failed observation with absence')
        o.update(carriage='absent_assessable', value=0, callability='assessable',
            callability_evidence_id=r['evidence_id'], rule_version=r['rule_version'])
    generated=Path(research)/'genome_callability.tsv'
    if generated.is_file():
        manifest=verify_artifact(research,'research_summary.json')
        if 'genome_callability.tsv' not in manifest['outputs']:raise ValueError('Callability not in research inventory')
        for r in read_table(generated):
            if r['carriage']!='absent_assessable':continue
            matching=[o for o in tables['observations'] if o['entity_type']=='genome' and o['entity_id']==r['genome_id'] and o['feature_unit']=='votu' and o['feature_id']==r['votu_id']]
            if len(matching)!=1:raise ValueError('Generated callability endpoint mismatch')
            o=matching[0];oid=o['observation_id']
            if oid in seen:raise ValueError('Duplicate manual/generated negative source')
            if r['catalog_version']!=c['catalog_version'] or r['scope']!='whole_genome_target' or r['run_status']!='success':raise ValueError('Wrong generated absence scope')
            if r['genome_sha256']!=gs.get(r['genome_id'],{}).get('fasta_sha256'):raise ValueError('Generated absence source genome differs')
            if r['policy_status']!=c['status'] or not r['calibration_sha256'] or not r['assembly_evidence_sha256']:raise ValueError('Generated absence policy/evidence pending')
            if not proof(ev.get(r['evidence_id'],{}),'observation',oid,'callability','assessable'):raise ValueError('Generated absence lacks scoped review')
            if o['carriage']!='absent_assessable' or o['value']!=0 or o['rule_version']!=r['rule_version']:raise ValueError('Generated observation mismatch')
            seen.add(oid)
    for o in tables['observations']:
        if o['feature_unit']=='votu' and o['carriage']=='absent_assessable' and o['observation_id'] not in seen:
            raise ValueError('Whole-genome absence lacks scoped hash-bound import record')
    errors, warnings = study.relationships(cfg, tables)
    errors += study.feature_relationships(tables)
    if errors: raise ValueError('Merged study contract: ' + json.dumps(errors[:8]))
    return cfg, tables, dict(source_contract_sha256=sha(contract), evidence=evidence,
        warnings=warnings, absence_imports=len(seen), absence_evidence_truth_verified=False,
        absence_producer_implemented=generated.is_file())


def select_subjects(t, c, files, out):
    gs = unique(t['genomes'], 'genome_id'); iso = unique(t['isolates'], 'isolate_id')
    people = unique(t['participants'], 'participant_id'); studies = unique(t['studies'], 'study_id')
    ev = unique(t['evidence'], 'evidence_id'); selected = []; excluded = []; candidates = []
    for r in read_table(files['subjects.tsv'], ['genome_id','species','selection_rank','covariates_json','lineage','lineage_rule_version']):
        g = gs.get(r['genome_id']); reason = ''
        if g is None: raise ValueError('Subject genome absent from study contract')
        safe_id(r['genome_id']); i = iso.get(g['isolate_id'], {}); p = people.get(i.get('participant_id'), {})
        s = studies.get(i.get('study_id'), {})
        if not g['master_member'] or g['role'] != 'discovery': reason = 'not_discovery_master'
        elif r['species'] != c['species']: reason = 'outside_declared_species'
        elif not p or not proof(ev.get(p.get('identity_evidence_id'), {}), 'participant', p['participant_id'], 'identity', p['participant_id']): reason = 'independence_unconfirmed'
        elif s.get('role') != 'development': reason = 'not_development'
        elif i['body_site'] != c['body_site']: reason = 'outside_declared_body_site'
        rank = r['selection_rank']
        if not re.fullmatch(r'[0-9]+', rank): raise ValueError('Selection rank must be a predeclared nonnegative integer')
        row = dict(genome_id=r['genome_id'], participant_id=i.get('participant_id') or '', isolate_id=g['isolate_id'] or '',
            study_id=i.get('study_id') or '', species=r['species'], body_site=i.get('body_site') or '',
            CRC_source=1 if i.get('disease_status') == 'CRC_confirmed' else 0 if i.get('disease_status') == 'control_qualified' else '',
            selection_rank=int(rank), lineage=r['lineage'], lineage_rule_version=r['lineage_rule_version'], covariates_json=r['covariates_json'])
        if row['lineage'] and not row['lineage_rule_version']: raise ValueError('Lineage without declared rule')
        if reason: excluded.append(dict(row, reason=reason))
        else: candidates.append(row)
    unique(candidates + excluded, 'genome_id')
    # Choose representatives before disease eligibility/covariates/features. Never switch after exclusion.
    groups = collections.defaultdict(list)
    for row in candidates: groups[row['participant_id']].append(row)
    for pid, rows in sorted(groups.items()):
        ordered = sorted(rows, key=lambda r:(r['selection_rank'], r['genome_id']))
        chosen = ordered[0]
        for other in ordered[1:]: excluded.append(dict(other, reason='prespecified_other_assembly'))
        gid = chosen['genome_id']; reason = ''
        qualifying = [e for e in t['eligibility'] if e['entity_type'] == 'genome' and e['entity_id'] == gid and e['endpoint'] == 'disease']
        if len(qualifying) != 1 or qualifying[0]['status'] != 'eligible': reason = 'disease_eligibility_missing'
        elif not proof(ev.get(qualifying[0]['evidence_id'], {}), 'genome', gid, 'eligibility:disease', 'eligible'): reason = 'eligibility_review_missing'
        elif chosen['CRC_source'] == '': reason = 'clinical_label_not_qualified'
        try:
            cov = study.strict_json(chosen['covariates_json'])
            if not isinstance(cov, dict) or set(cov) != {v['name'] for v in c['covariates']}: raise ValueError('Covariate names differ')
            for spec in c['covariates']:
                v = cov[spec['name']]
                if v is None: reason = reason or 'missing_covariate'
                elif spec['kind'] == 'numeric' and (isinstance(v, bool) or not isinstance(v, (int,float)) or not math.isfinite(v)): raise ValueError('Invalid numeric covariate')
                elif spec['kind'] == 'categorical' and v not in spec['levels']: raise ValueError('Undeclared factor level')
        except (ValueError, TypeError) as exc: raise ValueError('Subject covariates: ' + str(exc))
        if reason: excluded.append(dict(chosen, reason=reason))
        else: selected.append(chosen)
    selected.sort(key=lambda r:r['genome_id'])
    write_table(out/'subjects.tsv', selected, SUBJECT)
    write_table(out/'subject_exclusions.tsv', excluded, SUBJECT+['reason'])
    save(out/'subject_selection.json',dict(declared_genomes=len(selected)+len(excluded),
        selected_independent_subjects=len(selected),declared_participants_before_exclusion=len({r['participant_id'] for r in selected+excluded if r['participant_id']}),
        exclusion_counts=dict(collections.Counter(r['reason'] for r in excluded)),selection_rule_version=c['selection_rule_version'],
        selection_order='declared rank before clinical eligibility, covariate and feature exclusion; no replacement',
        subjects_sha256=sha(files['subjects.tsv']),clinical_truth_verified=False))
    support=collections.Counter((r['study_id'], r['CRC_source'], r['body_site'], r['lineage'] or 'unknown') for r in selected)
    write_table(out/'support.tsv', [dict(study_id=k[0],CRC_source=k[1],body_site=k[2],lineage=k[3],independent_participants=v) for k,v in sorted(support.items())], ['study_id','CRC_source','body_site','lineage','independent_participants'])
    return selected


def resource_estimate(n, length, budget):
    # Conservative engineering ceiling, not a measured scalability guarantee.
    estimate = 512*1024**2 + 128*n*n + 8*n*length
    if estimate > budget*1024**3: raise ValueError('Estimated matrices/copies/alignment exceed declared memory budget')
    return dict(N=n,L=length,one_dense_matrix_bytes=8*n*n,estimated_peak_bytes=estimate,
        memory_budget_gib=budget,performance_validated=False,algorithm='pair loop O(N^2*L), dense MDS O(N^3)')


def distances(mask, tree, subjects, c, out):
    mask=Path(mask); tree=Path(tree)
    m=small_json(mask/'mask_manifest.json'); ts=small_json(tree/'tree_status.json')
    alignment_hash=sha(mask/'alignment.fna')
    if ts.get('input_alignment_sha256') != alignment_hash:
        if ts.get('alignment_relationship')!='Gubbins_branch_specific_inference_vs_global_union_distance_mask' or ts.get('distance_alignment_sha256')!=alignment_hash:
            raise ValueError('Final tree belongs to another alignment')
        parent=tree.parent
        rm=verify_artifact_recursive(parent,'recombination_manifest.json')
        if sha(parent/'gubbins_input.fna')!=ts['input_alignment_sha256'] or rm.get('applicability',{}).get('status')!='passed':
            raise ValueError('Gubbins tree/distance source relationship unverified')
    tips=unique(read_table(mask/'tree_tips.tsv'), 'genome_id')
    finaltips=unique(read_table(tree/'tree_tips.tsv'), 'genome_id')
    if tips != finaltips or any(r['genome_id']!=r['tip_id'] or r['genome_id']!=r['source_host_genome'] for r in tips.values()):
        raise ValueError('Host tip/source identity mismatch')
    # Preflight input size before loading: file upper bound on string payload, with allocation ceiling.
    if (mask/'alignment.fna').stat().st_size*8 > c['memory_gib']*1024**3: raise ValueError('Alignment input exceeds memory budget')
    seq=fasta(mask/'alignment.fna')
    if set(seq)!=set(tips): raise ValueError('Alignment/tip identities differ')
    lens={len(s) for s in seq.values()}
    if len(lens)>1: raise ValueError('Unequal alignment lengths')
    length=next(iter(lens),0); ids=[s['genome_id'] for s in subjects]
    if not set(ids)<=set(seq): raise ValueError('Selected independent subject missing in alignment')
    synthetic_cap=240000 if ts.get('alignment_relationship')=='Gubbins_branch_specific_inference_vs_global_union_distance_mask' else 20000
    if c['status']=='synthetic' and (sum(map(len,seq.values()))>synthetic_cap or len(ids)>100): raise ValueError('Synthetic input cap exceeded')
    if c.get('sensitivity',{}).get('lmm'):
        if len(ids)*length*96 + len(ids)**2*128 > c['memory_gib']*1024**3:
            raise ValueError('Kinship input/temporary memory estimate exceeds allocation')
        write_table(out/'kinship_input.tsv',[dict(genome_id=g,sequence=seq[g]) for g in ids],['genome_id','sequence'])
    res=resource_estimate(len(seq),length,c['memory_gib']); pairs=[]; matrix={g:{} for g in ids}; missing=False
    for i,g in enumerate(ids):
        matrix[g][g]=0
        for h in ids[i+1:]:
            comparable=0; mismatch=0
            for a,b in zip(seq[g],seq[h]):
                if a in 'ACGT' and b in 'ACGT': comparable+=1; mismatch+=a!=b
            ok=length>0 and comparable>=c['min_pair_sites'] and comparable/length>=c['min_pair_fraction']
            value=mismatch/comparable if ok and comparable else None
            missing |= value is None
            matrix[g][h]=matrix[h][g]=value
            pairs.append(dict(genome_a=g,genome_b=h,comparable_sites=comparable,mismatches=mismatch,
                total_columns=length,distance=value,reason='' if value is not None else 'insufficient_comparable_sites'))
    write_table(out/'pair_denominators.tsv',pairs,['genome_a','genome_b','comparable_sites','mismatches','total_columns','distance','reason'])
    write_table(out/'distance.tsv',[dict(genome_id=g,**matrix[g]) for g in ids],['genome_id']+ids)
    md=dict(method=c['distance_method'],alignment_sha256=alignment_hash,mask_manifest_sha256=sha(mask/'mask_manifest.json'),
        final_tree_manifest_sha256=sha(tree/'tree_status.json'),tree_tips_sha256=sha(mask/'tree_tips.tsv'),
        sample_order=ids,resource_estimate=res,missing_distances=missing,
        status='ready' if len(ids)>c['mds_k'] and not missing else 'not_estimable',
        reason='insufficient_subjects' if len(ids)<=c['mds_k'] else 'missing_distances' if missing else '',
        MDS_policy=c['mds_population_policy'],recombination_is_not_mobile_mask=True)
    save(out/'distance_manifest.json',md)
    return md


def family_rows(path):
    rows=read_table(path, KEY+['model_id','family_id']); seen=set()
    for r in rows:
        key=tuple(r[k] for k in KEY)
        if key in seen or any(not x for x in key) or r['model_id']!=MODEL or r['family_id']!='F1':
            raise ValueError('Invalid/duplicate F1 species x vOTU hypothesis')
        for k in KEY: safe_id(r[k])
        seen.add(key)
    return rows


def observations(t, research, subjects, hypotheses, c, out):
    obs={(r['entity_id'],r['feature_id']):r for r in t['observations'] if r['entity_type']=='genome'}
    features=unique(t['features'],'feature_id'); rows=[]; exclusions=[]; links=[]
    for h in hypotheses:
        f=features.get(h['votu_id'],{})
        if f.get('feature_unit')!='votu': raise ValueError('F1 target must be a catalog vOTU')
        for s in subjects:
            o=obs.get((s['genome_id'],h['votu_id']),{}); why=''
            if o.get('feature_unit')!='votu' or o.get('value_axis')!='carriage_binary': why='endpoint_missing_or_mismatched'
            elif o.get('run_status')!='success': why='not_success'
            elif o.get('callability')!='assessable' or o.get('carriage') not in ('present','absent_assessable') or o.get('value') not in (0,1): why='unknown_or_unassessable'
            if why: exclusions.append(dict(h,genome_id=s['genome_id'],participant_id=s['participant_id'],reason=why))
            else: rows.append(dict(h,**s,feature_carriage=o['value'],observation_id=o['observation_id']))
    cross=read_table(Path(research)/'locus_feature_crosswalk.tsv')
    unique(cross,'locus_feature_id')
    for x in cross:
        f=features.get(x['locus_feature_id'],{})
        o=obs.get((x['genome_id'],x['locus_feature_id']),{})
        if f.get('feature_unit')!='locus' or f.get('genome_id')!=x['genome_id'] or not x['locus_group']:
            raise ValueError('Invalid cross-genome locus crosswalk')
        links.append(dict(catalog_version=c['catalog_version'],locus_group=x['locus_group'],genome_id=x['genome_id'],
            locus_feature_id=x['locus_feature_id'],carriage=o.get('carriage'),value=o.get('value'),
            endpoint_scope='local_locus_only',included_in_F1=False))
    write_table(out/'model_input.tsv', rows, KEY+['model_id','family_id']+[k for k in SUBJECT if k!='species']+['feature_carriage','observation_id'])
    write_table(out/'feature_exclusions.tsv', exclusions, KEY+['genome_id','participant_id','reason'])
    write_table(out/'locus_observations.tsv',links,['catalog_version','locus_group','genome_id','locus_feature_id','carriage','value','endpoint_scope','included_in_F1'])
    return rows


def bh(values):
    order=sorted(range(len(values)),key=lambda i:values[i]); out=[None]*len(values); prev=1.0
    for rank,i in reversed(list(enumerate(order,1))):
        prev=min(prev,values[i]*len(values)/rank);out[i]=prev
    return out


def aggregate(family, sources, out, allow_incomplete=False):
    expected=family_rows(family); wanted={tuple(r[k] for k in KEY):r for r in expected}; results={}; hashes=[];scopes=set()
    for source in sources:
        source=Path(source); m=verify_artifact(source,'host_source_manifest.json')
        scopes.add(m['scope'])
        if len(scopes)>1 or m['scope'] not in ('synthetic','frozen'):raise ValueError('Cannot mix synthetic and formal F1 results')
        if m['family_manifest_sha256']!=sha(family): raise ValueError('Different frozen family manifest')
        hashes.append(dict(path=str(source),model_results_sha256=sha(source/'model_results.tsv'),config_sha256=m['config_sha256']))
        for r in read_table(source/'model_results.tsv',RESULT):
            k=tuple(r[x] for x in KEY)
            if k in results or k not in wanted or r['model_id']!=MODEL: raise ValueError('Duplicate/outside-F1 result')
            if r['status'] not in ('estimated','not_estimable','not_identifiable'):raise ValueError('Unknown model result state')
            results[k]=r
    missing=set(wanted)-set(results)
    if missing and not allow_incomplete: raise ValueError('Missing model rows; refuse to shrink F1 denominator')
    rows=[]; ps=[]
    for k,h in wanted.items():
        r=results.get(k,dict(h,status='not_received',reason='missing_species_result',p_raw=''))
        if r['status']=='estimated':
            try:p=float(r['p_raw'])
            except (ValueError,TypeError):raise ValueError('Estimated result lacks numeric P')
            if not math.isfinite(p) or not 0<=p<=1:raise ValueError('Invalid P')
        else:
            if r.get('p_raw') not in ('','NA',None):raise ValueError('Failed model cannot carry a P')
            p=1.0
        rows.append(dict(h,status=r['status'],reason=r['reason'],p_raw=r.get('p_raw',''),p_for_correction=p));ps.append(p)
    # BH within comparison_tier, never across. Pooling tiers would let a tier3
    # comparison (niche confounded with disease) share a denominator with a
    # tier1 one. See protocol section 9.2b and section 9.4 F1.
    if missing:
        for r in rows:
            r['q_BH']=''
    else:
        untiered=[r for r in rows if not str(r.get('comparison_tier','')).strip()]
        if untiered:
            raise ValueError('%d F1 rows carry no comparison_tier; refusing to correct, '
                             'because a default tier would silently pool strata'%len(untiered))
        by_tier={}
        for i,r in enumerate(rows):
            by_tier.setdefault(r['comparison_tier'],[]).append(i)
        for tier,idx in by_tier.items():
            tq=bh([rows[i]['p_for_correction'] for i in idx])
            for i,q in zip(idx,tq):
                rows[i]['q_BH']=q
        tier_sizes={t:len(i) for t,i in sorted(by_tier.items())}
    write_table(out/'F1_results.tsv',rows,KEY+['model_id','family_id','status','reason','p_raw','p_for_correction','q_BH'])
    save(out/'F1_manifest.json',dict(status='pending_other_species' if missing else 'complete' if rows else 'empty_family',
        family_sha256=sha(family),family_size=len(rows),missing_rows=len(missing),policy='retain_failed_with_p1',
        scope=next(iter(scopes),'not_assessed'),
        sources=hashes,correction='BH_within_comparison_tier',
        tier_sizes=(tier_sizes if not missing else {}),
        tier_policy='tiers never pooled; no selecting the most significant tier; tier frozen before disease results (protocol 9.2b)',
        clinical_truth_verified=False,scientific_calibration=False,
        outputs={'F1_results.tsv':sha(out/'F1_results.tsv')}))


def execute(a):
    if not os.environ.get('SLURM_JOB_ID'): raise ValueError('Host-source computation requires controlled Slurm allocation')
    c,_=config(a.config,execute=True)
    allocation=spark_runtime.allocation(require_lock=True)
    if c['memory_gib']*1024>allocation['memory_mib']:raise ValueError('Declared memory exceeds actual allocation')
    destination=Path(a.outdir).resolve()
    for raw in [a.bundle,a.research,a.mask,a.tree,*([a.lineages] if getattr(a,'lineages',None) else []),*([a.recombination_result] if getattr(a,'recombination_result',None) else [])]:
        protected=Path(raw).resolve()
        if destination==protected or destination.is_relative_to(protected) or protected.is_relative_to(destination):
            raise ValueError('Output overlaps an input artifact')
    p,m,files,contract=bundle(a.bundle)
    if not contract: raise ValueError('Clinical metadata contract required')
    if c['status']=='frozen' and m.get('status')!='frozen': raise ValueError('Production metadata bundle must be frozen')
    research=Path(a.research); rm=verify_artifact(research,'research_summary.json')
    if rm['catalog_version']!=c['catalog_version']:raise ValueError('Research/config catalog mismatch')
    if c['status']=='frozen' and rm.get('policy_status')!='frozen':raise ValueError('Production research catalog not frozen')
    cfg,t,clinical=merged_contract(contract,research,c,files)
    out=new_output(a.outdir)
    for subset in c.get('sensitivity_subsets',[]):
        name=subset['evidence_file']
        if name not in files or sha(files[name])!=subset['evidence_sha256']:
            raise ValueError('Sensitivity subset lacks hash-bound prespecified evidence')
        evidence=small_json(files[name])
        if evidence.get('subset_name')!=subset['name'] or set(evidence.get('genome_ids',[]))!=set(subset['genome_ids']):
            raise ValueError('Sensitivity subset evidence identity mismatch')
        if c['status']=='frozen' and not all(evidence.get(k) for k in ('reviewer','date','rule_version')):
            raise ValueError('Formal sensitivity subset requires rule review')
    save(out/'config.json',c);save(out/'clinical_merge.json',clinical)
    subjects=select_subjects(t,c,files,out)
    if getattr(a,'lineages',None):
        lr=Path(a.lineages);lm=verify_artifact(lr,'lineage_manifest.json')
        if c['status']=='frozen' and lm.get('config_status')!='frozen':raise ValueError('Formal model requires frozen lineage rule')
        # Require the exact tree used for the structure result, not a similarly named tree.
        tree_files=[f for f in Path(a.tree).rglob('*') if f.is_file() and f.suffix in ('.nwk','.tree','.treefile','.tre')]
        if not any(sha(f)==lm.get('tree_sha256') for f in tree_files):
            raise ValueError('Lineage assignment does not match the supplied final tree')
        assignments=unique(read_table(lr/'lineages.tsv'),'genome_id')
        if not {r['genome_id'] for r in subjects}<=set(assignments):raise ValueError('Lineage assignments lack model subjects')
        for row in subjects:
            v=assignments[row['genome_id']]
            if row['lineage'] not in ('','unknown',v['lineage']):raise ValueError('Declared and produced lineage conflict')
            row.update(lineage=v['lineage'],lineage_rule_version=v['lineage_rule_version'])
        write_table(out/'subjects.tsv',subjects,SUBJECT)
        save(out/'lineage_input.json',dict(manifest_sha256=sha(lr/'lineage_manifest.json'),tree_sha256=lm['tree_sha256']))

    hypotheses=family_rows(files['family.tsv'])
    ours=[h for h in hypotheses if h['species']==c['species'] and h['catalog_version']==c['catalog_version']]
    if not ours:raise ValueError('No declared F1 hypotheses for current species/catalog')
    write_table(out/'hypotheses.tsv',ours,KEY+['model_id','family_id'])
    inputs=observations(t,research,subjects,ours,c,out)
    md=distances(a.mask,a.tree,subjects,c,out)
    if getattr(a,'recombination_result',None):
        rr=Path(a.recombination_result);rm=verify_artifact_recursive(rr,'recombination_manifest.json')
        if rm.get('status')!='completed' or rm.get('scope')!='declared_recent_common_ancestor_group':
            raise ValueError('Unqualified recombination producer')
        if sha(rr/'masked/alignment.fna')!=md['alignment_sha256']:
            raise ValueError('Recombination result does not match model alignment')
        review=rm.get('applicability_review',{})
        ready=rm.get('config_status')=='frozen' and rm.get('applicability',{}).get('status')=='passed' and bool(review.get('evidence')) and all(review.get(k) for k in ('reviewer','date','decision_ref'))
        rec=dict(status='reviewed' if ready else 'not_assessed',method=rm['method'],applicability=rm['scope'],
            alignment_sha256=md['alignment_sha256'],evidence_sha256=sha(rr/'recombination_manifest.json'),review=review,
            scientific_calibration=False)
    else:rec=recombination(c,files,md['alignment_sha256'])
    if c['status']=='frozen' and rec['status']!='reviewed':raise ValueError('Formal model requires reviewed recombination control')
    save(out/'recombination.json',rec)
    _,item,prefix,lock,activation=spark_runtime.environment_definition('r-sensitivity' if any(c.get('sensitivity',{}).get(k,False) for k in ('lmm','firth_on_separation')) else 'association')
    env=dict(os.environ,**activation);env.update(PATH=str(prefix/'bin')+':'+env.get('PATH',''),CONDA_PREFIX=str(prefix),
        R_HOME=str(prefix/'lib/R'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    if any(c.get('sensitivity',{}).get(k,False) for k in ('lmm','firth_on_separation')):
        spark_runtime.stable_identity('r-sensitivity',item,spark_runtime.read(spark_runtime.DEPENDENCIES))
    env['CRC_HOST_SENSITIVITY_SCRIPT']=str(Path(a.sensitivity_script).resolve())
    command=[str(prefix/'bin/Rscript'),'--vanilla',str(Path(a.r_script).resolve()),str(out.resolve())]
    with (out/'R.stdout.log').open('w') as stdout,(out/'R.stderr.log').open('w') as stderr:
        try:result=subprocess.run(command,env=env,stdout=stdout,stderr=stderr,timeout=a.timeout)
        except subprocess.TimeoutExpired:
            save(out/'R_command.json',dict(argv=command,exit_code=124,timeout_seconds=a.timeout));raise
    save(out/'R_command.json',dict(argv=command,exit_code=result.returncode,timeout_seconds=a.timeout))
    if result.returncode:raise RuntimeError('CRC-source R consumer failed: '+str(result.returncode))
    results=read_table(out/'model_results.tsv',RESULT)
    if {tuple(r[k] for k in KEY) for r in results}!={tuple(r[k] for k in KEY) for r in ours} or len(results)!=len(ours):raise ValueError('Missing/duplicate model result')
    for r in results:r.update(config_sha256=sha(a.config),structure_sha256=sha(out/'distance_manifest.json'))
    write_table(out/'model_results.tsv',results,RESULT)
    record=dict(schema_version=1,status='completed',model_id=MODEL,catalog_version=c['catalog_version'],species=c['species'],
        config_sha256=sha(a.config),research_manifest_sha256=sha(research/'research_summary.json'),bundle_manifest_sha256=sha(p/'manifest.json'),
        family_manifest_sha256=sha(files['family.tsv']),subject_count=len(subjects),model_input_rows=len(inputs),hypotheses=len(results),
        estimated=sum(r['status']=='estimated' for r in results),scope=c['status'],recombination_status=rec['status'],
        sensitivity_script_sha256=sha(a.sensitivity_script),r_script_sha256=sha(a.r_script),environment=dict(prefix=str(prefix),manifest_sha256=item['manifest_sha256'],lock_sha256=sha(lock)),
        whole_genome_absence_producer=clinical.get('absence_producer_implemented',False),scientific_calibration=False,nextflow_runtime_validated=False,
        outputs={f.name:sha(f) for f in out.iterdir() if f.is_file()})
    nested={str(f.relative_to(out)):sha(f) for f in out.glob('sensitivity_*/*') if f.is_file()}
    if nested:
        save(out/'sensitivity_inventory.json',dict(outputs=nested))
        record['outputs']['sensitivity_inventory.json']=sha(out/'sensitivity_inventory.json')
    save(out/'host_source_manifest.json',record)
    # A core run may contain one species: never publish partial-species q-values as global F1.
    aggregate(files['family.tsv'],[out],out,allow_incomplete=True)
    record['outputs'].update({n:sha(out/n) for n in ('F1_results.tsv','F1_manifest.json')})
    # F1 references immutable model table/config hashes, not its own containing manifest.
    save(out/'host_source_manifest.json',record)


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('check');q.add_argument('--config',required=True);q.add_argument('--bundle');q.add_argument('--production',action='store_true')
    q=s.add_parser('run')
    for n in ('config','bundle','research','mask','tree','outdir'):q.add_argument('--'+n,required=True)
    q.add_argument('--recombination-result')
    q.add_argument('--lineages')
    q.add_argument('--sensitivity-script',default=str(ROOT/'bin/host_sensitivity.R'))
    q.add_argument('--r-script',default=str(ROOT/'bin/host_source_models.R'));q.add_argument('--timeout',type=int,default=3600)
    q=s.add_parser('aggregate-f1');q.add_argument('--family',required=True);q.add_argument('--sources',nargs='+',required=True);q.add_argument('--outdir',required=True)
    a=p.parse_args()
    if a.action=='check':print(json.dumps(inspect(a.config,a.bundle,a.production),indent=2))
    elif a.action=='aggregate-f1':aggregate(a.family,a.sources,new_output(a.outdir))
    else:execute(a)


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,RuntimeError,KeyError,subprocess.TimeoutExpired) as e:
        print('ERROR: '+str(e),file=sys.stderr);raise SystemExit(2)
