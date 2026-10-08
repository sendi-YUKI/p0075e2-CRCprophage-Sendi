"""Bounded metadata contracts only. No subprocess, network, sequence or model calls."""
import ast
import csv
import hashlib
import io
import json
import re
from collections import defaultdict
from pathlib import Path
from jsonschema import Draft202012Validator, FormatChecker
from phageflow import validate_id, unique_index

REPO = Path(__file__).resolve().parents[1]
SCHEMA_PATH = REPO / 'assets/aimbeta/v1/contract.schema.json'
LIMIT = 1024 * 1024
TABLE_LIMIT = 10000  # metadata only; distributed synthetic cases stay <=100 rows/table
ENTITY_TABLE = dict(study='studies', participant='participants', isolate='isolates',
                    genome='genomes', cohort='cohorts', sample='samples',
                    library='libraries', run='runs', observation='observations')
KEYS = {v:k+'_id' for k,v in ENTITY_TABLE.items()}
KEYS.update(evidence='evidence_id', eligibility='eligibility_id', decisions='decision_id', features='feature_id', feature_links='link_id')


def strict_json(text):
    def pairs(items):
        d = {}
        for k,v in items:
            if k in d: raise ValueError('Duplicate JSON key: '+k)
            d[k] = v
        return d
    def constant(v): raise ValueError('Nonfinite JSON number: '+v)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=constant)


def read_small(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size > LIMIT:
        raise ValueError('Missing/nonregular/over-1MiB metadata: '+str(path))
    with path.open('rb') as h: data = h.read(LIMIT+1)
    if len(data)>LIMIT: raise ValueError('Metadata grew beyond limit')
    return data.decode('utf-8-sig')


def local_path(raw, base):
    if re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:',raw) or '\\' in raw:
        raise ValueError('Expected Spark path, not URI/Windows path: '+raw)
    p = Path(raw).expanduser()
    return (p if p.is_absolute() else base/p).resolve()


def metadata_path(raw, base):
    p = local_path(raw,base)
    if not any(p.is_relative_to(r) for r in (REPO, Path('/srv/CRC-PHIRE/projects/crc/p0075e2_beta'))):
        raise ValueError('Metadata must be inside Spark repo or existing research audit directory: '+str(p))
    return p


def schema_validator():
    schema = strict_json(read_small(SCHEMA_PATH))
    def safe(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k in ('$ref','$dynamicRef') and not v.startswith('#/'):
                    raise ValueError('Remote schema references forbidden')
                safe(v)
        elif isinstance(x,list):
            for v in x:safe(v)
    safe(schema)
    Draft202012Validator.check_schema(schema)
    return schema, Draft202012Validator(schema,format_checker=FormatChecker())


def validate_schema(validator, value, label):
    errors = sorted(validator.iter_errors(value),key=lambda e:str(list(e.path)))
    if errors:
        raise ValueError(label+': '+'; '.join('/'.join(map(str,e.path))+': '+e.message for e in errors[:8]))


def row_value(raw, spec):
    if raw=='': return None
    branches=spec.get('anyOf',[spec]); types={x.get('type') for x in branches}
    if 'boolean' in types:
        if raw not in ('true','false'): raise ValueError('TSV booleans must be true/false')
        return raw=='true'
    if types & {'number','integer'}:
        v=strict_json(raw)
        if isinstance(v,bool) or not isinstance(v,(int,float)):raise ValueError('Expected numeric cell')
        return v
    return raw


def load_contract(path):
    path=Path(path).resolve(); schema,validator=schema_validator()
    cfg=strict_json(read_small(path));validate_schema(validator,cfg,'contract')
    tables={}; evidence=[];protected=[path,SCHEMA_PATH]
    def record(p,kind):
        data=read_small(p)
        evidence.append(dict(path=str(p),kind=kind,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
        protected.append(p);return data
    for field in ('source_documents','identity_references'):
        for ref in cfg[field]:
            p=metadata_path(ref['path'],path.parent);record(p,field)
            if evidence[-1]['sha256']!=ref['sha256']:raise ValueError('Reference checksum mismatch: '+str(p))
    for name,raw in cfg['tables'].items():
        tables[name]=[]
        if raw is None:continue
        p=metadata_path(raw,path.parent);text=record(p,name)
        rows=csv.reader(io.StringIO(text),delimiter='\t'); header=next(rows,None)
        spec=schema['$defs']['rows'][name];expected=list(spec['properties'])
        if header is None or len(set(header))!=len(header) or set(header)!=set(expected):
            raise ValueError(name+': exact unique contract columns required; no silently ignored fields')
        row_validator=Draft202012Validator(spec,format_checker=FormatChecker())
        for i,cells in enumerate(rows,2):
            if i>TABLE_LIMIT+1:raise ValueError('Metadata row cap exceeded: '+name)
            if len(cells)!=len(header):raise ValueError(name+': ragged row '+str(i))
            row={k:row_value(v,spec['properties'][k]) for k,v in zip(header,cells)}
            validate_schema(row_validator,row,name+':'+str(i));tables[name].append(row)
        evidence[-1].update(rows=len(tables[name]),headers=header)
    for name,keys in (('genomes',('fasta_path',)),('runs',('read1_path','read2_path'))):
        for r in tables[name]:
            for k in keys:
                if r[k]:protected.append(local_path(r[k],path.parent)) # never open sequence inputs
    return cfg,tables,evidence,protected


def relationships(cfg,t):
    errors=[]; warnings=[];idx={}
    def error(code,msg):errors.append(dict(code=code,message=msg))
    def warn(code,msg):warnings.append(dict(code=code,message=msg))
    for name,rows in t.items():
        try:
            for r in rows:validate_id(r[KEYS[name]])
            idx[name]=unique_index(rows,KEYS[name],name)
        except ValueError as exc:error('unique_id',str(exc));idx[name]={r[KEYS[name]]:r for r in rows}
    def fk(name,key,target,optional=True):
        for r in t[name]:
            v=r[key]
            if v is None and optional:continue
            if v not in idx[target]:error('foreign_key',f'{name}.{key}={v}: missing {target}')
    for args in [('participants','identity_evidence_id','evidence'),('isolates','participant_id','participants'),('isolates','study_id','studies'),('genomes','isolate_id','isolates'),('genomes','pair_evidence_id','evidence'),('cohorts','study_id','studies'),('samples','participant_id','participants'),('samples','cohort_id','cohorts'),('libraries','sample_id','samples'),('runs','library_id','libraries'),('isolates','evidence_id','evidence'),('samples','evidence_id','evidence'),('eligibility','evidence_id','evidence'),('observations','callability_evidence_id','evidence')]:fk(*args)
    for name in ('evidence','eligibility','observations'):
        for r in t[name]:
            if r['entity_id'] not in idx[ENTITY_TABLE[r['entity_type']]]:error('foreign_key',f'{name}: unknown {r["entity_type"]} {r["entity_id"]}')
    def proof(eid,entity_type,entity_id,field,value):
        e=idx['evidence'].get(eid,{})
        return (e.get('entity_type'),e.get('entity_id'),e.get('field'),e.get('value'))==(entity_type,entity_id,field,value) and e.get('level') in ('A','B') and e.get('review_state')=='accepted' and bool(e.get('reviewer')) and bool(e.get('review_date'))
    for e in t['evidence']:
        if e['review_state']=='accepted' and (not e['reviewer'] or not e['review_date']):error('evidence_review','Accepted evidence requires reviewer/date')
    assertions=defaultdict(set)
    for e in t['evidence']:
        if e['review_state']=='accepted' and e['level'] in ('A','B'):
            assertions[(e['entity_type'],e['entity_id'],e['field'])].add(e['value'])
    for name,entity in (('isolates','isolate'),('samples','sample')):
        for r in t[name]:
            if len(assertions[(entity,r[KEYS[name]],'disease_status')])>1 and r['disease_status']!='conflict':error('disease_conflict','Conflicting reviewed clinical sources must remain conflict')
            if r['disease_status'] not in ('unknown','conflict') and not proof(r['evidence_id'],entity,r[KEYS[name]],'disease_status',r['disease_status']):error('disease_evidence',f'{name} {r[KEYS[name]]}: disease requires matching reviewed A/B evidence')
            if r['participant_id'] is None:warn('independence_unconfirmed',f'{name} {r[KEYS[name]]}: participant unknown; no substitute ID generated')
    # Same asserted identity/recruitment across held-out roles is a conflict, never automatic deduplication.
    memberships=defaultdict(set)
    for s in t['studies']:
        if s['recruitment_id']:memberships['recruitment:'+s['recruitment_id']].add(s['role'])
        elif s['role'] in ('development','validation'):warn('independence_unconfirmed','Study '+s['study_id']+' recruitment identity missing')
    for c in t['cohorts']:
        s=idx['studies'].get(c['study_id'],{})
        if s and c['role']!=s['role']:error('cohort_role','Cohort/study roles disagree: '+c['cohort_id'])
    for name in ('isolates','samples'):
        for r in t[name]:
            s=idx['studies'].get(r['study_id'],{}) if name=='isolates' else idx['studies'].get(idx['cohorts'].get(r['cohort_id'],{}).get('study_id'),{})
            pid=r['participant_id'];p=idx['participants'].get(pid,{})
            role=s.get('role','unknown')
            if pid:
                memberships['participant:'+pid].add(role)
                if not proof(p.get('identity_evidence_id'),'participant',pid,'identity',pid):warn('independence_unconfirmed','Participant '+pid+' identity evidence not confirmed')
                if p.get('recruitment_id'):memberships['recruitment:'+p['recruitment_id']].add(role)
                if p.get('recruitment_id') and s.get('recruitment_id') and p['recruitment_id']!=s['recruitment_id']:error('recruitment_mismatch','Participant/study recruitment differs: '+pid)
    for key,roles in memberships.items():
        if {'development','validation'}<=roles:error('holdout_overlap',key+' occurs in development and validation')
    accession={};pairs=defaultdict(list);isolate_genomes=defaultdict(list)
    for g in t['genomes']:
        gid=g['genome_id']
        if g['role']!='discovery' and g['master_member']:error('external_denominator',gid+' external/reference/control cannot enter discovery master')
        if g['role']=='discovery' and g['isolate_id'] is None:warn('isolate_unknown',gid+' cultured isolate relationship not established')
        if g['assembly_accession']:
            acc=g['assembly_accession']
            if acc in accession:error('duplicate_assembly',acc+' repeated as '+accession[acc]+' and '+gid)
            accession[acc]=gid
        if g['assembly_pair_id']:
            if not proof(g['pair_evidence_id'],'genome',gid,'official_assembly_pair',g['assembly_pair_id']):error('pair_evidence',gid+' official pairing evidence missing/mismatched')
            pairs[g['assembly_pair_id']].append(g)
        elif g['pair_evidence_id']:error('pair_evidence',gid+' pair evidence without pair ID')
        if g['master_member'] and g['isolate_id']:isolate_genomes[g['isolate_id']].append(gid)
    for pair,gs in pairs.items():
        if len({g['isolate_id'] for g in gs if g['isolate_id']})>1:error('pair_conflict',pair+' paired records name different isolates')
        if sum(g['master_member'] for g in gs)>1:error('pair_duplicate',pair+' multiple official paired assemblies in master')
    for isolate,gs in isolate_genomes.items():
        if len(gs)>1:warn('master_choice_pending',isolate+' multiple assemblies require human version/pair review: '+','.join(gs))
    selected=defaultdict(list);merges=defaultdict(set);read_paths={}
    for s in t['samples']:
        if s['baseline_selection']=='selected' and s['participant_id']:selected[s['participant_id']].append(s['sample_id'])
    for pid,ss in selected.items():
        if len(ss)>1:error('baseline_conflict',pid+' multiple selected baselines; no automatic choice')
    for r in t['runs']:
        lib=idx['libraries'].get(r['library_id'],{})
        if lib and r['layout']!=lib['layout']:error('layout',r['run_id']+' run/library PE/SE mismatch')
        if r['layout']=='SE' and r['read2_path']:error('layout','SE must not declare read2')
        if r['layout']=='PE' and bool(r['read1_path'])!=bool(r['read2_path']):error('layout','PE declarations require both mates or neither in draft')
        for field in ('read1_path','read2_path'):
            if r[field]:
                key=r[field]
                if key in read_paths:error('duplicate_read_path','Read file declared more than once: '+key)
                read_paths[key]=r['run_id']
        if r['merge_group']:merges[r['merge_group']].add((r['library_id'],lib.get('sample_id'),r['layout']))
    for group,units in merges.items():
        if len(units)>1:error('merge_units',group+' would merge different library/sample/layout entities')
    for r in t['eligibility']:
        if r['status'] in ('eligible','ineligible') and (not r['rule_version'] or not r['evidence_id']):error('eligibility_evidence',r['eligibility_id']+' decided eligibility requires rule/evidence')
        if r['status']=='eligible' and r['endpoint'] in ('catalog','function') and r['entity_type']=='genome':
            g=idx['genomes'].get(r['entity_id'],{})
            if g.get('role')!='discovery' or not g.get('master_member'):error('external_denominator','Catalog/function eligibility outside discovery master')
    seen=set()
    for o in t['observations']:
        key=(o['entity_type'],o['entity_id'],o['feature_unit'],o['feature_id'])
        if key in seen:error('duplicate_observation',str(key))
        seen.add(key)
        if o['value'] is not None:
            if o['value_axis']=='carriage_binary':
                if o['value_unit']!='binary' or o['value'] not in (0,1) or o['carriage'] not in ('present','absent_assessable'):error('value_semantics','Carriage values require binary units and assessed state')
            elif o['value_axis']=='measurement':
                if o['value_unit'] in (None,'binary') or o['measurement'] not in ('detected','not_detected_at_this_depth'):error('value_semantics','Measurement value requires measurement state/nonbinary unit')
                if o['value_unit'] in ('breadth','relative_abundance') and o['value']>1:error('value_semantics','Fraction-valued measurement exceeds one')
                if o['value_unit'] in ('fragments','copies') and not float(o['value']).is_integer():error('value_semantics','Count unit requires an integer')
            else:error('value_semantics','Numeric value requires explicit value axis')
        if o['run_status']!='success' and (o['value'] is not None or o['carriage'] not in (None,'unknown') or o['measurement'] not in (None,'unassessable')):error('state_conflict','Failed/not assessed cannot carry measured result: '+o['observation_id'])
        if o['carriage']=='unknown' and o['value'] is not None and o['value_axis']=='carriage_binary':error('unknown_not_zero','Unknown carriage cannot have numeric value')
        if o['measurement']=='unassessable' and o['value'] is not None and o['value_axis']=='measurement':error('unknown_not_zero','Unassessable measurement cannot have numeric value')
        if o['carriage']=='absent_assessable':
            if o['callability']!='assessable' or not o['rule_version'] or not proof(o['callability_evidence_id'],'observation',o['observation_id'],'callability','assessable'):error('absence_evidence','Absence requires explicit reviewed callability evidence and rule')
            if o['value_axis']=='carriage_binary' and o['value'] not in (None,0):error('state_conflict','Absent carriage cannot have positive value')
        if o['value_axis']=='measurement' and o['measurement']=='detected' and o['value'] is not None and o['value']<=0:error('state_conflict','Detected cannot have zero value')
        if o['value_axis']=='measurement' and o['measurement']=='not_detected_at_this_depth' and o['value'] not in (None,0):error('state_conflict','Nondetection cannot have positive value')
        if o['value_axis']=='carriage_binary' and o['carriage']=='present' and o['value']==0:error('state_conflict','Present carriage cannot have zero value')
    for key,p in cfg['parameters'].items():
        d=idx['decisions'].get(p['decision_id']) if p['decision_id'] else None
        if p['decision_id'] and not d:error('parameter_decision','Missing decision for '+key)
        if p['state']=='frozen' and (p['value'] is None or d is None):error('parameter_decision','Frozen value requires value and decision: '+key)
        if d:
            try:value=strict_json(d['value_json'])
            except ValueError:error('parameter_decision','Invalid value_json: '+key);continue
            if (d['parameter_key'],d['state'],value,d['unit'])!=(key,p['state'],p['value'],p['unit']):error('parameter_decision','Decision contradicts parameter '+key)
    for key,r in cfg['rules'].items():
        if r['state']=='frozen' and (not r['version'] or not r['decision_ref']):error('rule_decision','Frozen rule requires version and decision reference: '+key)
    for a,b in (('high_completeness_min','completeness_min'),('contamination_max','high_contamination_max')):
        av=cfg['parameters'][a]['value'];bv=cfg['parameters'][b]['value']
        if av is not None and bv is not None and av<bv:error('parameter_conflict',a+' cannot be below '+b)
    return errors,warnings


def feature_relationships(t):
    errors=[];idx={r['feature_id']:r for r in t['features']};genomes={r['genome_id'] for r in t['genomes']}
    def err(msg):errors.append(dict(code='feature_relationship',message=msg))
    for f in t['features']:
        if f['genome_id'] and f['genome_id'] not in genomes:err('Unknown source genome: '+f['feature_id'])
        if f['candidate_id']:
            c=idx.get(f['candidate_id'],{})
            if c.get('feature_unit')!='candidate' or c.get('genome_id')!=f['genome_id']:err('Candidate/source genome mismatch: '+f['feature_id'])
        if (f['start'] is None)!=(f['end'] is None):err('Incomplete coordinates: '+f['feature_id'])
        if f['start'] is not None:
            if f['coordinate_system']=='one_based_closed' and f['start']<1:err('One-based start must be positive')
            if f['coordinate_system']=='zero_based_half_open' and f['end']<=f['start']:err('Half-open interval must have positive length')
            if f['end']<f['start'] or not f['genome_id'] or not f['contig_id'] or not f['strand']:err('Coordinates require ordered bounds, source genome, contig and strand: '+f['feature_id'])
        if bool(f['source_artifact'])!=bool(f['source_row_id']):err('Source artifact and row identity must be paired')
    kinds={'votu_member':('votu','candidate'),'candidate_gene':('candidate','gene'),'system_gene':('system','gene'),'quant_group_member':('quantification_group','votu'),'locus_candidate':('locus','candidate')}
    seen=set()
    for r in t['feature_links']:
        key=(r['parent_feature_id'],r['member_feature_id'],r['relation'])
        if key in seen:err('Duplicate feature membership: '+str(key))
        seen.add(key);a=idx.get(r['parent_feature_id'],{});b=idx.get(r['member_feature_id'],{})
        if (a.get('feature_unit'),b.get('feature_unit'))!=kinds[r['relation']]:err('Missing feature or incompatible link units: '+r['link_id'])
        if r['relation'] in ('candidate_gene','system_gene','locus_candidate') and a.get('genome_id')!=b.get('genome_id'):err('Physical membership spans different source genomes')
    for o in t['observations']:
        if idx.get(o['feature_id'],{}).get('feature_unit')!=o['feature_unit']:err('Observation feature missing/unit mismatch: '+o['observation_id'])
    return errors


def parameter_plan(cfg):
    """Extract literal declarations via AST/regex; never execute configs or launchers."""
    tree=ast.parse(read_small(REPO/'bin/run_core.py'));defaults={}
    for node in ast.walk(tree):
        if isinstance(node,ast.For) and isinstance(node.target,ast.Tuple):
            if [getattr(x,'id',None) for x in node.target.elts]==['option','default']:
                try:defaults.update({k.replace('-','_'):v for k,v in ast.literal_eval(node.iter)})
                except (ValueError,TypeError):pass
    conf=read_small(REPO/'conf/core.config')
    core={k:strict_json(v) for k,v in re.findall(r'^\s*(\w+)\s*=\s*([0-9.]+)\s*$',conf,re.M)}
    profile=strict_json(read_small(REPO/'assets/profiles/cohort_measurement_dev_v1.json'))
    # 2026-10-01: dev profiles carry research_eligible=false. Crosswalking against them
    # is allowed (they are the current engineering values) but a study contract must not
    # inherit them as frozen research thresholds. See docs/ELEMENT_CLASS_SPEC_20261001.md
    # and assets/profiles/cohort_measurement_research_v1.json.
    if profile.get('research_eligible') is False:
        profile=dict(profile,_research_eligible_note='engineering profile; not valid for a frozen study contract')
    planned={};unresolved=[]
    for key,p in cfg['parameters'].items():
        if key in ('identity_min','aligned_query_fraction_min','mapq_min'):
            current=profile.get(key);source='assets/profiles/cohort_measurement_dev_v1.json';agrees=current is not None
        else:
            current=defaults.get(key);source='bin/run_core.py argparse -> explicit Nextflow CLI; conf/core.config literal';agrees=current is not None and key in core and core[key]==current
        if not agrees:unresolved.append(key+': literals missing or disagree; no effective-value claim')
        planned[key]=dict(current_declared_default=current,source=source,literals_agree=agrees,research_value=p['value'],research_state=p['state'],unit=p['unit'],proposed_preview=p['value'],activated=False,check_value=p['value'],check_purpose='type/unit/consistency only; no analysis',effective_runtime_value=None)
    nf=read_small(REPO/'nextflow.config');spark=read_small(REPO/'conf/spark_native.config')
    def m(pattern,text):
        found=re.search(pattern,text);return found.group(1) if found else None
    return dict(status='declared_static_not_runtime_resolved',parameters=planned,unresolved=unresolved,overrides_supplied=[],precedence=['launcher argparse default or future explicit CLI','launcher passes supported numeric options as explicit Nextflow parameters','Spark native config loaded last for resource/runtime overrides','Nextflow CLI parameters have priority; invocation and config execution intentionally not performed'],pipeline_version=m(r"version\s*=\s*'([^']+)'",nf),nextflow_engine_constraint=m(r"nextflowVersion\s*=\s*'([^']+)'",nf),nextflow_engine_actual='not_reprobed_this_phase',spark_config_contains_processbuilder='ProcessBuilder' in spark,resources=dict(slurm_declared='6 CPU / 48 GiB',nextflow_memory_literal=m(r"memory\s*=\s*'([^']+)'",spark),task_resource_limit='46.GB',cohort_budget_gb=46,queue_size=1,unit_status='Nextflow v25.10.4 MemoryUnit source uses 1024**3 for GB, so 46 GB means 46 GiB; runtime engine/allocation not verified here'),runtime_unresolved=['future explicit CLI','profile/config evaluation and environment interpolation','Nextflow engine and actual Slurm allocation','database payload readiness not rechecked'])


def evaluate(cfg,t):
    errors,warnings=relationships(cfg,t);errors.extend(feature_relationships(t))
    plan=parameter_plan(cfg);reasons=[]
    crosswalk=strict_json(read_small(REPO/'assets/aimbeta/v1/field_crosswalk.json'))
    if not cfg['source_documents']:reasons.append('source_documents_missing')
    if not cfg['identity_references']:reasons.append('code_model_identity_references_missing')
    if cfg['master_universe_version'] is None:reasons.append('master_universe_not_versioned')
    if not any(g['master_member'] for g in t['genomes']):reasons.append('discovery_master_empty')
    required={'core':['genomes','isolates'],'within_species':['genomes','isolates','participants','eligibility'],'cohort':['cohorts','samples','libraries','runs'],'validation':['cohorts','samples','participants'],'experiments':[],'preparation':[]}
    for name in required[cfg['stage']]:
        if not t[name]:reasons.append('stage_table_empty:'+name)
    if cfg['status']!='frozen':reasons.append('contract_not_frozen')
    for k,p in cfg['parameters'].items():
        if p['state']!='frozen' or p['value'] is None:reasons.append('parameter_not_frozen:'+k)
    for k,r in cfg['rules'].items():
        if r['state']!='frozen':reasons.append('rule_not_frozen:'+k)
    if warnings:reasons.extend(sorted({w['code'] for w in warnings}))
    if cfg['stage'] in ('within_species','cohort','validation'):
        reasons.extend('statistic_unspecified:'+k for k,v in cfg['statistics'].items() if v is None)
    if cfg['stage']=='experiments' and any(v is None for v in cfg['experiments'].values()):reasons.append('experiment_endpoint_or_lock_missing')
    for g in t['genomes']:
        if g['master_member'] and (not g['fasta_path'] or not g['fasta_sha256']):reasons.append('genome_input_identity_incomplete:'+g['genome_id'])
    for r in t['runs']:
        if not r['read1_path']:reasons.append('run_input_declaration_incomplete:'+r['run_id'])
    unsupported=[]
    s=cfg['statistics']
    if s['outcome']=='CRC' or s['exposure']=='feature':unsupported.append('CRC-source logistic/MDS opt-in consumer implemented; separate frozen host-source config/bundle and scientific validation required; old B2 remains feature~disease')
    if s['lineage_method']=='kinship':unsupported.append('Opt-in allelic kinship/rrBLUP LMM implemented; its configuration and scientific calibration are separate from lineage')
    if s['family']=='global_F1_F2':unsupported.append('Opt-in global F1 and separate community F2/F3 implemented; complete family and model contracts required; legacy B2 unchanged')
    if s['correction']=='Holm':unsupported.append('Independent holdout/Holm implemented in community_analysis; candidate/model/family freeze and independent recruitment required')
    if cfg['stage']=='within_species':unsupported.extend(['opt-in masked core and conditional reference-coordinate recombination/kinship implemented; applicability freeze and scientific validation pending','opt-in anchor/fragment and locus evidence implemented; scientific freeze/calibration and whole-genome absence pending'])
    if cfg['stage'] in ('cohort','validation'):unsupported.extend(['conserved nonmobile host core','quantification group measurement and calibrated detection limits'])
    if cfg['real_execution_authorized']:warnings.append(dict(code='authorization_ignored',message='File true is audit information only; this checker cannot execute anything'))
    return dict(schema_version='aimbeta-report-1.0',mode='metadata_check_only',schema_valid=True,relationships_consistent=not errors,contract_valid=not errors,errors=errors,warnings=warnings,declared_stage=cfg['stage'],stage_input_declarations_complete=not any('empty' in x or 'missing' in x or 'incomplete' in x or 'unspecified' in x for x in reasons),stage_readiness='not_ready',readiness_reasons=sorted(set(reasons+plan['unresolved']+['biological_payloads_not_inspected','scientific_calibration_and_stage_acceptance_not_performed'])),field_mapping=dict(version=crosswalk['version'],fields=crosswalk['fields'],status='explicit_contract_crosswalk_not_full_analysis_wiring'),analysis_consumers=dict(existing_core_mapping='partial_draft_only',new_research_methods='opt_in_host_core_KOfam_catalog_CRC_source_MDS_F1_kinship_Firth_and_community_F2_F3_holdout_implemented; scientific_readiness_separate',unsupported=unsupported),evidence_truth_verified=False,independent_people_count_verified=False,entity_counts={k:len({r[KEYS[k]] for r in rows}) for k,rows in t.items()},independent_unit_note='Counts are distinct declared IDs, not verified independent people; runs and assemblies never inflate people counts.',parameter_plan=plan,execution_authorized=False,nextflow_started=False,tasks_submitted=0,tools_executed=0,network_requests=0,real_data_validated=False,scientifically_calibrated=False)
