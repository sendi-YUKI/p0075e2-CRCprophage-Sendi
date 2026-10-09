#!/usr/bin/env python3
"""V1 occurrence/sequence catalog adapter; no bacterial-genome QC inference."""
import argparse
from collections import defaultdict
import json
import os
import importlib.util
from pathlib import Path
import shutil
from types import SimpleNamespace

import phageflow as pf
import core_votu as cv

SCHEMA = 'viral_catalog_v1'
SCOPES = {'host_resolved_prophage', 'virome_discovery', 'combined_exploratory'}
DISCOVERY_SOURCE_KINDS = {'viral_assembly', 'metagenome_assembly', 'virome_assembly', 'independent_fasta'}
VIRUS_SCOPES = {'bacteriophage', 'archaeal_virus', 'eukaryotic_virus', 'unknown'}
POLICY = dict(version='analysis1_v4_primary_length5000_medium_quality_v1', min_length_bp=5000,
              min_checkv_completeness_percent=50, max_checkv_contamination_percent=10,
              identity='retained_geNomad_candidate', missing_quality='unknown',
              scientifically_calibrated=False, source_slice='verified_this_run_required', host_scope='requires_host_QC_eligible',
              boundary='provirus_locus_or_full_contig_extent_required')
OCC_FIELDS = ['source_occurrence_id', 'viral_sequence_id', 'candidate_id', 'legacy_candidate_id',
              'assembly_id', 'biological_sample_id', 'library_id', 'material_type', 'source_kind',
              'source_host_genome', 'predicted_host', 'sample_host_evidence', 'contig_id',
              'original_contig_id', 'start0', 'end0', 'strand', 'caller_strand', 'candidate_type',
              'boundary_version', 'sequence_sha256', 'slice_sha256', 'source_fasta_sha256',
              'slice_verification', 'integration_evidence', 'boundary_confidence',
              'source_genome_eligibility', 'host_resolved_prophage', 'virus_scope',
              'virus_scope_evidence', 'virus_scope_evidence_version', 'checkv_status',
              'completeness', 'contamination', 'checkv_quality', 'length', 'caller', 'viral_identity_confidence', 'virus_score', 'fdr', 'n_hallmarks', 'taxonomy']
SEQ_FIELDS = ['viral_sequence_id', 'candidate_id', 'length', 'sequence_sha256', 'virus_scope', 'viral_identity_confidence',
              'virus_scope_evidence', 'virus_scope_evidence_version', 'n_source_occurrences',
              'checkv_quality', 'completeness', 'contamination', 'viral_catalog_eligibility',
              'eligibility_reason', 'eligibility_policy', 'catalog_scope', 'boundary_confidence']
MEMBER_FIELDS = ['catalog_version', 'catalog_scope', 'votu_id', 'viral_sequence_id', 'candidate_id',
                 'representative_id', 'membership_status', 'matching_representative_ids']


_IPC_REJECT = None


def _native_input(path):
    if os.environ.get('CRC_PHAGE_RUNTIME') == 'spark-native':
        global _IPC_REJECT
        if _IPC_REJECT is None:
            repo = os.environ.get('CRC_SPARK_REPO')
            helper = Path(repo) / 'bin/spark_fs_guard.py' if repo else Path(__file__).resolve().with_name('spark_fs_guard.py')
            spec = importlib.util.spec_from_file_location('_crc_viral_fs_guard', helper)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            _IPC_REJECT = module.reject_ipc_storage
        _IPC_REJECT([path])


def table(path, fields=()):
    _native_input(path)
    return pf.read_tsv(Path(path), fields)[1]


def write(path, rows, fields):
    pf.write_tsv(Path(path), fields, rows)


def stable(prefix, payload):
    return prefix + pf.sha256_bytes(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode())[:32]


def guard(out, inputs):
    out = Path(out).resolve()
    for path in inputs:
        path = Path(path).resolve()
        if out == path or out.is_relative_to(path) or path.is_relative_to(out):
            raise ValueError(f'Output overlaps protected input/database: {path}')
    out.mkdir(parents=True, exist_ok=True)
    return out


def load_units(path, verify=True):
    source = Path(path).resolve()
    _native_input(source)
    doc = json.loads(source.read_text())
    if doc.get('schema_version') not in {1, SCHEMA}:
        raise ValueError('Unsupported/missing assembly manifest schema_version')
    if not isinstance(doc.get('units'), list):
        raise ValueError('Manifest requires units array; an empty array explicitly means zero assemblies')
    units, seen = [], set()
    for raw in doc['units']:
        row = dict(raw)
        uid = pf.validate_id(row.get('assembly_id', row.get('unit_id', '')))
        if uid in seen:
            raise ValueError('Duplicate assembly ID: ' + uid)
        seen.add(uid)
        p = Path(row['assembly_fasta']).expanduser()
        p = (p if p.is_absolute() else source.parent / p).resolve()
        _native_input(p)
        if not isinstance(row.get('sha256'), str) or len(row['sha256']) != 64 or any(c not in '0123456789abcdef' for c in row['sha256']):
            raise ValueError('Assembly SHA256 must be lowercase 64-character hex')
        if verify and (not p.is_file() or pf.file_hash(p) != row.get('sha256')):
            raise ValueError('Missing assembly or SHA256 mismatch: ' + str(p))
        kind = row.get('input_kind', 'viral_assembly')
        if kind not in {'viral_assembly', 'metagenome_assembly', 'virome_assembly', 'independent_fasta', 'bacterial_genome'}:
            raise ValueError('Unsupported manifest input_kind: ' + kind)
        scope = row.get('virus_scope', 'unknown')
        if scope not in VIRUS_SCOPES:
            raise ValueError('Unsupported virus_scope')
        if scope != 'unknown' and not all(row.get(k) for k in ['virus_scope_evidence', 'virus_scope_evidence_version']):
            raise ValueError('Known virus scope needs explicit evidence and version')
        if row.get('contig_map'):
            mp = Path(row['contig_map'])
            mp = (mp if mp.is_absolute() else source.parent / mp).resolve()
            _native_input(mp)
            if verify and not mp.is_file():
                raise ValueError('Missing assembly contig map: ' + str(mp))
            row['contig_map'] = str(mp)
            if verify:
                row['contig_map_sha256'] = pf.file_hash(mp)
        row.update(assembly_id=uid, unit_id=uid, assembly_fasta=str(p), input_kind=kind,
                   biological_sample_id=row.get('biological_sample_id') or row.get('sample_id') or 'unknown',
                   library_id=row.get('library_id') or 'unknown', material_type=row.get('material_type') or 'unknown',
                   virus_scope=scope, virus_scope_evidence=row.get('virus_scope_evidence', 'not_assessed'),
                   virus_scope_evidence_version=row.get('virus_scope_evidence_version', 'not_assessed'))
        units.append(row)
    return sorted(units, key=lambda r: r['assembly_id'])


def validate_manifest(args):
    units = load_units(args.manifest)
    out = guard(args.outdir, [args.manifest] + [u['assembly_fasta'] for u in units])
    pf.write_json(out / 'assembly_manifest.json', dict(schema_version=SCHEMA, units=units,
                  source_statuses=json.loads(Path(args.manifest).read_text()).get('source_statuses', []),
                  input_manifest_sha256=pf.file_hash(Path(args.manifest))))
    write(out / 'assemblies.tsv', units, ['assembly_id', 'assembly_fasta', 'sha256'])


def prepare(args):
    _native_input(args.fasta)
    # Existing A0 routines supply tested FASTA normalization only. Their internal
    # genome_id is a computational namespace, never host or CheckM2 evidence.
    if args.expected_sha256 and pf.file_hash(Path(args.fasta)) != args.expected_sha256:
        raise ValueError('Assembly SHA256 mismatch before geNomad')
    pf.prepare(SimpleNamespace(genome_id=args.assembly_id, fasta=args.fasta, outdir=args.outdir))
    out = Path(args.outdir)
    rows = table(out / 'contig_map.tsv')
    write(out / 'assembly_contig_map.tsv', [dict(r, assembly_id=args.assembly_id) for r in rows],
          ['assembly_id', 'contig_id', 'original_contig_id', 'original_header', 'length', 'sequence_sha256'])
    pf.write_json(out / 'namespace_semantics.json', dict(internal_genome_id='assembly namespace only',
                  source_host_genome='unknown', source_genome_eligibility='not_applicable',
                  bacterial_genome_qc_performed=False))


def number(value):
    try:
        n = float(value)
        return n if 0 <= n <= 100 else None
    except (ValueError, TypeError):
        return None


def eligibility(row, policy=POLICY):
    if int(row['length']) < policy['min_length_bp']:
        return 'ineligible', 'below_primary_length_5000bp_retained_in_raw_catalog'
    if row['virus_scope'] == 'eukaryotic_virus':
        return 'unknown', 'CheckV_prokaryotic_model_applicability_unsupported'
    q, c = number(row['completeness']), number(row['contamination'])
    if q is None or c is None:
        return 'unknown', 'CheckV_estimate_missing_or_inconsistent'
    if q < policy['min_checkv_completeness_percent'] or c > policy['max_checkv_contamination_percent']:
        return 'ineligible', 'below_development_quality_or_contamination_threshold'
    return 'eligible', 'development_policy_only' if row['virus_scope'] != 'unknown' else 'provisional_prokaryotic_virus_assumption'


def applicable(scope, tool):
    if scope == 'eukaryotic_virus':
        return 'unsupported' if tool in {'CheckV', 'DefenseFinder', 'AntiDefenseFinder', 'PHROGs'} else 'provisional'
    if scope == 'unknown':
        return 'provisional'
    if scope == 'archaeal_virus' and tool in {'PHROGs', 'DefenseFinder', 'AntiDefenseFinder'}:
        return 'provisional'
    return 'applicable'


def make_occurrence(hit, sequence, unit, source_records=None, legacy=False):
    digest = pf.sequence_hash(sequence)
    if digest != hit['sequence_sha256'] or len(sequence) != int(hit['length']):
        raise ValueError('Candidate sequence/hash mismatch: ' + hit['candidate_id'])
    start, end = int(hit['start0']), int(hit['end0'])
    if not 0 <= start < end or end - start != len(sequence):
        raise ValueError('Invalid/wrapping source coordinates')
    verified = 'legacy_record_preserved_not_revalidated'
    if source_records is not None:
        if hit['contig_id'] not in source_records:
            raise ValueError('Missing source contig: ' + hit['contig_id'])
        source = source_records[hit['contig_id']][1]
        sliced = source[start:end]
        if end > len(source) or hit['strand'] not in {'+', '-'}:
            raise ValueError('Out-of-range source coordinates or unsupported orientation')
        if hit['strand'] == '-':
            sliced = pf.reverse_complement(sliced)
        if sliced != sequence:
            raise ValueError('Source slice does not reproduce viral sequence')
        verified = 'verified_this_run'
    host = hit.get('source_host_genome', '') if legacy or unit['input_kind'] == 'bacterial_genome' else ''
    host = host or 'unknown'
    resolved = (host.lower() not in {'unknown', 'not_assessed', 'not_applicable', ''} and hit['candidate_type'] == 'provirus_locus'
                and hit.get('integration_evidence', 'unknown').lower() not in {'', 'unknown', 'not_assessed', 'not_applicable', 'none', 'absent', 'false'} and verified == 'verified_this_run')
    oid = stable('occ_', [unit['assembly_id'], hit['contig_id'], start, end, hit['strand'], hit['boundary_version'], digest])
    return dict(source_occurrence_id=oid, viral_sequence_id='vs_' + digest,
                candidate_id=hit['candidate_id'], legacy_candidate_id=hit['candidate_id'] if legacy else '',
                assembly_id=unit['assembly_id'], biological_sample_id=unit['biological_sample_id'],
                library_id=unit['library_id'], material_type=unit['material_type'],
                source_kind='legacy_bacterial_genome' if legacy else unit['input_kind'],
                source_host_genome=host, predicted_host='not_assessed', sample_host_evidence='not_assessed',
                source_fasta_sha256=unit['sha256'], slice_verification=verified, slice_sha256=digest,
                source_genome_eligibility=hit.get('main_catalog_eligible', 'unknown') if legacy else 'not_assessed' if unit['input_kind'] == 'bacterial_genome' else 'not_applicable',
                host_resolved_prophage=str(resolved).lower(), virus_scope=unit['virus_scope'],
                virus_scope_evidence=unit['virus_scope_evidence'], virus_scope_evidence_version=unit['virus_scope_evidence_version'],
                **{k: hit.get(k, 'unknown') for k in ['contig_id', 'original_contig_id', 'start0', 'end0',
                   'strand', 'caller_strand', 'candidate_type', 'boundary_version', 'sequence_sha256',
                   'integration_evidence', 'boundary_confidence', 'checkv_status', 'completeness',
                   'contamination', 'checkv_quality', 'length', 'caller', 'viral_identity_confidence', 'virus_score', 'fdr', 'n_hallmarks', 'taxonomy']})


def sequence_identity_confidence(occurrences):
    values = {row.get('viral_identity_confidence') or 'unknown' for row in occurrences}
    return next(iter(values)) if len(values) == 1 else 'mixed_source_evidence_see_occurrences' if values else 'unknown'


def reject_research_legacy(path):
    if not path:
        return
    p = Path(path) / 'catalog_manifest.json'
    if p.is_file():
        metadata = json.loads(p.read_text())
        if str(metadata.get('catalog_version', '')).startswith('research_'):
            raise ValueError('Research anchor/supplemental catalog cannot be silently reclustered by the legacy V1 importer; dedicated versioned crosswalk is pending')


def assemble(args):
    reject_research_legacy(args.legacy_catalog)
    units = load_units(args.manifest, verify=False)
    by_unit = {u['assembly_id']: u for u in units}
    protected = list(args.sample_dirs) + list(args.prepared_dirs) + [args.manifest]
    if args.legacy_catalog:
        protected.append(args.legacy_catalog)
    out = guard(args.outdir, protected)
    occurrences, seqs, statuses, seen = [], {}, [], set()
    prep = {}
    for directory in args.prepared_dirs:
        d = Path(directory)
        doc = json.loads((d / 'input_manifest.json').read_text())
        if doc['genome_id'] in prep:
            raise ValueError('Duplicate prepared assembly')
        prep[doc['genome_id']] = (d, pf.read_fasta(d / 'genome.fna'))
    for directory in args.sample_dirs:
        d = Path(directory)
        state = table(d / 'sample_status.tsv')
        if len(state) != 1 or state[0]['sample_status'] not in {'zero_candidates', 'completed_with_candidates'}:
            raise ValueError('Failed/incomplete assembly must not become empty viral catalog')
        aid = state[0]['genome_id']
        if aid in seen or aid not in by_unit or aid not in prep:
            raise ValueError('Assembly status foreign-key or duplication error')
        seen.add(aid)
        hits = table(d / 'prophage_master.tsv', pf.MASTER_FIELDS)
        records = pf.read_fasta(d / 'candidates.fna', allow_empty=True)
        if {h['candidate_id'] for h in hits} != set(records) or int(state[0]['candidate_count']) != len(hits):
            raise ValueError('Candidate table/FASTA/status disagreement')
        for hit in hits:
            seq = records[hit['candidate_id']][1]
            occ = make_occurrence(hit, seq, by_unit[aid], prep[aid][1])
            occurrences.append(occ)
            seqs[occ['viral_sequence_id']] = seq
        statuses.append(dict(assembly_id=aid, status=state[0]['sample_status'], candidate_count=len(hits),
                             checkv_status=state[0]['checkv_status'], source_genome_eligibility='not_applicable',
                             biological_sample_id=by_unit[aid]['biological_sample_id'], library_id=by_unit[aid]['library_id'],
                             upstream_status='success', genomad_status='completed', reason=''))
    if seen != set(by_unit):
        raise ValueError('Missing final assembly results; never treat a failed/missing unit as zero')
    for state in json.loads(Path(args.manifest).read_text()).get('source_statuses', []):
        if state.get('material_type') not in {'bulk_metagenome', 'vlp'}:
            continue
        aid = state['unit_id']
        if aid in seen:
            if state['status'] != 'success':
                raise ValueError('Upstream assembly status contradicts supplied successful unit')
            continue
        if state['status'] not in {'not_assessed', 'zero_contigs'}:
            raise ValueError('Failed/missing successful upstream assembly cannot become a zero viral result')
        statuses.append(dict(assembly_id=aid, status='not_assessed', candidate_count=0, checkv_status='skipped_upstream_not_assessed',
                             source_genome_eligibility='not_applicable', biological_sample_id=state.get('biological_sample_id') or 'unknown',
                             library_id=state.get('library_id') or 'unknown', upstream_status=state['status'], genomad_status='not_assessed', reason=state.get('reason', state['status'])))
    legacy_files = []
    if args.legacy_catalog:
        legacy = Path(args.legacy_catalog)
        _native_input(legacy / 'candidates.fna')
        records = pf.read_fasta(legacy / 'candidates.fna', allow_empty=True)
        hits = table(legacy / 'prophage_master.tsv', pf.MASTER_FIELDS)
        if set(records) != {h['candidate_id'] for h in hits}:
            raise ValueError('Legacy candidate FASTA/master foreign-key mismatch')
        staged_sources = getattr(args, 'legacy_source_fastas', None)
        source_list = load_units(args.legacy_source_manifest, verify=staged_sources is None) if args.legacy_source_manifest else []
        if staged_sources is not None:
            if len(staged_sources) != len(source_list):
                raise ValueError('Legacy source FASTA count mismatch')
            for u, p in zip(source_list, staged_sources):
                _native_input(p)
                if pf.file_hash(Path(p)) != u['sha256']:
                    raise ValueError('Legacy source checksum mismatch')
                u['assembly_fasta'] = str(Path(p).resolve())
        source_units = {u['assembly_id']: u for u in source_list}
        for hit in hits:
            aid = hit['genome_id']
            unit = source_units.get(aid, dict(assembly_id=aid, biological_sample_id='unknown', library_id='unknown',
                    material_type='bacterial_genome', input_kind='bacterial_genome', sha256='not_provided',
                    virus_scope='unknown', virus_scope_evidence='not_assessed', virus_scope_evidence_version='not_assessed'))
            src = pf.read_fasta(Path(unit['assembly_fasta'])) if aid in source_units else None
            occ = make_occurrence(hit, records[hit['candidate_id']][1], unit, src, legacy=True)
            occurrences.append(occ)
            seqs[occ['viral_sequence_id']] = records[hit['candidate_id']][1]
        # Preserve every original table, along with sequence/JSON evidence, byte for byte.
        for p in sorted(legacy.rglob('*')):
            if p.is_file() and p.suffix in {'.tsv', '.csv', '.json', '.fna', '.faa', '.gff', '.ffn'}:
                _native_input(p)
                rel = p.relative_to(legacy)
                dest = out / 'legacy' / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(p, dest)
                legacy_files.append(dict(path=str(rel), sha256=pf.file_hash(p), bytes=p.stat().st_size))
    if len({r['source_occurrence_id'] for r in occurrences}) != len(occurrences):
        raise ValueError('Duplicate source occurrences; import each source once')
    by_seq = defaultdict(list)
    for row in occurrences:
        by_seq[row['viral_sequence_id']].append(row)
    masters = []
    for sid, seq in sorted(seqs.items()):
        occs = by_seq[sid]
        known_scope = {r['virus_scope'] for r in occs} - {'unknown'}
        scope = next(iter(known_scope)) if len(known_scope) == 1 else 'unknown'
        estimates = [number(r['completeness']) for r in occs]
        contamination = [number(r['contamination']) for r in occs]
        row = dict(viral_sequence_id=sid, candidate_id=sid, length=len(seq), sequence_sha256=pf.sequence_hash(seq),
                   virus_scope=scope, viral_identity_confidence=sequence_identity_confidence(occs), virus_scope_evidence=';'.join(sorted({r['virus_scope_evidence'] for r in occs})),
                   virus_scope_evidence_version=';'.join(sorted({r['virus_scope_evidence_version'] for r in occs})),
                   n_source_occurrences=len(occs), checkv_quality=sorted(occs, key=lambda r: cv.rank_candidate(dict(r, candidate_id=sid)))[-1]['checkv_quality'],
                   completeness=min(estimates) if all(v is not None for v in estimates) else 'unknown',
                   contamination=max(contamination) if all(v is not None for v in contamination) else 'unknown',
                   eligibility_policy=POLICY['version'], catalog_scope=args.catalog_scope,
                   boundary_confidence='multiple_occurrences_preserved')
        row['viral_catalog_eligibility'], row['eligibility_reason'] = eligibility(row)
        if not any(o['caller'] == 'geNomad' for o in occs):
            row.update(viral_catalog_eligibility='unknown', eligibility_reason='required_geNomad_identity_evidence_missing')
        if not any(o['candidate_type'] in {'provirus_locus', 'full_contig_virus'} for o in occs):
            row.update(viral_catalog_eligibility='unknown', eligibility_reason='candidate_extent_uncertain')
        if not any(o['slice_verification'] == 'verified_this_run' for o in occs):
            row.update(viral_catalog_eligibility='unknown', eligibility_reason='source_slice_not_revalidated')
        in_scope = (args.catalog_scope == 'combined_exploratory' or
                    (args.catalog_scope == 'virome_discovery' and any(o['source_kind'] in DISCOVERY_SOURCE_KINDS for o in occs)) or
                    (args.catalog_scope == 'host_resolved_prophage' and any(o['host_resolved_prophage'] == 'true' for o in occs)))
        if not in_scope:
            row.update(viral_catalog_eligibility='ineligible', eligibility_reason='outside_explicit_catalog_scope')
        if args.catalog_scope == 'host_resolved_prophage' and in_scope:
            qualifying = [o for o in occs if o['host_resolved_prophage'] == 'true']
            if not any(o['source_genome_eligibility'].lower() in {'true','eligible','pass','1','yes'} for o in qualifying):
                unknown = any(o['source_genome_eligibility'].lower() in {'unknown','not_assessed','not_applicable',''} for o in qualifying)
                row.update(viral_catalog_eligibility='unknown' if unknown else 'ineligible', eligibility_reason='source_genome_QC_not_eligible_for_host_catalog')
        if len(known_scope) > 1:
            row.update(viral_catalog_eligibility='unknown', eligibility_reason='conflicting_virus_scope_evidence')
        masters.append(row)
    write(out / 'source_occurrences.tsv', sorted(occurrences, key=lambda r: r['source_occurrence_id']), OCC_FIELDS)
    write(out / 'host_resolved_prophage.tsv', [r for r in occurrences if r['host_resolved_prophage'] == 'true'], OCC_FIELDS)
    write(out / 'viral_sequence_master.tsv', masters, SEQ_FIELDS)
    write(out / 'candidate_sequence_crosswalk.tsv', occurrences,
          ['candidate_id', 'legacy_candidate_id', 'source_occurrence_id', 'viral_sequence_id', 'boundary_version', 'sequence_sha256'])
    write(out / 'assembly_status.tsv', statuses, ['assembly_id', 'biological_sample_id', 'library_id', 'upstream_status', 'status', 'candidate_count', 'genomad_status', 'checkv_status', 'source_genome_eligibility', 'reason'])
    applicability = [dict(viral_sequence_id=r['viral_sequence_id'], virus_scope=r['virus_scope'], tool=t,
                         applicability=applicable(r['virus_scope'], t), evidence='scope_declared_evidence_or_explicit_provisional_unknown',
                         policy='viral_tool_applicability_v1', calibration='not_performed')
                     for r in masters for t in ['geNomad', 'CheckV', 'Pyrodigal-gv', 'PHROGs', 'DefenseFinder', 'AntiDefenseFinder']]
    write(out / 'tool_applicability.tsv', applicability, ['viral_sequence_id', 'virus_scope', 'tool', 'applicability', 'evidence', 'policy', 'calibration'])
    pf.write_fasta(out / 'viral_sequences.fna', sorted(seqs.items()))
    pf.write_json(out / 'discovery_manifest.json', dict(schema_version=SCHEMA, catalog_scope=args.catalog_scope,
                  units=units, assembly_manifest_sha256=pf.file_hash(Path(args.manifest)), eligibility_policy=POLICY,
                  legacy_files=legacy_files, orientation_policy='source_orientation_no_RC_or_DTR_or_rotation_normalization',
                  coordinate_system='0-based-half-open_nonwrapping', sequences=len(seqs), occurrences=len(occurrences),
                  real_data_validated=False, scientifically_calibrated=False))


def cluster(args):
    source = Path(args.catalog)
    out = guard(args.outdir, [source, args.blast])
    rows = table(source / 'viral_sequence_master.tsv', SEQ_FIELDS)
    seqs = {k: s for k, (_, s) in pf.read_fasta(source / 'viral_sequences.fna', allow_empty=True).items()}
    if {r['viral_sequence_id'] for r in rows} != set(seqs):
        raise ValueError('Sequence/master foreign-key mismatch')
    for row in rows:
        if row['sequence_sha256'] != pf.sequence_hash(seqs[row['viral_sequence_id']]):
            raise ValueError('Changed sequence')
    grouped = cv.parse_blast(args.blast, seqs)
    pairs, directions = cv.pairwise(grouped, seqs, args.ani_percent, args.af_shorter_percent)
    eligible = [dict(r, genome_id='not_applicable') for r in rows if r['viral_catalog_eligibility'] == 'eligible']
    reps, members = cv.cluster(eligible, pairs)
    by_id = {r['candidate_id']: r for r in members}
    for row in rows:
        if row['candidate_id'] not in by_id:
            members.append(dict(candidate_id=row['candidate_id'], representative_id='', votu_id='',
                                membership_status='excluded_viral_eligibility_' + row['viral_catalog_eligibility'],
                                matching_representative_ids=''))
    manifest = json.loads((source / 'discovery_manifest.json').read_text())
    scope = manifest['catalog_scope']
    policy = dict(algorithm=cv.ALGORITHM, representative_rule='defer_short_bridges_recheck_all_final_representatives',
                  ani_percent=args.ani_percent, af_shorter_percent=args.af_shorter_percent,
                  scope=scope, eligibility=POLICY, sequences=rows, pairs=pairs, membership=members)
    version = stable('viral_', policy)
    for m in members:
        m.update(catalog_version=version, catalog_scope=scope, viral_sequence_id=m['candidate_id'])
    members.sort(key=lambda r: r['candidate_id'])
    by_id = {r['candidate_id']: r for r in members}
    occurrences = table(source / 'source_occurrences.tsv')
    mappings = [dict(o, **{k: by_id[o['viral_sequence_id']][k] for k in
                    ['catalog_version', 'catalog_scope', 'votu_id', 'representative_id', 'membership_status', 'matching_representative_ids']}) for o in occurrences]
    for p in source.iterdir():
        if p.is_file():
            shutil.copy2(p, out / p.name)
    if (source / 'legacy').is_dir():
        shutil.copytree(source / 'legacy', out / 'legacy')
    write(out / 'vOTU_members.tsv', members, MEMBER_FIELDS)
    write(out / 'vOTU_representatives.tsv', [m for m in members if m['membership_status'] == 'representative'], MEMBER_FIELDS)
    write(out / 'occurrence_votu_membership.tsv', mappings, OCC_FIELDS + MEMBER_FIELDS[:3] + MEMBER_FIELDS[5:])
    write(out / 'votu_pairwise.tsv', pairs, ['candidate_a', 'candidate_b', 'ani_percent', 'af_shorter_percent', 'ani_threshold_percent', 'af_threshold_percent', 'passes', 'status'])
    write(out / 'votu_directional_ani.tsv', directions, ['query_id', 'subject_id', 'ani_percent', 'af_shorter_percent', 'aligned_columns', 'identical_columns', 'query_covered_bp', 'subject_covered_bp', 'accepted_hsps'])
    pf.write_fasta(out / 'vOTU_representatives.fna', [(rep, seqs[rep]) for rep in reps])
    pf.write_json(out / 'catalog_manifest.json', dict(manifest, status='completed', catalog_version=version, version_payload=policy,
                  output_files={p.name: pf.file_hash(p) for p in out.iterdir() if p.is_file()},
                  discovery_samples=sorted({u['biological_sample_id'] for u in manifest['units']} - {'unknown'}),
                  representative_fasta_id_semantics='viral_sequence_id; use representative table to obtain votu_id',
                  mapping_reference_id='created_independently_by_B1', representatives=len(reps),
                  code_sha256=pf.file_hash(Path(__file__)), blast_sha256=pf.file_hash(Path(args.blast))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    r = sub.add_parser('validate-manifest'); r.add_argument('--manifest', required=True); r.add_argument('--outdir', required=True)
    r = sub.add_parser('prepare'); r.add_argument('--assembly-id', required=True); r.add_argument('--fasta', required=True); r.add_argument('--outdir', required=True); r.add_argument('--expected-sha256')
    r = sub.add_parser('assemble'); r.add_argument('--manifest', required=True)
    r.add_argument('--sample-dirs', nargs='*', default=[]); r.add_argument('--prepared-dirs', nargs='*', default=[])
    r.add_argument('--legacy-catalog'); r.add_argument('--legacy-source-manifest'); r.add_argument('--legacy-source-fastas', nargs='*', default=None)
    r.add_argument('--catalog-scope', choices=sorted(SCOPES), required=True); r.add_argument('--outdir', required=True)
    r = sub.add_parser('cluster'); r.add_argument('--catalog', required=True); r.add_argument('--blast', required=True)
    r.add_argument('--outdir', required=True); r.add_argument('--ani-percent', type=float, default=95)
    r.add_argument('--af-shorter-percent', type=float, default=85)
    args = p.parse_args()
    if args.command == 'cluster' and not (50 <= args.ani_percent <= 100 and 1 < args.af_shorter_percent <= 100):
        p.error('ANI/AF must use percentages, e.g. 95/85')
    {'validate-manifest': validate_manifest, 'prepare': prepare, 'assemble': assemble, 'cluster': cluster}[args.command](args)


if __name__ == '__main__':
    main()
