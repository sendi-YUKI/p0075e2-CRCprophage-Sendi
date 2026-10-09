"""Add provisional occurrence-level prophage evidence without rerunning callers.

Requires the accepted ten-donor report and catalog. Source sequence coordinates,
quality estimates and flanks are kept as separate evidence dimensions. Neither a
vOTU membership nor a reference-genome occurrence establishes integration in a
different donor. All production work must run inside Slurm.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import math
import os

STATUS = 'completed'
FLANK_BP = 5000
UNKNOWN = {'', 'na', 'nan', 'none', 'unknown', 'not_assessed', 'not-determined'}
RC = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')
QUALITY_FIELDS = ['contig_length', 'touches_contig_start', 'touches_contig_end',
                  'caller_topology', 'virus_score', 'fdr', 'n_genes', 'n_hallmarks',
                  'marker_enrichment', 'integration_evidence', 'boundary_confidence',
                  'checkv_status', 'checkv_provirus', 'checkv_proviral_length',
                  'checkv_gene_count', 'checkv_viral_genes', 'checkv_host_genes',
                  'checkv_quality', 'checkv_miuvig_quality', 'completeness',
                  'completeness_method', 'contamination', 'checkv_kmer_freq', 'checkv_warnings']
OCC_FIELDS = ['source_occurrence_id', 'viral_sequence_id', 'candidate_id', 'assembly_id',
              'biological_sample_id', 'votu_id', 'membership_status', 'source_kind',
              'candidate_type', 'original_contig_id', 'contig_id', 'start0', 'end0',
              'strand', 'length', 'sequence_sha256', 'source_fasta', 'source_fasta_sha256',
              'source_slice_verification', *QUALITY_FIELDS,
              'left_flank_start0', 'left_flank_end0', 'left_flank_bp', 'left_flank_sha256',
              'right_flank_start0', 'right_flank_end0', 'right_flank_bp', 'right_flank_sha256',
              'flank_annotation_status', 'flank_annotation_file', 'left_flank_gene_count',
              'right_flank_gene_count', 'left_host_marker_count', 'right_host_marker_count',
              'left_host_marker_evidence', 'right_host_marker_evidence',
              'two_sided_host_marker_evidence', 'provisional_evidence_label',
              'att_site_evidence', 'junction_read_evidence', 'inducibility_evidence']
GENE_FIELDS = ['source_occurrence_id', 'viral_sequence_id', 'candidate_id', 'assembly_id',
               'side', 'source_gene_id', 'contig_id', 'original_contig_id', 'gene_start0',
               'gene_end0', 'strand', 'flank_overlap_bp', 'fully_within_flank',
               'marker', 'specificity_class', 'spm_c', 'spm_p', 'spm_v', 'uscg', 'metadata_uscg',
               'uscg_observed', 'host_marker_rule_passed',
               'host_marker_rule', 'source_annotation_file', 'source_annotation_sha256',
               'raw_annotation_json']
VOTU_FIELDS = ['votu_id', 'has_reference_predicted_provirus',
               'has_metagenome_predicted_provirus', 'has_metagenome_two_sided_host_markers',
               'evidence_label', 'n_reference_provirus_occurrences',
               'n_metagenome_provirus_occurrences', 'n_metagenome_two_sided_marker_occurrences',
               'n_source_occurrences', 'n_distinct_sequences', 'source_metagenome_donors',
               'evidence_transfer_policy']


def require(value, message):
    if not value:
        raise RuntimeError(message)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def seqsha(sequence):
    return hashlib.sha256(sequence.upper().encode('ascii')).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def checked_path(root, relative):
    path = root / relative
    require(path.resolve().is_relative_to(root.resolve()), 'Unsafe relative output path: ' + str(path))
    return path


def reuse_completed(out, marker_metadata_path=None):
    """Revalidate a completed result without rewriting it or extracting flanks."""
    receipt_path = out / 'evidence_manifest.json'
    require(receipt_path.is_file(),
            'Existing evidence directory is partial (no completion receipt); inspect before explicit recovery: ' + str(out))
    receipt = read(receipt_path)
    require(receipt.get('status') == STATUS and receipt.get('analysis_type') == 'prophage_evidence',
            'Existing evidence is not a completed compatible result; inspect before recovery: ' + str(out))
    require(receipt.get('script_sha256') == sha(__file__),
            'Evidence script changed since completion; preserve existing result and inspect before a new run')
    if marker_metadata_path is not None:
        previous_marker = receipt.get('marker_metadata', {}).get('path')
        require(previous_marker and Path(previous_marker).resolve() == Path(marker_metadata_path).resolve(),
                'Requested marker metadata path differs from completed evidence; refusing silent parameter changes')
    require(bool(receipt.get('input_sha256')) and bool(receipt.get('outputs')),
            'Existing evidence receipt lacks input or output checksums; refusing reuse')
    for name, expected in receipt['input_sha256'].items():
        path = Path(name)
        require(path.is_absolute() and path.is_file(), 'Evidence input is missing or not absolute: ' + name)
        require(sha(path) == expected,
                'Evidence input changed; refusing reuse or overwrite: ' + name)
    for relative, expected in receipt['outputs'].items():
        path = checked_path(out, relative)
        require(path.is_file(), 'Completed evidence output is missing: ' + str(path))
        require(sha(path) == expected,
                'Evidence output changed; refusing reuse or overwrite: ' + str(path))
    print(json.dumps(dict(status='reused_verified_completed_evidence', output=str(out),
                          occurrences=receipt['occurrences'], votus=receipt['votus'])), flush=True)
    return receipt


def load_required_contigs(path, needed):
    """Retain only relevant contigs while checking duplicate requested IDs."""
    with Path(path).open('rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    result = {}
    identifier = None
    parts = []

    def finish():
        if identifier in needed:
            require(identifier not in result, 'Duplicate requested FASTA identifier: ' + identifier)
            sequence = ''.join(parts).upper()
            require(bool(sequence), 'Empty requested FASTA sequence: ' + identifier)
            result[identifier] = sequence

    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                finish()
                identifier = line[1:].split()[0]
                parts = []
            elif identifier in needed:
                parts.append(line)
        finish()
    require(set(result) == set(needed), 'Missing source contigs: ' + ','.join(sorted(set(needed)-set(result))[:5]))
    return result


def verify_slice_and_flanks(source, occurrence, master):
    start, end = int(occurrence['start0']), int(occurrence['end0'])
    require(0 <= start < end <= len(source), 'Source slice outside contig')
    require(int(master['contig_length']) == len(source), 'Source contig length differs from primary master')
    for key in ['candidate_type', 'start0', 'end0', 'strand', 'sequence_sha256']:
        require(str(occurrence[key]) == str(master[key]), 'Primary/occurrence mismatch: ' + key)
    sequence = source[start:end]
    if occurrence['strand'] == '-':
        sequence = sequence.translate(RC)[::-1]
    else:
        require(occurrence['strand'] == '+', 'Unsupported source orientation')
    require(seqsha(sequence) == occurrence['sequence_sha256'], 'Candidate/source slice SHA256 mismatch')
    require(len(sequence) == int(occurrence['length']) == int(master['length']), 'Candidate length mismatch')
    require(master['touches_contig_start'] == str(start == 0).lower()
            and master['touches_contig_end'] == str(end == len(source)).lower(), 'Contig-end flags differ from source')
    left_start, right_end = max(0, start-FLANK_BP), min(len(source), end+FLANK_BP)
    return dict(left=(left_start, start, source[left_start:start]),
                right=(end, right_end, source[end:right_end]))


def bool_any(values, empty='unknown'):
    values = list(values)
    if 'true' in values:
        return 'true'
    if not values or 'unknown' in values:
        return empty if not values else 'unknown'
    require(set(values) <= {'false'}, 'Unexpected evidence truth values')
    return 'false'


def combine_sides(left, right):
    if left == right == 'true':
        return 'true'
    # An unavailable or zero-length flank does not establish a negative result.
    if 'unknown' in (left, right):
        return 'unknown'
    return 'false'


def summarize_votus(occurrences, memberships):
    by_votu = defaultdict(list)
    for row in occurrences:
        if row['votu_id']:
            by_votu[row['votu_id']].append(row)
    output = []
    for votu in sorted({m['votu_id'] for m in memberships if m['votu_id']}):
        rows = by_votu[votu]
        require(bool(rows), 'vOTU lacks source occurrences: ' + votu)
        ref = [r for r in rows if r['source_kind'] == 'legacy_bacterial_genome' and r['candidate_type'] == 'provirus_locus']
        meta = [r for r in rows if r['source_kind'] == 'metagenome_assembly' and r['candidate_type'] == 'provirus_locus']
        marker = bool_any(r['two_sided_host_marker_evidence'] for r in meta)
        label = ('provisional_metagenome_provirus_with_caller_flank_markers' if marker == 'true'
                 else 'provisional_metagenome_predicted_provirus' if meta
                 else 'provisional_reference_predicted_provirus_only' if ref
                 else 'integration_unknown_viral_sequence')
        output.append(dict(votu_id=votu, has_reference_predicted_provirus=str(bool(ref)).lower(),
                           has_metagenome_predicted_provirus=str(bool(meta)).lower(),
                           has_metagenome_two_sided_host_markers=marker, evidence_label=label,
                           n_reference_provirus_occurrences=len(ref), n_metagenome_provirus_occurrences=len(meta),
                           n_metagenome_two_sided_marker_occurrences=sum(r['two_sided_host_marker_evidence']=='true' for r in meta),
                           n_source_occurrences=len(rows), n_distinct_sequences=len({r['viral_sequence_id'] for r in rows}),
                           source_metagenome_donors=';'.join(sorted({r['biological_sample_id'] for r in meta})),
                           evidence_transfer_policy='source_occurrence_only_not_integration_in_every_mapped_donor'))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--marker-metadata', type=Path,
                        default=Path('/srv/CRC-PHIRE/databases/genomad/1.9/genomad_db/genomad_marker_metadata.tsv'))
    args = parser.parse_args()
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Production evidence runs require Slurm')
    root = args.root.resolve()
    report = root / '12_reports/pilot10'
    catalog = root / '06_viral_catalog/pilot10_combined'
    out = report / 'prophage_evidence'
    if out.exists():
        reuse_completed(out, args.marker_metadata)
        return
    inputs = {}

    def record(path):
        path = Path(path).resolve()
        digest = sha(path)
        previous = inputs.get(str(path))
        require(previous in (None, digest), 'Input changed during evidence collection: ' + str(path))
        inputs[str(path)] = digest
        return digest

    def verify_manifest(base, values):
        require(bool(values), 'Empty upstream checksum manifest')
        for relative, digest in values.items():
            path = checked_path(base, relative)
            require(record(path) == digest, 'Upstream checksum mismatch: ' + str(path))

    acceptance = report / 'pilot10_final_acceptance.json'
    accepted = read(acceptance)
    require(accepted['status'] == 'completed_balanced_engineering_pilot' and accepted['donors'] == 10,
            'Ten-donor final report is not accepted')
    verify_manifest(report, accepted['output_sha256'])
    record(acceptance)
    state_path = catalog / 'pipeline_info/stage_state.json'
    from pilot10_catalog_acceptance import verified_catalog_completion
    publication = verified_catalog_completion(root, catalog)
    if publication is not None:
        record(publication)
    record(state_path)
    terminal = catalog / 'pipeline_info/terminal_output_manifest.json'
    terminal_outputs = read(terminal)
    verify_manifest(catalog, terminal_outputs)
    record(terminal)
    required = ['catalog/source_occurrences.tsv', 'catalog/occurrence_votu_membership.tsv', 'catalog/vOTU_members.tsv']
    require(all(path in terminal_outputs for path in required), 'Evidence catalog tables are not in terminal acceptance')
    occurrences = table(catalog / required[0])
    membership_rows = table(catalog / required[1])
    members = table(catalog / required[2])
    membership = {r['source_occurrence_id']: r for r in membership_rows}
    require(len(membership) == len(membership_rows) == len(occurrences)
            and set(membership) == {r['source_occurrence_id'] for r in occurrences}, 'Occurrence/membership identity mismatch')
    by_assembly = defaultdict(list)
    for row in occurrences:
        require(row['source_kind'] in {'legacy_bacterial_genome', 'metagenome_assembly'}, 'Unexpected occurrence source kind')
        require(row['slice_verification'] == 'verified_this_run', 'Catalog occurrence lacks source slice verification')
        by_assembly[row['assembly_id']].append(row)

    reference = root / '03_reference/reference258_raw_v1'
    integration_path = reference / 'integration_manifest.json'
    integration = read(integration_path)
    verify_manifest(reference, integration['files'])
    record(integration_path)
    reference_master = table(reference / 'prophage_master.tsv')
    require(len(reference_master) == 812, 'Unexpected frozen reference candidate count')
    require(sum(r['source_kind']=='legacy_bacterial_genome' for r in occurrences) == len(reference_master),
            'Catalog lost or duplicated reference occurrences')
    ref_index = {(r['genome_id'], r['candidate_id']): r for r in reference_master}
    require(len(ref_index) == len(reference_master), 'Duplicate reference candidates')
    genome_path = root / '01_manifests/genome_manifest.tsv'
    genomes = {r['genome_id']: r for r in table(genome_path)}
    require(len(genomes) == 258, 'Reference panel is not the frozen 258 genomes')
    record(genome_path)
    imported_path = root / '01_manifests/pilot10_imported_primary_assemblies.json'
    imported = read(imported_path)
    units = {u['assembly_id']: u for u in imported['units']}
    require(len(units) == len(imported['units']) == 10, 'Expected ten patient assemblies')
    record(imported_path)
    # Nextflow sorts/normalizes and JSON-reserializes the assembly manifest.
    # Its staged byte hash is therefore not the external JSON file hash.
    catalog_units = {u['assembly_id']: u for u in read(catalog / 'catalog/catalog_manifest.json')['units']}
    require(set(catalog_units) == set(units), 'Catalog/imported assembly identities differ')
    for aid, unit in units.items():
        require(all(catalog_units[aid][key] == unit[key] for key in
                    ['sha256', 'biological_sample_id', 'input_kind', 'material_type', 'primary_final_dir', 'primary_prepared_dir']),
                'Imported primary assembly differs from the accepted catalog source: ' + aid)
    frozen_path = root / '01_manifests/pilot10_frozen_assemblies.json'
    frozen = read(frozen_path)
    require(Counter(frozen['groups'].values()) == {'CRC':5, 'control':5}
            and {u['biological_sample_id'] for u in imported['units']} == set(frozen['groups']),
            'Evidence inputs differ from the ten-donor five-versus-five freeze')
    require({u['assembly_id']: u['sha256'] for u in frozen['units']}
            == {aid: unit['sha256'] for aid, unit in units.items()}, 'Frozen/imported assembly SHA256 mismatch')
    record(frozen_path)
    import_acceptance_path = root / '13_handoff/pilot10_primary_import_acceptance.json'
    import_acceptance = read(import_acceptance_path)
    require(import_acceptance['status'] == 'primary_files_verified_pending_collector_source_slice_checks',
            'Primary import receipt is missing')
    imported_files = {}
    for entry in import_acceptance['primary_sources']:
        imported_files.update(entry['files'])
    record(import_acceptance_path)
    marker_metadata, marker_metadata_info = load_marker_metadata(args.marker_metadata, record)

    out.mkdir()
    evidence_rows, gene_rows, annotation_receipts = [], [], []
    try:
        with gzip.open(out / 'source_flanks.fna.gz', 'wt', encoding='ascii', newline='\n') as flank_output:
            for aid, rows in sorted(by_assembly.items()):
                reference_side = rows[0]['source_kind'] == 'legacy_bacterial_genome'
                require(all((r['source_kind']=='legacy_bacterial_genome') == reference_side for r in rows), 'Mixed assembly source kinds')
                if reference_side:
                    source = Path(genomes[aid]['fasta'])
                    expected_source_sha = genomes[aid]['sha256']
                    master_index = {candidate: data for (genome, candidate), data in ref_index.items() if genome == aid}
                    gene_roots = [Path(s['source']) / f'samples/{aid}/raw/genomad' for s in integration['sources']]
                else:
                    unit = units[aid]
                    source = Path(unit['assembly_fasta'])
                    expected_source_sha = unit['sha256']
                    primary_master = Path(unit['primary_final_dir']) / 'prophage_master.tsv'
                    master_rows = table(primary_master)
                    master_index = {r['candidate_id']: r for r in master_rows}
                    require(len(master_index) == len(master_rows), 'Duplicate metagenome candidate IDs')
                    require(record(primary_master) == imported_files.get(str(primary_master)),
                            'Primary master changed since its accepted import: ' + str(primary_master))
                    primary_base = Path(unit['primary_final_dir']).parents[3]
                    gene_roots = [primary_base / f'samples/{aid}/raw/genomad', primary_base / f'assemblies/{aid}/geNomad']
                source_hash = record(source)
                require(source_hash == expected_source_sha, 'Original source FASTA has changed: ' + aid)
                source_contigs = load_required_contigs(source, {r['original_contig_id'] for r in rows})
                annotations, annotation_info = load_annotations(gene_roots, {r['contig_id'] for r in rows}, record)
                annotation_receipts.append(dict(assembly_id=aid, **annotation_info))
                for occurrence in rows:
                    master = master_index[occurrence['candidate_id']]
                    require(master['genome_id'] == aid, 'Primary master assembly mismatch')
                    contig = source_contigs[occurrence['original_contig_id']]
                    require(all(gene['end0'] <= len(contig) for gene in annotations.get(occurrence['contig_id'], [])),
                            'Full-contig annotation coordinates exceed original contig length')
                    flanks = verify_slice_and_flanks(contig, occurrence, master)
                    mem = membership[occurrence['source_occurrence_id']]
                    require(mem['viral_sequence_id'] == occurrence['viral_sequence_id'], 'Membership sequence mismatch')
                    ev = {field: occurrence.get(field, 'unknown') for field in OCC_FIELDS}
                    ev.update({field: master.get(field, 'unknown') for field in QUALITY_FIELDS})
                    local_annotation = dict(annotation_info)
                    if occurrence['contig_id'] not in annotations and annotation_info['status'] == 'assessed_caller_annotations':
                        local_annotation['status'] = 'unknown_contig_not_in_annotation'
                    if not marker_metadata and local_annotation['status'] == 'assessed_caller_annotations':
                        local_annotation['status'] = 'unknown_marker_metadata_not_available'
                    ev.update(votu_id=mem['votu_id'], membership_status=mem['membership_status'],
                              source_fasta=str(source), source_fasta_sha256=source_hash,
                              source_slice_verification='verified_this_run',
                              att_site_evidence='not_assessed', junction_read_evidence='not_assessed',
                              inducibility_evidence='not_assessed', flank_annotation_status=local_annotation['status'],
                              flank_annotation_file=annotation_info.get('path', ''))
                    for side, (start, end, sequence) in flanks.items():
                        ev.update({side+'_flank_start0': start, side+'_flank_end0': end,
                                   side+'_flank_bp': len(sequence), side+'_flank_sha256': seqsha(sequence) if sequence else ''})
                        if sequence:
                            flank_output.write('>'+occurrence['source_occurrence_id']+'|'+side+'|'+str(start)+'-'+str(end)+'\n'+sequence+'\n')
                        selected = flank_gene_evidence(occurrence, side, start, end,
                                                       annotations.get(occurrence['contig_id'], []), annotation_info, marker_metadata)
                        gene_rows.extend(selected)
                        assessed = local_annotation['status'] == 'assessed_caller_annotations' and end > start
                        nmarkers = sum(r['host_marker_rule_passed'] == 'true' for r in selected)
                        ev[side+'_flank_gene_count'] = len(selected) if assessed else 'unknown'
                        ev[side+'_host_marker_count'] = nmarkers if assessed else 'unknown'
                        ev[side+'_host_marker_evidence'] = (bool_any((r['host_marker_rule_passed'] for r in selected), empty='false')
                                                          if assessed else 'unknown')
                    ev['two_sided_host_marker_evidence'] = combine_sides(ev['left_host_marker_evidence'], ev['right_host_marker_evidence'])
                    ev['provisional_evidence_label'] = ('caller_provirus_with_two_sided_caller_marker_context' if occurrence['candidate_type']=='provirus_locus' and ev['two_sided_host_marker_evidence']=='true'
                                                       else 'caller_predicted_provirus_only' if occurrence['candidate_type']=='provirus_locus'
                                                       else 'viral_sequence_integration_unknown')
                    evidence_rows.append(ev)
        evidence_rows.sort(key=lambda r: r['source_occurrence_id'])
        votus = summarize_votus(evidence_rows, members)
        write_table(out / 'occurrence_evidence.tsv', evidence_rows, OCC_FIELDS)
        write_table(out / 'votu_evidence.tsv', votus, VOTU_FIELDS)
        write_table(out / 'flank_genes.tsv', gene_rows, GENE_FIELDS)
        manifest = dict(schema_version=1, status=STATUS, analysis_type='prophage_evidence', created_utc=datetime.now(timezone.utc).isoformat(),
                        job_id=os.environ['SLURM_JOB_ID'], input_sha256=inputs, script_sha256=sha(__file__),
                        outputs={p.name: sha(p) for p in out.iterdir() if p.is_file()},
                        occurrences=len(evidence_rows), votus=len(votus), flank_genes=len(gene_rows),
                        occurrence_labels=dict(Counter(r['provisional_evidence_label'] for r in evidence_rows)),
                        annotation_sources=annotation_receipts, marker_rule=MARKER_RULE,
                        marker_metadata=marker_metadata_info,
                        flank_window_bp=FLANK_BP, flank_coordinates='0-based_half-open_original_contig_orientation',
                        evidence_dimensions=['caller_integration_prediction', 'CheckV_quality_and_completeness', 'contig_end_truncation', 'geometric_flanks', 'caller_flank_markers'],
                        not_assessed=['attL_attR', 'junction_read_support', 'inducibility', 'infectivity', 'independent_host_flank_taxonomy'],
                        limitations=[
                            '541 reference provirus instances are caller predictions, not confirmed inducible phages; this analysis does not change that status.',
                            'Geometric flanks are not automatically bacterial sequence. Marker context derives from the same geNomad annotation and is not independent caller validation.',
                            'A missing annotation table or unavailable flank is unknown, not absence. False marker flags mean no qualifying marker observed in an assessed window only.',
                            'CheckV No is not a veto of integration: CheckV input is an already extracted candidate which may lack host flanks.',
                            'CheckV Complete or high completeness does not demonstrate induction, viable particles or infectivity.',
                            'A vOTU evidence label summarizes its source occurrences. It must not be transferred into an integration claim for each donor whose reads map to that vOTU.',
                            'No caller, gene predictor, mapping task or database download is performed by this evidence layer.',
                        ])
        (out / 'evidence_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        print(json.dumps(dict(status=STATUS, occurrences=len(evidence_rows), votus=len(votus), output=str(out))), flush=True)
    except Exception as error:
        (out / 'failed_evidence.json').write_text(json.dumps(dict(status='failed_preserved', error=str(error), input_sha256=inputs), indent=2)+'\n', encoding='utf-8')
        raise


# The annotation parser is deliberately separate so the field semantics and
# coordinate conversion can be checked against the installed geNomad schema.
MARKER_RULE = ('Fully contained gene in the available 0-5 kb original-contig flank, with an observed geNomad marker whose '
               'current pinned metadata has SPM_C > SPM_P and SPM_C > SPM_V; strict maximum, no universal cutoff. '
               'USCG is retained independently from the geNomad gene flag and metadata USCG identifier; it does not substitute '
               'for the SPM rule or establish taxonomy. Project reporting rule using the same caller, not independent validation.')


def load_marker_metadata(path, record):
    if not path.is_file():
        return {}, dict(status='unknown_missing_marker_metadata', path=str(path))
    digest = record(path)
    metadata = {}
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = ['MARKER', 'SPECIFICITY_CLASS', 'SPM_C', 'SPM_P', 'SPM_V', 'USCG']
        require(set(fields) <= set(reader.fieldnames or []), 'geNomad marker metadata schema is not supported')
        for row in reader:
            require(None not in row and row['MARKER'] not in metadata, 'Malformed or duplicate marker metadata')
            metadata[row['MARKER']] = {field: row[field] for field in fields[1:]}
    version_path = path.parent / 'version.txt'
    info = dict(status='loaded', path=str(path), sha256=digest, markers=len(metadata),
                version=version_path.read_text().strip() if version_path.is_file() else 'unknown',
                version_sha256=record(version_path) if version_path.is_file() else 'unknown',
                interpretation='SPM metadata reports chromosome/plasmid/virus preference; not an independent host taxonomy classifier')
    return metadata, info


def load_annotations(roots, contigs, record):
    # Only the full-genome annotation table has the validated complete header.
    # The similarly named provirus table can have a different schema and must
    # never be treated as a full-contig or flank annotation source.
    paths = sorted({p for root in roots if root.is_dir()
                    for p in root.glob('*/genome_annotate/genome_genes.tsv')})
    info = dict(searched_roots=[str(p) for p in roots])
    if not paths:
        return {}, dict(info, status='unknown_annotation_not_available')
    hashes = {str(path): record(path) for path in paths}
    if len(set(hashes.values())) != 1:
        return {}, dict(info, status='unknown_ambiguous_annotation_sources', candidate_sha256=hashes)
    path = paths[0]
    annotations = defaultdict(list)
    seen = set()
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'gene', 'start', 'end', 'length', 'strand', 'marker', 'uscg'}
        if not required <= set(reader.fieldnames or []):
            return {}, dict(info, status='unknown_unsupported_annotation_schema', path=str(path), sha256=hashes[str(path)])
        for row in reader:
            require(None not in row and all(value is not None for value in row.values()), 'Malformed full geNomad gene annotation row')
            contig, separator, ordinal = row['gene'].rpartition('_')
            if not separator or contig not in contigs:
                continue
            require(ordinal.isdigit() and row['gene'] not in seen, 'Invalid or duplicate source gene identifier')
            seen.add(row['gene'])
            start, end = int(row['start'])-1, int(row['end'])
            require(0 <= start < end and end-start == int(row['length']), 'Invalid geNomad 1-based-inclusive gene coordinates')
            require(row['strand'] in {'1', '-1'}, 'Unsupported geNomad gene strand')
            annotations[contig].append(dict(row, start0=start, end0=end,
                                           normalized_strand='+' if row['strand']=='1' else '-'))
    for genes in annotations.values():
        genes.sort(key=lambda r: (r['start0'], r['end0'], r['gene']))
    return dict(annotations), dict(info, status='assessed_caller_annotations', path=str(path),
                                   sha256=hashes[str(path)], duplicate_sources=hashes,
                                   coordinate_system='input_1-based_inclusive_converted_to_0-based_half-open',
                                   relevant_genes=len(seen), source='same_geNomad_caller_not_independent_validation')


def flank_gene_evidence(occurrence, side, start, end, genes, annotation_info, marker_metadata):
    result = []
    for gene in genes:
        overlap = max(0, min(end, gene['end0'])-max(start, gene['start0']))
        if not overlap:
            continue
        contained = start <= gene['start0'] < gene['end0'] <= end
        marker = gene['marker']
        meta = marker_metadata.get(marker)
        spms = ['unknown']*3
        passed = 'false'
        if marker.lower() not in UNKNOWN:
            if meta is None:
                passed = 'unknown' if contained else 'false'
            else:
                try:
                    spms = [float(meta[key]) for key in ['SPM_C', 'SPM_P', 'SPM_V']]
                    require(all(math.isfinite(value) for value in spms), 'Nonfinite marker specificity values')
                    passed = str(contained and spms[0] > spms[1] and spms[0] > spms[2]).lower()
                except ValueError:
                    spms = ['unknown']*3
                    passed = 'unknown' if contained else 'false'
        result.append(dict(source_occurrence_id=occurrence['source_occurrence_id'],
                           viral_sequence_id=occurrence['viral_sequence_id'], candidate_id=occurrence['candidate_id'],
                           assembly_id=occurrence['assembly_id'], side=side, source_gene_id=gene['gene'],
                           contig_id=occurrence['contig_id'], original_contig_id=occurrence['original_contig_id'],
                           gene_start0=gene['start0'], gene_end0=gene['end0'], strand=gene['normalized_strand'],
                           flank_overlap_bp=overlap, fully_within_flank=str(contained).lower(), marker=marker,
                           specificity_class=meta['SPECIFICITY_CLASS'] if meta else 'unknown',
                           spm_c=spms[0], spm_p=spms[1], spm_v=spms[2], uscg=gene['uscg'],
                           metadata_uscg=meta['USCG'] if meta else 'unknown',
                           uscg_observed='true' if gene['uscg']=='1' else 'false' if gene['uscg']=='0' else 'unknown',
                           host_marker_rule_passed=passed, host_marker_rule=MARKER_RULE,
                           source_annotation_file=annotation_info.get('path', ''),
                           source_annotation_sha256=annotation_info.get('sha256', 'unknown'),
                           raw_annotation_json=json.dumps({k:v for k,v in gene.items() if k not in {'start0','end0','normalized_strand'}}, sort_keys=True)))
    return result


if __name__ == '__main__':
    main()
