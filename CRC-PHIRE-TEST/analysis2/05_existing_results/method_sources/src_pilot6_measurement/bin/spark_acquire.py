#!/usr/bin/env python3
"""Native acquisition receipts and narrow generated-R0 input validation.

This acquisition entry is not a Nextflow workflow. It never invents a session UUID.
"""
import argparse
import csv
import datetime
import fcntl
import json
from pathlib import Path
import re
import signal
import sys
import uuid

import reads_io
import spark_runtime as runtime


def table(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact_hash(path, digest):
    path = Path(path).resolve()
    require(path.is_file() and runtime.sha(path) == digest, 'Generated/source file checksum mismatch: ' + str(path))
    return path


def validate_declared_reads_inputs(sheet, out):
    """Run before prepare: user-provided local inputs never live in this output."""
    sheet = Path(sheet).resolve()
    inputs = [sheet]
    for row in table(sheet):
        for field in ('reads_1', 'reads_2', 'genbank'):
            if row.get(field):
                path = Path(row[field]); inputs.append(path if path.is_absolute() else sheet.parent / path)
    runtime.validate_output_inputs(out, inputs)


def validate_prepared_reads(sheet, out, raw, merge_map=None):
    """Allow only receipt-bound output files generated/reused by this prepare call.

    The raw manifest binds exact sheet/merge inputs; each download and each merged
    FASTQ has its own checksum/source receipt. No general output-subtree exemption.
    """
    sheet = Path(sheet).resolve(); out = Path(out).resolve(); root = out / 'r0'
    validate_declared_reads_inputs(sheet, out)
    require(raw == runtime.read(root / 'raw_reads_manifest.json'), 'Returned/stored prepare manifest differs')
    require(raw.get('samplesheet_sha256') == runtime.sha(sheet), 'Prepare samplesheet checksum mismatch')
    require(raw.get('merge_map_sha256') == (runtime.sha(merge_map) if merge_map else None), 'Prepare merge-map checksum mismatch')
    rows = {row['readset_id']: row for row in table(sheet)}
    require(len(rows) == len(table(sheet)), 'Duplicate readset in samplesheet')
    source_hashes = {}; generated = {}
    for rid, row in rows.items():
        sources = {}
        if row.get('run_accession'):
            accession = row['run_accession']; folder = root / 'downloads' / accession
            require(re.fullmatch(r'[SED]RR\d+', accession), 'Invalid declared run accession')
            acquisition = runtime.read(folder / 'acquisition_manifest.json')
            require(acquisition.get('accession') == accession and acquisition.get('layout') == row['layout'], 'Acquisition receipt identity mismatch')
            files = acquisition.get('files', [])
            by_path = {str(Path(f['path']).resolve()): f['sha256'] for f in files}
            selected = acquisition.get('primary', []) + acquisition.get('singletons', [])
            require(selected and len(selected) == len(set(selected)), 'Acquisition primary/singleton list invalid')
            for name in selected:
                path = Path(name).resolve()
                require(path.is_relative_to(folder.resolve()) and path.is_relative_to(root.resolve()), 'Acquisition path escapes exact run folder')
                digest = by_path.get(str(path)); exact_hash(path, digest)
                generated[str(path)] = digest
                if name in acquisition['primary']: sources[str(path)] = digest
        else:
            for field in ('reads_1', 'reads_2'):
                if row.get(field):
                    path = Path(row[field]); path = (path if path.is_absolute() else sheet.parent / path).resolve()
                    sources[str(path)] = runtime.sha(path)
        source_hashes[rid] = sources
    for unit in raw.get('units', []):
        ids = unit.get('source_readsets', [])
        require(ids and all(rid in rows for rid in ids), 'Prepared unit has unknown source readset')
        merged = unit.get('merge_source_checksums')
        if merged is not None:
            expected = {rid: source_hashes[rid] for rid in ids}
            require(merged == expected, 'Merged source checksums differ from declared/downloaded inputs')
        else:
            require(len(ids) == 1 and unit.get('checksums') == source_hashes[ids[0]], 'Prepared readset was not bound to its own declared source')
        for name, digest in {**unit.get('checksums', {}), **unit.get('source_singleton_checksums', {})}.items():
            path = exact_hash(name, digest)
            if path.is_relative_to(out):
                if str(path) in generated:
                    require(generated[str(path)] == digest, 'Downloaded checksum mismatch')
                elif merged is not None:
                    matched = next((mate for mate in ('reads_1', 'reads_2') if unit.get(mate) and Path(unit[mate]).resolve() == path), None)
                    expected_path = root / 'merged' / (unit['unit_id'] + '.' + str(matched) + '.fastq.gz')
                    require(matched and path == expected_path.resolve() and path.is_relative_to((root / 'merged').resolve()), 'Unrecognized internal merged input')
                    receipt = runtime.read(expected_path.with_suffix('.merge.json'))
                    expected_inputs = []
                    for rid in sorted(ids):
                        source = rows[rid]
                        if source.get('run_accession'):
                            acq = runtime.read(root / 'downloads' / source['run_accession'] / 'acquisition_manifest.json')
                            source_path = acq['primary'][0 if matched == 'reads_1' else 1]
                        else:
                            source_path = Path(source[matched]); source_path = source_path if source_path.is_absolute() else sheet.parent / source_path
                        expected_inputs.append({'readset_id': rid, 'source_sha256': source_hashes[rid][str(Path(source_path).resolve())]})
                    require(receipt.get('output_sha256') == digest and receipt.get('source_manifest') ==
                            {'inputs': expected_inputs, 'namespace_policy': 'readset_prefix_v1', 'gzip_mtime': 0}, 'Merge receipt not bound to exact sources')
                else:
                    raise ValueError('Prepared input inside output lacks recognized acquisition/merge receipt: ' + str(path))
            else:
                runtime.validate_output_inputs(out, [path])
    return True


def verify_genomes(out, sheet, expected_hashes=None):
    out = Path(out).resolve(); rows = table(out / 'genomes.tsv')
    source = json.loads((out / 'source_manifest.json').read_text())
    require(isinstance(source, list) and rows == source, 'Genome table and source manifest disagree')
    expected = table(sheet)
    require(expected and len(rows) == len(expected), 'Acquisition is empty or lacks requested genomes')
    require(len({r['genome_id'] for r in rows}) == len(rows) and len({r['assembly_accession'] for r in rows}) == len(rows), 'Duplicate genome/accession identity')
    require({(r['genome_id'], r['assembly_accession']) for r in rows} ==
            {(r['genome_id'], r['assembly_accession']) for r in expected}, 'Acquired identity/version differs from explicit sheet')
    files = {}
    for row in rows:
        acc = row['assembly_accession']
        require(re.fullmatch(r'GC[AF]_\d+\.\d+', acc), 'Assembly version not explicit')
        folder = out / acc
        for field, digest_field in [('fasta', 'sha256'), ('genbank', None)]:
            if not row.get(field): continue
            path = Path(row[field]).resolve()
            require(path.is_relative_to(folder) and path.is_file(), 'Acquired file escapes exact accession folder')
            digest = runtime.sha(path)
            if digest_field: require(digest == row[digest_field], 'Acquired FASTA checksum differs')
            files[str(path)] = digest
        archive = folder / 'ncbi_dataset.zip'
        exact_hash(archive, row['source_archive_sha256']); files[str(archive)] = row['source_archive_sha256']
    for name in ('genomes.tsv', 'source_manifest.json'): files[str(out / name)] = runtime.sha(out / name)
    if expected_hashes is not None: require(files == expected_hashes, 'Acquisition output hashes changed')
    return files


def verify_completed_acquisition(out):
    out = Path(out).resolve(); path = out / 'pipeline_info/acquisition_state.json'
    state = runtime.read(path)
    require(state.get('status') == 'completed' and state.get('exit_code') == 0 and state.get('nextflow_workflow') is False,
            'Acquisition did not complete')
    require(runtime.sha(state['samplesheet']) == state['samplesheet_sha256'], 'Acquisition sheet changed')
    require(bool(state.get('output_sha256')), 'Acquisition lacks recorded output hashes')
    verify_genomes(out, state['samplesheet'], state.get('output_sha256'))
    return state


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    runtime.allocation()
    require(runtime.enabled(), 'Use bin/spark-run fetch-genomes inside the approved allocation')
    sheet = args.samplesheet.resolve(); out = args.outdir.resolve()
    runtime.validate_output_inputs(out, [sheet])
    info = out / 'pipeline_info'; state_path = info / 'acquisition_state.json'
    previous = runtime.read(state_path) if state_path.exists() else None
    if previous:
        require(args.resume, 'Prior acquisition exists; use --resume or a new outdir')
        require(previous.get('samplesheet') == str(sheet) and previous.get('samplesheet_sha256') == runtime.sha(sheet), 'Resume cannot change acquisition inputs')
    else:
        require(not args.resume, 'Resume needs an existing acquisition receipt')
        require(not out.exists() or not any(out.iterdir()), 'Protect unbound pre-existing acquisition output')
    info.mkdir(parents=True, exist_ok=True)
    with (out / '.acquisition.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
        if previous: runtime.save(info / 'attempts' / (stamp + '.previous.json'), previous)
        state = dict(schema_version=1, target='SPARK_PORTABILITY_BUILD', status='running', nextflow_workflow=False,
                     execution_scope='explicit_versioned_genome_acquisition_not_analysis',
                     samplesheet=str(sheet), samplesheet_sha256=runtime.sha(sheet),
                     real_data_validated=False, scientifically_calibrated=False, resume_requested=args.resume,
                     code_sha256=runtime.sha(__file__), reads_io_sha256=runtime.sha(Path(reads_io.__file__)))
        runtime.save(state_path, state)
        def interrupted(signum, frame):
            raise RuntimeError('Acquisition interrupted by signal ' + str(signum))
        prior_signals = {signum: signal.signal(signum, interrupted) for signum in (signal.SIGINT, signal.SIGTERM)}
        try:
            runtime.native_record(info, ['reads-fetch'])
            if previous and previous.get('status') == 'completed':
                files = verify_genomes(out, sheet, previous.get('output_sha256'))
                state['reused_verified_outputs'] = True
            else:
                reads_io.genomes(sheet, out)
                files = verify_genomes(out, sheet)
                state['reused_verified_outputs'] = False
            state.update(status='completed', exit_code=0, output_sha256=files,
                         genome_count=len(table(out / 'genomes.tsv')))
        except Exception as exc:
            state.update(status='failed', exit_code=1, error=str(exc))
        finally:
            runtime.save(state_path, state)
            for signum, handler in prior_signals.items(): signal.signal(signum, handler)
        print(json.dumps(state, indent=2))
        return state['exit_code']


if __name__ == '__main__':
    try: sys.exit(main())
    except Exception as exc:
        print('Acquisition error: ' + str(exc), file=sys.stderr); sys.exit(1)
