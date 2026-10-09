#!/usr/bin/env python3
"""Read-only core readiness inspection. Never launches tools, jobs or analyses."""
import argparse
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
from phageflow import (file_hash, load_samplesheet, read_fasta, read_tsv,
                       unique_index)
import spark_runtime as runtime
from core_genomes import validate_reference_manifest
from spark_fs_guard import reject_ipc_storage

REPO = Path(__file__).resolve().parents[1]
CORE_ENVS = ('python', 'core-genomad', 'core-checkv', 'core-checkm2',
             'core-fastani', 'core-blast', 'core-functions', 'core-callers')
DB_KEYS = ('genomad', 'checkv', 'checkm2', 'functions')
PAUSED = ('iPHoP', 'DRAM-v', 'VirSorter2')


def bounded_json(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size > 16 * 1024 * 1024:
        raise ValueError('Missing or oversized metadata JSON: ' + str(path))
    return json.loads(path.read_text())


def table(path, required):
    if Path(path).stat().st_size > 1024 * 1024:
        raise ValueError('Table exceeds 1 MiB inspection limit')
    return read_tsv(Path(path), required)[1]


def fresh_result_dir(path):
    path = Path(path).expanduser()
    if path.is_symlink() or path.exists():
        raise ValueError('Choose a new, nonexistent directory: ' + str(path))
    path = path.resolve()
    if path.parent != (REPO / 'results').resolve():
        raise ValueError('Use a new direct child of the repository results/ directory')
    reject_ipc_storage([path])
    return path


def within(path, root):
    path = Path(path).expanduser().resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError('Path escapes expected root: ' + str(path))
    reject_ipc_storage([path])
    return path


def nonempty(path, size=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0 or not os.access(path, os.R_OK):
        raise ValueError('Missing/empty/unreadable component: ' + str(path))
    if size is not None and path.stat().st_size != size:
        raise ValueError('Component size differs from manifest: ' + str(path))


def metadata_hash(path, expected):
    bounded_json(path)
    observed = file_hash(Path(path))
    if not expected or observed != expected:
        raise ValueError('Metadata SHA256 mismatch: ' + str(path))
    return observed


def check_database(key, entry):
    root = Path.home() / 'databases'
    if entry.get('status') != 'ready':
        raise ValueError('Database not registered ready: ' + key)
    path = within(entry['path'], root)
    if not path.exists():
        raise ValueError('Missing database path: ' + str(path))
    checked = {}
    for field in ('inventory', 'receipt', 'index_receipt', 'runtime_manifest',
                  'native_profile_receipt'):
        if field in entry:
            p = Path(entry[field])
            if not (p.resolve().is_relative_to(root) or p.resolve().is_relative_to(REPO)):
                raise ValueError('Database evidence is outside Spark project/storage')
            checked[field] = metadata_hash(p, entry.get(field + '_sha256'))
    if not {'inventory', 'receipt'}.issubset(checked):
        raise ValueError('Database inventory/receipt identities missing')
    inventory = bounded_json(entry['inventory'])
    files = inventory.get('files', []) if isinstance(inventory, dict) else inventory
    if isinstance(files, dict):
        files = [dict(path=k, **v) for k, v in files.items()]
    if key == 'checkm2' and not files:
        files = [dict(path=inventory['database_file'],
                      bytes=inventory['database_bytes'],
                      sha256=inventory['database_sha256'])]
        if files != entry['files']:
            raise ValueError('CheckM2 inventory/registry identity mismatch')
        if Path(inventory['actual_path']).resolve() != path:
            raise ValueError('CheckM2 inventory path mismatch')
    if not files:
        raise ValueError('Empty database component inventory')
    # Extraction inventories are rooted at version directories; functions at DB root.
    base = Path(entry['inventory']).parent
    for item in files:
        component = within(base / item['path'], base)
        if not component.is_file() or not os.access(component, os.R_OK):
            raise ValueError('Missing/unreadable inventory file: ' + str(component))
        if 'bytes' not in item or component.stat().st_size != item['bytes']:
            raise ValueError('Component size differs from manifest: ' + str(component))
    if key == 'checkv':
        receipt = bounded_json(entry['index_receipt'])
        nonempty(path / 'genome_db/checkv_reps.dmnd', receipt['bytes'])
        lock = runtime.environment_definition('core-checkv')[1]['lock_sha256']
        if receipt['container'] != 'native-lock:core-checkv:' + lock:
            raise ValueError('CheckV index/environment identity mismatch')
    if key == 'checkm2':
        nonempty(path, entry['files'][0]['bytes'])
    return dict(path=str(path), version=entry['version'], metadata_sha256=checked,
                components_stat_checked=len(files), payload_hashes_recomputed=False,
                scope='metadata_hashes_and_component_presence_size; historical full hashes retained')


def check_platform():
    if platform.system() != 'Linux' or platform.machine() not in {'aarch64', 'arm64'}:
        raise ValueError('Requires Linux ARM64 on Spark')
    if REPO != Path('/srv/CRC-PHIRE/projects/crc-pipeline'):
        raise ValueError('Requires the unique Spark execution repository')
    return dict(platform='linux-aarch64', repository=str(REPO))


def inspect(args):
    report = dict(schema_version=1, mode='core_preflight_only', checks=[],
                  real_data_validated=False, scientifically_calibrated=False,
                  nextflow_started=False, tasks_submitted=0, tools_executed=0,
                  payload_hashes_recomputed=False)
    def add(name, state, message, evidence=None):
        report['checks'].append(dict(check=name, state=state, message=message,
                                     evidence=evidence))
    def attempt(name, action):
        try:
            evidence = action()
            add(name, 'pass', 'Read-only contract checked', evidence)
            return evidence
        except (ValueError, OSError, KeyError, TypeError, RuntimeError, ImportError) as exc:
            add(name, 'blocker', str(exc))
            return None
    attempt('platform', check_platform)
    add('resources', 'info', 'Configured budget; no Slurm allocation requested',
        dict(cpus=6, memory_gib=48, task_memory_gib=46, queue_size=1))
    for key in CORE_ENVS:
        def env_check(key=key):
            _, item, prefix, lock, _ = runtime.environment_definition(key)
            return dict(prefix=str(prefix), lock_sha256=file_hash(lock),
                        tools=item.get('tools'), native_binaries_rehashed=False)
        attempt('environment:' + key, env_check)
    db = attempt('database_registry', lambda: bounded_json(REPO / 'assets/spark_databases.json'))
    if db:
        for key in DB_KEYS:
            attempt('database:' + key, lambda key=key: check_database(key, db['databases'][key]))
    if getattr(args,'enable_kofam',False) or getattr(args,'host_config',None) or getattr(args,'catalog_policy',None) or getattr(args,'host_association_config',None):
        from mainline_preflight import inspect_extensions
        extension=attempt('mainline_extensions',lambda:inspect_extensions(args))
        if extension and extension.get('host_source',{}).get('pending'):
            add('host_source_scientific_choices','pending_decision',','.join(extension['host_source']['pending']))
        if extension and extension.get('catalog_research',{}).get('status')=='configuration_proposed':
            add('catalog_scientific_choices','pending_decision','Catalog/locus evidence preparation only; thresholds and calibration pending')
        if extension and extension['host_phylogeny']=='configuration_proposed':
            add('host_scientific_choices','pending_decision','Host proposed configuration is not authorization or readiness for real execution')
    records = {}
    inputs = []
    if args.input:
        def check_sheet():
            table(args.input, ['genome_id', 'fasta'])
            return load_samplesheet(args.input)
        inputs = attempt('samplesheet', check_sheet) or []
        remaining = args.max_total_bases
        for row in inputs:
            def check_fasta(row=row):
                p = Path(row['fasta'])
                reject_ipc_storage([p])
                if p.stat().st_size > 64 * 1024 * 1024:
                    raise ValueError('FASTA file exceeds 64 MiB inspection limit')
                parsed = read_fasta(p, max_bases=remaining)
                records[row['genome_id']] = {k: v[1] for k, v in parsed.items()}
                return dict(path=str(p), sha256=file_hash(p), contigs=len(parsed),
                            bases=sum(len(v[1]) for v in parsed.values()))
            evidence = attempt('fasta:' + row['genome_id'], check_fasta)
            if evidence:
                remaining -= evidence['bases']
    else:
        add('samplesheet', 'not_provided', 'Environment readiness only; research input not selected')
    ids = {x['genome_id'] for x in inputs}
    protected = [x['fasta'] for x in inputs] + [args.input, args.metadata, args.genbanks,
                                                args.reference_dir, Path.home() / 'databases']
    if args.metadata:
        def check_metadata():
            rows = table(args.metadata, ['genome_id'])
            index = unique_index(rows, 'genome_id', 'genome metadata')
            if set(index) - ids:
                raise ValueError('Metadata includes unknown genome IDs')
            return dict(rows=len(rows), missing_genomes=sorted(ids-set(index)),
                        clinical_independence_verified=False)
        attempt('metadata', check_metadata)
    else:
        add('metadata', 'not_provided', 'Declared taxonomy/context unknown; no clinical labels inferred')
    if args.genbanks:
        def check_genbanks():
            rows = table(args.genbanks, ['genome_id', 'genbank'])
            index = unique_index(rows, 'genome_id', 'GenBank mapping')
            if set(index) - ids:
                raise ValueError('GenBank mapping includes unknown genome IDs')
            checked = []
            for gid, row in index.items():
                if not row['genbank']:
                    checked.append(dict(genome_id=gid, state='not_provided'))
                    continue
                p = Path(row['genbank']).expanduser()
                p = p if p.is_absolute() else args.genbanks.parent/p
                p = p.resolve()
                protected.append(p)
                reject_ipc_storage([p])
                nonempty(p)
                if p.stat().st_size > 64*1024*1024:
                    raise ValueError('GenBank exceeds inspection limit')
                expected = row.get('sha256') or row.get('genbank_sha256')
                actual = file_hash(p)
                if expected and expected != actual:
                    raise ValueError('GenBank checksum mismatch: ' + gid)
                if gid not in records:
                    raise ValueError('Cannot check GenBank without valid matching FASTA: ' + gid)
                from core_callers import validate_genbank
                checked.append(dict(genome_id=gid, sha256=actual,
                                    contigs=validate_genbank(p, records[gid])))
            return dict(rows=checked, missing_genomes=sorted(ids-set(index)))
        attempt('genbanks', check_genbanks)
    else:
        add('genbanks', 'conditional_not_provided', 'PhiSpy not assessed without matching annotated GenBank; core continues')
    if args.reference_dir:
        def check_references():
            root = args.reference_dir.resolve()
            manifest = validate_reference_manifest(bounded_json(root/'reference_manifest.json'))
            refs = manifest['references']
            if not refs:
                return dict(references=0, taxonomy='unresolved', content_checked=False)
            unique_index(refs, 'reference_id', 'reference panel')
            seen = set()
            checked = []
            budget = args.max_total_bases
            for ref in refs:
                for field in ('reference_id', 'fasta', 'fasta_sha256', 'species', 'species_taxid'):
                    if not ref.get(field):
                        raise ValueError('Missing reference field: ' + field)
                # Production directory override stages a flat reference catalog.
                if ref['fasta'] != ref['reference_id'] + '.fna':
                    raise ValueError('Reference FASTA must be <reference_id>.fna')
                p = within(root/ref['fasta'], root)
                if p in seen:
                    raise ValueError('Duplicate reference FASTA')
                seen.add(p)
                nonempty(p)
                if p.stat().st_size > 64*1024*1024:
                    raise ValueError('Reference exceeds inspection limit')
                if file_hash(p) != ref['fasta_sha256']:
                    raise ValueError('Reference checksum mismatch')
                seqs = read_fasta(p, max_bases=budget)
                budget -= sum(len(v[1]) for v in seqs.values())
                checked.append(dict(path=str(p), trust_source=ref.get('trust_source') or None,
                                    eligibility='declared_trusted' if ref.get('trust_source') else 'excluded_untrusted'))
            return dict(references=len(refs), files=checked, scientific_panel_suitability='not_assessed')
        attempt('reference_panel', check_references)
    else:
        add('reference_panel', 'conditional_not_provided', 'Empty default catalog; FastANI taxonomy unresolved')
    if args.planned_outdir:
        def check_output():
            out = fresh_result_dir(args.planned_outdir)
            if out == args.report_dir:
                raise ValueError('Planned analysis output must differ from preflight report')
            runtime.validate_output_inputs(out, protected)
            return dict(path=str(out), created=False)
        attempt('planned_output', check_output)
    else:
        add('planned_output', 'not_provided', 'No future analysis destination selected')
    attempt('report_isolation', lambda: runtime.validate_output_inputs(args.report_dir, protected))
    for tool in ('BACPHLIP', 'VIBRANT', 'vConTACT3'):
        add(tool, 'optional_not_enabled', 'Prepared opt-in V1 evidence; actual tool/workflow validation pending')
    add('CoverM', 'optional_not_enabled', 'Installed comparison tool; B1 fragment metrics remain independent')
    for tool in PAUSED:
        add(tool, 'paused', 'User-paused capability; excluded from readiness blockers')
    add('PHROGs_equivalence', 'scientific_limitation', 'Historical cross-platform strict differences unresolved; no rerun')
    add('real_validation', 'not_executed', 'No biological analysis/calibration or Nextflow tasks performed')
    report['blockers'] = sum(x['state']=='blocker' for x in report['checks'])
    report['status'] = 'blocked' if report['blockers'] else ('input_contract_ready_real_validation_pending' if inputs else 'environment_metadata_ready_input_not_provided')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir', type=Path, required=True)
    for name in ('input', 'metadata', 'genbanks', 'reference-dir', 'planned-outdir'):
        parser.add_argument('--'+name, type=Path)
    parser.add_argument('--host-association-config',type=Path)
    parser.add_argument('--host-association-bundle',type=Path,default=REPO/'assets/host_source_empty')
    parser.add_argument('--catalog-policy',type=Path)
    parser.add_argument('--catalog-evidence',type=Path,default=REPO/'assets/catalog_evidence_empty')
    parser.add_argument('--enable-kofam',action='store_true')
    parser.add_argument('--kofam-db',type=Path,default=Path('/srv/CRC-PHIRE/databases/kofam/2026-08-02'))
    parser.add_argument('--host-config',type=Path)
    parser.add_argument('--host-recombination-config',type=Path)
    parser.add_argument('--whole-genome-search-config',type=Path)
    parser.add_argument('--host-mask-bundle',type=Path,default=REPO/'assets/host_mask_empty')
    parser.add_argument('--max-total-bases' , type=int, default=20_000_000)
    args = parser.parse_args(argv)
    try:
        if not 1 <= args.max_total_bases <= 20_000_000:
            raise ValueError('Inspection total must be 1..20,000,000 bases')
        args.report_dir = fresh_result_dir(args.report_dir)
        for key in ('input','metadata','genbanks','reference_dir','planned_outdir'):
            if getattr(args,key):
                setattr(args,key,getattr(args,key).expanduser().resolve())
        if not args.input and (args.metadata or args.genbanks):
            raise ValueError('Metadata/GenBank checks require --input')
        report = inspect(args)
        # Isolation failure must not write even a report in a protected location.
        if any(x['check']=='report_isolation' and x['state']=='blocker' for x in report['checks']):
            raise ValueError('Report location overlaps protected inputs')
        args.report_dir.mkdir(exist_ok=False)
        (args.report_dir/'preflight.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
        lines = ['# Core readiness: '+report['status'], '', 'Read-only inspection; no analysis or calibration.', '']
        lines += [f"- {x['check']}: {x['state']} — {x['message']}" for x in report['checks']]
        (args.report_dir/'REPORT.md').write_text('\n'.join(lines)+'\n')
        print(json.dumps(dict(status=report['status'],blockers=report['blockers'],report=str(args.report_dir))))
        return 2 if report['blockers'] else 0
    except (ValueError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
