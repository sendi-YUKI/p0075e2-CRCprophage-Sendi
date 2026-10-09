#!/usr/bin/env python3
"""Stable native task cache inputs: small sources/manifests/models, never large DB scans."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def stable_database(item):
    # Deliberately exclude job/output/log/receipt timestamps and mutable execution evidence.
    result = {key: item[key] for key in ['version', 'inventory_sha256', 'manifest_sha256',
              'fasta_sha256', 'amrfinder_source_lock_sha256'] if key in item}
    if 'files' in item:
        result['files'] = sorted([{key: row[key] for key in ['path', 'bytes', 'sha256'] if key in row}
                                  for row in item['files']], key=lambda row: row['path'])
    if 'symlinks' in item:
        result['symlinks'] = sorted(item['symlinks'], key=lambda row: row['path'])
    if item.get('index_receipt'):
        index = json.loads(Path(item['index_receipt']).read_text())
        result['native_index'] = {key: index[key] for key in
            ['sha256', 'bytes', 'source_sha256', 'diamond_binary_sha256',
             'diamond_version', 'container', 'platform'] if key in index}
    if item.get('runtime_manifest'):
        manifest = json.loads(Path(item['runtime_manifest']).read_text())
        result['package_models'] = {key: manifest[key] for key in
                                   ['environment_lock_sha256', 'installed_package_files'] if key in manifest}
    return result


def scientific_sources(root):
    paths = set(root.glob('main*.nf')) | {root/'nextflow.config'}
    for name in ['bin', 'workflows', 'modules', 'conf']:
        paths.update(path for path in (root/name).rglob('*') if path.is_file()
                     and '__pycache__' not in path.parts
                     and (path.suffix in ('.py', '.R', '.sh', '.nf', '.config', '.yml', '.yaml')
                          or (name == 'bin' and path.suffix == '')))
    for name in ['profiles', 'schemas']:
        paths.update(path for path in (root/'assets'/name).rglob('*') if path.is_file())
    for pattern in ['*tools.json', 'required_assets.json', '*policy*.json']:
        paths.update((root/'assets').glob(pattern))
    return {path.relative_to(root).as_posix(): sha(path) for path in sorted(paths)}


def genomad_models(environment):
    prefix = Path(environment['prefix']).resolve()
    roots = {path.resolve() for path in prefix.glob('lib/python*/site-packages/genomad/data')}
    if len(roots) != 1: raise ValueError('Cannot uniquely resolve pinned geNomad model directory')
    folder = roots.pop()
    if not folder.is_relative_to(prefix): raise ValueError('geNomad models escape pinned prefix')
    required = {'decision_forest.ubj', 'nn_classifier.h5', 'provirus_tagger.crfsuite',
                'rbs_categories.tsv', 'score_calibration_weights.npz'}
    files = sorted(path for path in folder.rglob('*') if path.is_file() and '__pycache__' not in path.parts)
    if not required.issubset({path.name for path in files}): raise ValueError('Missing geNomad package model')
    return {'environment_lock_sha256': environment['lock_sha256'],
            'required_assets': sorted(required),
            'files': [{'path': path.relative_to(folder).as_posix(), 'bytes': path.stat().st_size,
                       'sha256': sha(path)} for path in files]}


def phageboost_models(environment,contract):
    prefix=Path(environment['prefix']).resolve()
    roots={path.resolve() for path in prefix.glob('lib/python*/site-packages/PhageBoost/models')}
    if len(roots)!=1:raise ValueError('Cannot uniquely resolve pinned PhageBoost models')
    folder=roots.pop()
    if not folder.is_relative_to(prefix):raise ValueError('PhageBoost models escape pinned prefix')
    rows=[]
    for expected in sorted(contract['models'],key=lambda row:row['relative_path']):
        name=expected['relative_path']
        if Path(name).name!=name:raise ValueError('Unsafe PhageBoost model filename')
        path=folder/name
        if not path.resolve().is_relative_to(folder):raise ValueError('PhageBoost model escapes package')
        if not path.is_file() or path.stat().st_size!=expected['bytes'] or sha(path)!=expected['sha256']:
            raise ValueError('Missing or changed pinned PhageBoost model: '+name+'; run bootstrap/spark/restore_phageboost_model.py inside the approved allocation')
        rows.append(dict(relative_path=name,bytes=expected['bytes'],sha256=expected['sha256']))
    return dict(source_commit=contract['git_commit'],environment_lock_sha256=environment['lock_sha256'],files=rows)


def build(root=ROOT):
    dependencies = json.loads((root/'assets/spark_dependencies.json').read_text())
    databases = json.loads((root/'assets/spark_databases.json').read_text())
    material = {'schema_version': 1, 'scientific_sources': scientific_sources(root),
                'databases': {key: stable_database(item) for key, item in sorted(databases['databases'].items())},
                'genomad_models': genomad_models(dependencies['environments']['core-genomad']),
                'phageboost_models': phageboost_models(dependencies['environments']['core-callers'],
                    json.loads((root/'assets/core_callers_tools.json').read_text())['tools']['PhageBoost'])}
    return {'schema_version': 1, 'fingerprint': digest(material), 'material': material}


if __name__ == '__main__': print(json.dumps(build(), sort_keys=True))
