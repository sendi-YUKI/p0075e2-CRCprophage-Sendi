#!/usr/bin/env python3
"""Restore the explicitly reviewed AMRFinder database; never choose a latest release."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import urllib.request

VERSION = '2026-08-07.1'
OFFICIAL = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/' + VERSION + '/'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def validate_lock(lock):
    if lock.get('version') != VERSION or lock.get('official_base_url') != OFFICIAL:
        raise ValueError('Only the reviewed AMRFinder 2026-08-07.1 official release is allowed')
    names = set()
    for row in lock['source_files']:
        name = row['name']
        if Path(name).name != name or name in {'.', '..'} or name in names:
            raise ValueError('Unsafe or duplicate AMRFinder source filename')
        names.add(name)
        if row['url'] != OFFICIAL + name or not re.fullmatch('[a-f0-9]{64}', row['sha256']) or not isinstance(row['bytes'], int) or row['bytes'] < 0:
            raise ValueError('Invalid pinned AMRFinder file identity')
    required = {'AMR.LIB', 'AMRProt.fa', 'AMRProt-mutation.tsv', 'AMRProt-suppress.tsv',
                'AMRProt-susceptible.fa', 'AMRProt-susceptible.tsv', 'AMR_CDS.fa',
                'database_format_version.txt', 'fam.tsv', 'taxgroup.tsv', 'version.txt', 'changes.txt'}
    if not required.issubset(names):
        raise ValueError('Incomplete official AMRFinder source set')


def restore(database, lock):
    validate_lock(lock)
    database = Path(database).resolve()
    database.mkdir(parents=True, exist_ok=True)
    target = database / VERSION
    if target.is_symlink():
        raise ValueError('Version directory must not be a symlink')
    target.mkdir(exist_ok=True)
    for row in lock['source_files']:
        path = target / row['name']
        if path.is_symlink():
            raise ValueError('Source file must not be a symlink')
        if path.exists():
            if path.stat().st_size != row['bytes'] or sha(path) != row['sha256']:
                raise ValueError('Existing pinned AMRFinder source changed: ' + row['name'])
            continue
        temporary = path.with_name(path.name + '.pinned.part')
        if temporary.is_symlink():
            raise ValueError('Download temporary file must not be a symlink')
        with urllib.request.urlopen(row['url'], timeout=120) as response, temporary.open('wb') as output:
            for block in iter(lambda: response.read(1024 * 1024), b''):
                output.write(block)
        if temporary.stat().st_size != row['bytes'] or sha(temporary) != row['sha256']:
            raise ValueError('Official AMRFinder source checksum mismatch: ' + row['name'])
        temporary.replace(path)
    if (target / 'version.txt').read_text().strip() != VERSION:
        raise ValueError('Downloaded AMRFinder version marker disagrees')
    # Same index builder used by NCBI amrfinder_update; caller fixes the container.
    subprocess.run(['amrfinder_index', str(target)], check=True)
    latest = database / 'latest'
    if latest.exists() and not latest.is_symlink():
        raise ValueError('Refusing to replace a non-symlink latest directory')
    pending = database / 'latest.pinned.tmp'
    if pending.is_symlink():
        pending.unlink()
    if pending.exists():
        raise ValueError('Unexpected file at temporary latest link')
    pending.symlink_to(VERSION)
    pending.replace(latest)
    return dict(status='pinned_sources_verified_and_indexed', version=VERSION,
                official_base_url=OFFICIAL, source_files=len(lock['source_files']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', required=True)
    parser.add_argument('--lock', required=True)
    args = parser.parse_args()
    print(json.dumps(restore(args.database, json.loads(Path(args.lock).read_text()))))


if __name__ == '__main__':
    main()
