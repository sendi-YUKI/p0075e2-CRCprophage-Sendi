#!/usr/bin/env python3
"""Small deterministic interfaces shared by opt-in mainline extensions."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def read_table(path, required=()):
    with Path(path).open(newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        if not set(required).issubset(r.fieldnames or []):
            raise ValueError('Missing columns: ' + str(path))
        rows = list(r)
    if any(None in r or None in r.values() for r in rows):
        raise ValueError('Malformed table: ' + str(path))
    return rows


def write_table(path, rows, fields):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def unique(rows, key):
    result = {r[key]: r for r in rows}
    if len(result) != len(rows) or '' in result:
        raise ValueError('Duplicate/empty identity: ' + key)
    return result


def fasta(path):
    result = {}
    name = None
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError('Duplicate FASTA identity ' + name)
                result[name] = ''
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence without header')
                result[name] += line.strip().upper()
    if any(not seq for seq in result.values()):
        raise ValueError('Empty sequence')
    return result


def write_fasta(path, seqs):
    with Path(path).open('w') as f:
        for name, seq in seqs.items():
            f.write('>' + name + '\n' + seq + '\n')


def safe_id(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError('Unsafe identity: ' + value)
    return value


def run(argv, directory, label, stdout=None, timeout=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    env['PYTHONNOUSERSITE'] = '1'
    code = None
    try:
        with (directory / (label + '.stderr.log')).open('w') as err, Path(stdout or directory / (label + '.stdout.log')).open('w') as out:
            command = list(map(str, argv))
            if timeout is not None:
                # A separate timeout supervisor survives an enclosing helper being interrupted.
                command = ['/usr/bin/timeout', '--signal=TERM', '--kill-after=1s', str(timeout) + 's', *command]
            process = subprocess.Popen(command, stdout=out, stderr=err, env=env, start_new_session=True)
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                code = 124
    except subprocess.TimeoutExpired:
        code = 124
    finally:
        save(directory / (label + '.command.json'), dict(argv=list(map(str, argv)), exit_code=code,
             elapsed_seconds=round(time.monotonic() - start, 3), timeout_seconds=timeout))
    if code:
        raise RuntimeError(label + ' failed, exit ' + str(code))


def new_output(path):
    out = Path(path)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output must be new/empty: ' + str(out))
    out.mkdir(parents=True, exist_ok=True)
    return out
