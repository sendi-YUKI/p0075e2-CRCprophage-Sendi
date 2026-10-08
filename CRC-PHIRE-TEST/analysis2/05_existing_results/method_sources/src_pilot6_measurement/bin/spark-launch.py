#!/usr/bin/env python3
"""Submit or execute one output-bound pipeline inside a native ARM allocation."""
import argparse
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import uuid

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spark_runtime as runtime

REPO = runtime.REPO
ENTRIES = {'project':'run-project','a0': 'run', 'core': 'run-core', 'core-from-reads': 'run-core-from-reads',
           'fetch-genomes': 'spark-fetch-genomes',
           'reads': 'run-stage', 'assembly': 'run-stage', 'viral': 'run-stage',
           'cohort': 'run-stage', 'association': 'run-stage'}
PATH_OPTIONS = {'--conjscan-models-dir','--element-pfam-db','--config','--host-recombination-config','--whole-genome-search-config','--research-catalog','--research-source-manifest','--research-reference-config','--host-association-config', '--host-association-bundle', '--catalog-policy', '--catalog-evidence', '--host-config', '--host-mask-bundle', '--kofam-db', '--outdir', '--input', '--a1-result', '--database-root', '--work-dir',
                '--checkm2-db', '--functions-db', '--reference-dir', '--metadata', '--genbanks',
                '--assembly-result', '--samplesheet', '--merge-map', '--clean-reads-manifest',
                '--assembly-manifest', '--catalog-dir', '--legacy-catalog', '--legacy-source-manifest',
                '--canonical-gene-tables', '--preprocessing-profile', '--measurement-profile',
                '--detection-profiles', '--host-references', '--decoy-fasta', '--association-long',
                '--association-metadata', '--association-spec'}


def normalize_paths(args):
    result = []
    option = None
    for value in args:
        if value.startswith('--'):
            option = value.split('=', 1)[0]
            if '=' in value and option in PATH_OPTIONS:
                result.append(option + '=' + str(Path(value.split('=', 1)[1]).expanduser().resolve()))
                option = None
            else:
                result.append(value)
        elif option in PATH_OPTIONS:
            result.append(str(Path(value).expanduser().resolve()))
            if option != '--canonical-gene-tables':
                option = None
        else:
            result.append(value)
    return result


def declared_inputs(args):
    inputs = []
    option = None
    for value in args:
        if value.startswith('--'):
            option = value.split('=', 1)[0]
            if '=' in value and option in PATH_OPTIONS and option != '--outdir':
                inputs.append(value.split('=', 1)[1]); option = None
        elif option in PATH_OPTIONS and option != '--outdir':
            inputs.append(value)
            if option != '--canonical-gene-tables':
                option = None
    return inputs


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def output_path(args):
    values = [args[n + 1] for n, value in enumerate(args[:-1]) if value == '--outdir']
    values += [value.split('=', 1)[1] for value in args if value.startswith('--outdir=')]
    if len(values) != 1:
        raise ValueError('Pass exactly one explicit --outdir to the pipeline')
    path = Path(values[0]).expanduser().resolve()
    if any(x in str(path) for x in ('\n', '\r', "'")):
        raise ValueError('Output paths containing quotes/newlines are unsupported')
    campaign = Path('/srv/CRC-PHIRE/analysis1')
    if not path.is_relative_to(campaign):
        raise ValueError('analysis1 runtime outputs must stay inside its campaign')
    for protected in (Path.home() / 'databases', campaign / 'work', campaign / 'spark_preflight', REPO, runtime.DEPENDENCY_REPO / 'envs', REPO / 'baseline_evidence'):
        if path == protected or path.is_relative_to(protected) or protected.is_relative_to(path):
            raise ValueError('Output overlaps protected database/environment/work/evidence path')
    if path == REPO or REPO.is_relative_to(path):
        raise ValueError('Output must not contain the repository')
    runtime.validate_output_inputs(path, declared_inputs(args))
    return path


def command(entry, args):
    prefix = [str(REPO / 'bin' / ENTRIES[entry])]
    if ENTRIES[entry] == 'run-stage':
        prefix.append(entry)
    return prefix + args


def execute(entry, args, record):
    out = output_path(args)
    receipt_path = REPO / 'bootstrap/spark/executions' / (record + '.json')
    if receipt_path.exists():
        record += '-' + uuid.uuid4().hex[:12]
        receipt_path = receipt_path.with_name(record + '.json')
    state = {'schema_version': 1, 'target': 'SPARK_PORTABILITY_BUILD', 'status': 'waiting_for_resource_lock',
             'execution_id': record,
             'entry': entry, 'args': args, 'outdir': str(out), 'started_utc': now(),
             'runtime': 'native-locked-prefix',
             'scientifically_calibrated': False, 'real_data_validated': False,
             'resume_requested': '--resume' in args}
    runtime.save(receipt_path, state)
    lock_path = Path('/srv/CRC-PHIRE/analysis1/core-heavy.lock')
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    child = None
    def interrupted(signum, frame):
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signum)
        raise RuntimeError('Execution interrupted by signal ' + str(signum))
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        state['allocation'] = runtime.allocation(require_lock=False)
        with lock_path.open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            runtime.write_runtime_identities()
            env = dict(os.environ)
            env.update(CRC_PHAGE_RUNTIME='spark-native',
                       CRC_SPARK_HEAVY_LOCK_OWNER=str(os.getpid()),
                       CRC_SPARK_HEAVY_LOCK_FD=str(lock.fileno()),
                       CRC_SPARK_PREFLIGHT_ROOT=str(Path('/srv/CRC-PHIRE/analysis1/spark_preflight') / record),
                       CRC_SPARK_ALLOWED_WRITE_ROOTS='[]',
                       NXF_OPTS='-Xms256m -Xmx2g -XX:ActiveProcessorCount=' + str(min(6, state['allocation']['cpus'])))
            # Exact native config and prefix lock fingerprints are recorded even
            # when the pipeline later fails preflight.
            state.update(status='running', command=command(entry, args),
                         native_config_sha256=runtime.sha(REPO / 'conf/spark_native.config'),
                         dependencies_sha256=runtime.sha(runtime.DEPENDENCIES))
            runtime.save(receipt_path, state)
            child = subprocess.Popen(state['command'], cwd=REPO, env=env, start_new_session=True)
            result = child.wait()
            state.update(exit_code=result, status='completed' if result == 0 else 'failed')
            # Read UUID only from the relevant output-bound production receipt.
            candidate = out / 'pipeline_info' / ('acquisition_state.json' if entry == 'fetch-genomes' else ('stage_state.json' if ENTRIES[entry] == 'run-stage' else 'run_status.json'))
            if candidate.is_file():
                pipeline = runtime.read(candidate)
                if result and pipeline.get('status') in ('running', 'preflight', 'validating_outputs'):
                    pipeline.update(status='failed', spark_child_exit_code=result,
                                    error='Child exited without closing its running state; see spark execution receipt')
                    runtime.save(candidate, pipeline)
                state['pipeline_state_path'] = str(candidate)
                state['pipeline_state_sha256'] = runtime.sha(candidate)
                if entry == 'fetch-genomes':
                    state['nextflow_workflow'] = False
                    state['execution_scope'] = 'explicit_versioned_genome_acquisition_not_analysis'
                    if result == 0:
                        from spark_acquire import verify_completed_acquisition
                        verified = verify_completed_acquisition(out)
                        state['output_sha256'] = verified['output_sha256']
                else:
                    state['nextflow_session_id'] = pipeline.get('nextflow_session_id')
                if entry != 'fetch-genomes' and result == 0 and not re.fullmatch(r'[a-f0-9]{8}(?:-[a-f0-9]{4}){3}-[a-f0-9]{12}', state.get('nextflow_session_id') or ''):
                    raise RuntimeError('Successful pipeline lacks output-bound Nextflow UUID')
            elif result == 0:
                raise RuntimeError('Successful command lacks production pipeline state')
    except Exception as error:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL); child.wait()
        state.update(status='failed', exit_code=1, error=str(error))
    finally:
        state['finished_utc'] = now()
        runtime.save(receipt_path, state)
        if (out / 'pipeline_info').is_dir():
            runtime.save(out / 'pipeline_info/spark_execution.json', state)
    print(json.dumps(state, indent=2))
    return state['exit_code']


def submit(entry, args, time_limit):
    if os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Already inside Slurm; use bin/spark-run, not a nested submission')
    site = runtime.read(runtime.SITE)
    if not site.get('scheduler', {}).get('actual_smoke_verified'):
        raise ValueError('Site execution smoke is not verified')
    out = output_path(args)
    # Validate walltime without permitting additional sbatch options.
    if not re.fullmatch(r'(?:[0-9]{1,2}-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}', time_limit):
        raise ValueError('Use Slurm time syntax HH:MM:SS or D-HH:MM:SS')
    record = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:12]
    directory = REPO / 'bootstrap/spark/submissions' / record
    directory.mkdir(parents=True, exist_ok=False)
    script = directory / 'job.sh'
    script.write_text('#!/bin/bash\nset -euo pipefail\ncd ' + shlex.quote(str(REPO)) + '\nexec ' +
                      shlex.join([str(REPO / 'bin/spark-run'), '--record-id', record, entry, '--', *args]) + '\n')
    cmd = ['sbatch', '--parsable', '--account=lab', '--partition=main', '--nodes=1', '--ntasks=1',
           '--cpus-per-task=6', '--mem=48G', '--time=' + time_limit, '--dependency=singleton',
           '--job-name=crc-phage-spark', '--output=' + str(directory / 'slurm-%j.out'),
           '--error=' + str(directory / 'slurm-%j.err'), str(script)]
    receipt = {'schema_version': 1, 'target': 'SPARK_PORTABILITY_BUILD', 'status': 'submitting',
               'entry': entry, 'args': args, 'outdir': str(out), 'command': cmd,
               'created_utc': now(), 'execution_receipt': str(REPO / 'bootstrap/spark/executions' / (record + '.json')),
               'resource_limits': {'cpus': 6, 'memory_gib': 48, 'max_parallel_heavy': 1},
               'job_script_sha256': runtime.sha(script)}
    path = directory / 'submission.json'
    runtime.save(path, receipt)
    result = subprocess.run(cmd, capture_output=True, text=True)
    receipt.update(submission_exit_code=result.returncode, stderr=result.stderr, stdout=result.stdout)
    job = result.stdout.strip().split(';')[0]
    receipt.update(status='submitted' if result.returncode == 0 and job.isdigit() else 'submission_failed', job_id=job if job.isdigit() else None)
    runtime.save(path, receipt)
    print(json.dumps(receipt, indent=2))
    return 0 if receipt['status'] == 'submitted' else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['submit', 'run'])
    p.add_argument('--time', default='24:00:00')
    p.add_argument('--record-id')
    p.add_argument('entry', choices=sorted(ENTRIES))
    p.add_argument('args', nargs=argparse.REMAINDER)
    a = p.parse_args()
    args = normalize_paths(a.args[1:] if a.args[:1] == ['--'] else a.args)
    if a.action == 'submit':
        return submit(a.entry, args, a.time)
    record = a.record_id or datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:12]
    if not re.fullmatch(r'[A-Za-z0-9_-]+', record):
        p.error('Invalid record ID')
    return execute(a.entry, args, record)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)
