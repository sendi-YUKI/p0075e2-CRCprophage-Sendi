#!/usr/bin/env python3
"""Explicit native ARM runtime; Docker behavior remains the default elsewhere.

This is a locked prefix runtime, not a container. It does not claim filesystem
isolation or read-only mounts. Slurm cgroups enforce the allocation ceiling.
"""
import argparse
import contextlib
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
DEPENDENCY_REPO = Path('/srv/CRC-PHIRE/projects/crc-pipeline')
DEPENDENCIES = REPO / 'assets/spark_dependencies.json'
SITE = REPO / 'bootstrap/spark/site_execution.json'
IDENTITIES = REPO / 'assets/spark_runtime_identities.json'
PROCESS_ENV = {
    'HOST_RECOMBINATION':'recombination',
    'CORE_WHOLE_GENOME_EVIDENCE':'core-blast',
    'VIRAL_KOFAM_DB_VALIDATE':'kofam',
    'VIRAL_KOFAM':'kofam',
    'RESEARCH_V1_BRIDGE': 'core-blast',
    'HOST_SOURCE_ASSOCIATION': 'python',
    'CORE_RESEARCH_CATALOG': 'python',
    'CORE_LOCUS_EVIDENCE': 'core-blast',
    'CORE_RESEARCH_SUMMARY': 'python',
    'CORE_KOFAM_DB_VALIDATE': 'kofam',
    'CORE_KOFAM': 'kofam',
    'CORE_KOFAM_SUMMARY': 'python',
    'HOST_CORE_GFF': 'python',
    'HOST_CORE_PANAROO': 'host-core',
    'HOST_CORE_ALIGN': 'host-tree',
    'HOST_TREE_PRELIMINARY': 'host-tree',
    'HOST_MOBILE_MASK': 'python',
    'HOST_TREE_FINAL': 'host-tree',

    'GENOMAD_ENDTOEND': 'core-genomad', 'VIRAL_GENOMAD': 'core-genomad',
    'CHECKV_ENDTOEND': 'core-checkv', 'VIRAL_CHECKV': 'core-checkv',
    'READS_CLEAN': 'reads-qc', 'READS_COLLECT': 'reads-qc',
    'READS_ASSEMBLE': 'reads-assembly', 'READS_ANNOTATE': 'reads-bakta',
    'COHORT_INDEX': 'reads-qc', 'COHORT_MAPPABILITY': 'reads-qc',
    'COHORT_MEASURE': 'reads-qc', 'COHORT_METAPHLAN': 'reads-metaphlan',
    'ASSOCIATION_MODELS': 'association', 'CORE_GENES': 'core-checkm2',
    'CORE_CHECKM2': 'core-checkm2', 'CORE_FASTANI': 'core-fastani',
    'CORE_BLAST': 'core-blast', 'VIRAL_BLAST': 'core-blast',
    'CORE_FUNCTIONS': 'core-functions', 'CORE_FUNCTION_SUMMARY': 'core-functions',
    'VIRAL_FUNCTIONS': 'core-functions', 'VIRAL_FUNCTION_SUMMARY': 'core-functions',
    'CORE_PHISPY': 'core-callers', 'CORE_PHAGEBOOST': 'core-callers',
    'CORE_PHISPY_UNASSESSED': 'core-callers', 'CORE_COMPARE_CALLERS': 'core-callers',
    'VIRAL_BACPHLIP': 'bacphlip',
    'VIRAL_VIBRANT': 'vibrant',
    'VIRAL_VCONTACT3': 'vcontact3',
    'VIRAL_EVIDENCE_SUMMARY': 'python',
    'VIRAL_GENES': 'viral-gene',
}


def enabled():
    return os.environ.get('CRC_PHAGE_RUNTIME') == 'spark-native'


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.' + str(os.getpid()) + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def allocation(require_lock=True):
    """Reject login-shell execution and allocations beyond the authorized cap."""
    if platform.system() != 'Linux' or platform.machine() not in ('aarch64', 'arm64'):
        raise RuntimeError('Spark native execution requires Linux aarch64')
    job = os.environ.get('SLURM_JOB_ID', '')
    if not job.isdigit():
        raise RuntimeError('Use bin/spark-submit; heavy work must run inside Slurm')
    # Environment values alone are insufficient: verify the actual live allocation.
    output = subprocess.check_output(['scontrol', 'show', 'job', '-o', job], text=True)
    fields = dict(re.findall(r'(\w+)=([^\s]+)', output))
    if fields.get('JobState') not in ('RUNNING', 'COMPLETING'):
        raise RuntimeError('Slurm allocation is not running')
    if not fields.get('UserId', '').endswith('(' + str(os.getuid()) + ')'):
        raise RuntimeError('Slurm allocation belongs to another UID')
    cpus = int(fields.get('NumCPUs', '0'))
    tres = fields.get('AllocTRES', '')
    mem_match = re.search(r'(?:^|,)mem=([0-9]+)([KMGTP]?)', tres)
    if not mem_match:
        raise RuntimeError('Cannot verify Slurm allocated memory')
    memory = int(mem_match[1]) * {'': 1, 'K': 1 / 1024, 'M': 1, 'G': 1024, 'T': 1024 ** 2, 'P': 1024 ** 3}[mem_match[2]]
    profile=os.environ.get('CRC_SPARK_RESOURCE_PROFILE','standard')
    if profile not in ('standard','exclusive20'):raise RuntimeError('Unknown Spark resource profile')
    cap_cpu,cap_mem=(20,110*1024) if profile=='exclusive20' else (6,49152)
    if not 1 <= cpus <= cap_cpu or not 1 <= memory <= cap_mem:
        raise RuntimeError('Allocation exceeds selected Spark resource profile')
    if fields.get('Partition') != 'main' or fields.get('Account') != 'lab':
        raise RuntimeError('Unreviewed Slurm account/partition')
    if 'gres/gpu' in tres or os.environ.get('SLURM_JOB_GPUS'):
        raise RuntimeError('This CPU pipeline has not requested a GPU')
    cgroup = Path('/proc/self/cgroup').read_text().strip()
    if not re.search(r'(?:^|/)job_' + re.escape(job) + r'(?:/|$)', cgroup):
        raise RuntimeError('Current process is outside the verified Slurm job cgroup')
    if require_lock:
        owner = os.environ.get('CRC_SPARK_HEAVY_LOCK_OWNER', '')
        fd = os.environ.get('CRC_SPARK_HEAVY_LOCK_FD', '')
        lock = Path(os.environ.get('CRC_PHAGE_WORK_ROOT', str(Path.home() / 'scratch/crc-phage'))) / 'core-heavy.lock'
        if not owner.isdigit() or not fd.isdigit():
            raise RuntimeError('Missing parent resource lock; use bin/spark-run inside allocation')
        try:
            observed = Path('/proc') / owner / 'fd' / fd
            try:
                matched=observed.resolve(strict=True)==lock.resolve()
            except PermissionError:
                # Landlock/no_new_privs may deny another process's /proc fd link.
                # Accept only this task-shell's sealed inherited receipt, never a JSON path/env bypass.
                attested=os.environ.get('CRC_SPARK_LOCK_ATTEST_FD','')
                if not attested.isdigit():raise RuntimeError('Missing sealed task lock receipt')
                afd=int(attested);seals=fcntl.fcntl(afd,getattr(fcntl,"F_GET_SEALS",1034))
                required=15 # Linux WRITE|SHRINK|GROW|SEAL
                if seals & required != required:raise RuntimeError('Unsealed task lock receipt')
                if 'memfd:crc-spark-verified-parent-lock' not in os.readlink('/proc/self/fd/'+attested):
                    raise RuntimeError('Unexpected lock receipt descriptor')
                payload=json.loads(os.pread(afd,4096,0))
                st=lock.stat()
                matched=payload==dict(job_id=job,uid=os.getuid(),lock_device=st.st_dev,lock_inode=st.st_ino,parent_pid=owner)
                if not os.environ.get('CRC_SPARK_WRITE_GUARD','').startswith('landlock-abi'):
                    raise RuntimeError('Lock receipt fallback only inside guarded task')
            if not matched:
                raise RuntimeError('Parent lock does not identify the shared heavy lock')
            with lock.open('r') as probe:
                try:
                    fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    pass
                else:
                    fcntl.flock(probe, fcntl.LOCK_UN)
                    raise RuntimeError('Parent heavy lock is not held')
        except OSError as error:
            raise RuntimeError('Parent resource lock disappeared') from error
    return {'job_id': job, 'cpus': cpus, 'memory_mib': memory,
            'partition': 'main', 'account': 'lab', 'alloc_tres': tres,
            'cgroup': cgroup}


def environment_definition(key):
    """Read locked prefix metadata only; no allocation, execution or binary scan."""
    data = read(DEPENDENCIES)
    if data.get('schema_version') != 1 or data.get('platform') != 'linux-aarch64':
        raise RuntimeError('Dependencies must be a linux-aarch64 schema_version=1 manifest')
    item = data.get('environments', {}).get(key)
    if not item or item.get('status') != 'ready':
        raise RuntimeError('Native ARM environment is not ready: ' + key)
    prefix = Path(item['prefix']).expanduser().resolve()
    manifest = Path(item['manifest']).expanduser().resolve()
    if not prefix.is_relative_to(DEPENDENCY_REPO / 'envs') or not manifest.is_relative_to(DEPENDENCY_REPO / 'envs'):
        raise RuntimeError('Native prefixes/manifests must be inside this repository envs/')
    lock = Path(item.get('lock', str(manifest.parent / 'pixi.lock'))).resolve()
    if not prefix.is_dir() or not lock.is_file() or sha(lock) != item.get('lock_sha256'):
        raise RuntimeError('Missing prefix or changed lock for ' + key)
    if item.get('manifest_sha256') and (not manifest.is_file() or sha(manifest) != item['manifest_sha256']):
        raise RuntimeError('Native environment manifest changed after installation: ' + key)
    if not (prefix / 'bin/python').is_file() and not (prefix / 'bin/python3').is_file():
        raise RuntimeError('Environment requires its own Python: ' + key)
    activation = dict(item.get('activation_env', {}))
    if any(not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', k) or not isinstance(v, str)
           for k, v in activation.items()):
        raise RuntimeError('Invalid environment activation variables: ' + key)
    return data, item, prefix, lock, activation


def environment(key):
    data, item, prefix, lock, activation = environment_definition(key)
    identity = stable_identity(key, item, data)
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    item = dict(item, runtime_identity=identity, runtime_fingerprint=fingerprint)
    activation.update(PATH=str(prefix / 'bin') + ':' + os.environ.get('PATH', ''),
                      CONDA_PREFIX=str(prefix), PYTHONNOUSERSITE='1',
                      CRC_SPARK_REPO=str(REPO),
                      CRC_SPARK_ENVIRONMENT=key, CRC_SPARK_ENV_LOCK_SHA256=sha(lock),
                      CRC_SPARK_ENV_FINGERPRINT=fingerprint)
    return item, activation


def stable_identity(key, item, registry):
    """Exclude probe timestamps/log paths but include all non-Pixi binaries."""
    normalize = lambda value: value.replace(str(REPO), '${REPO}') if isinstance(value, str) else value
    identity = {'schema_version': 1, 'platform': 'linux-aarch64', 'environment': key,
                'manifest_sha256': item.get('manifest_sha256') or sha(item['manifest']),
                'lock_sha256': item['lock_sha256'],
                'activation_env': {k: normalize(v) for k, v in item.get('activation_env', {}).items()},
                'tools': {name: {k: tool.get(k) for k in ('expected_version', 'observed_version')}
                          for name, tool in item.get('tools', {}).items()}}
    if item.get('source_inventory'):
        source_lock=Path(item['source_packages_lock'])
        inventory=Path(item['source_inventory'])
        if sha(source_lock)!=item['source_packages_lock_sha256'] or sha(inventory)!=item['source_inventory_sha256']:
            raise RuntimeError('Changed source-package lock/inventory: '+key)
        source_root=Path(item['manifest']).parent.resolve()
        for row in read(inventory):
            path=(source_root/row['path']).resolve()
            if not path.is_relative_to(source_root) or not path.is_file() or sha(path)!=row['sha256']:
                raise RuntimeError('Changed installed source package: '+row['path'])
        identity.update(source_packages_lock_sha256=item['source_packages_lock_sha256'],source_inventory_sha256=item['source_inventory_sha256'])
    source = item.get('native_source')
    if source:
        identity['native_source'] = {k: source[k] for k in
            ('url', 'bytes', 'sha256', 'expected_md5', 'commit', 'architecture', 'binary_sha256') if k in source}
        binaries = source.get('binary_sha256', {})
        if isinstance(binaries, str):
            if key != 'core-checkv':
                raise RuntimeError('Unmapped native source binary: ' + key)
            binaries = {'diamond': binaries}
        if not binaries:
            raise RuntimeError('Native source has no executable checksums: ' + key)
        for name, expected in binaries.items():
            binary = Path(item['prefix']) / 'bin' / name
            if not binary.is_file() or sha(binary) != expected:
                raise RuntimeError('Native executable changed: ' + str(binary))
        if key == 'core-checkv':
            # The current record may point to a later version-probe-only attempt.
            # Bind commands to the historical build receipt for this exact binary.
            commands = []
            for receipt in sorted((DEPENDENCY_REPO / 'results/spark_portability/dependencies').glob('*/core-checkv/diamond_native_source.json')):
                built = read(receipt)
                if built.get('binary_sha256') != source['binary_sha256'] or built.get('sha256') != source['sha256']:
                    continue
                pair = []
                for name in ('diamond_configure.log', 'diamond_build.log'):
                    log = receipt.parent / name
                    if not log.is_file():
                        break
                    argv = json.loads(log.read_text().splitlines()[0])
                    if not isinstance(argv, list) or not all(isinstance(v, str) for v in argv):
                        raise RuntimeError('Malformed native build argv receipt')
                    pair.append([normalize(v) for v in argv])
                if len(pair) == 2 and pair not in commands:
                    commands.append(pair)
            if not commands:
                raise RuntimeError('DIAMOND native build commands missing for installed binary')
            identity['build_commands'] = sorted(commands, key=lambda value: json.dumps(value))
            identity['build_tools_lock_sha256'] = registry['environments']['build-tools']['lock_sha256']
        elif key == 'core-blast':
            identity['installation'] = 'official_aarch64_archive_extracted; prefix_bin_symlinks; no_recompile'
    return identity


def write_runtime_identities():
    registry = read(DEPENDENCIES)
    identities = {}
    for key, item in sorted(registry['environments'].items()):
        if item.get('status') == 'ready':
            actual, _ = environment(key)
            identities[key] = {'fingerprint': actual['runtime_fingerprint'], 'identity': actual['runtime_identity']}
    payload = json.dumps({'schema_version': 1, 'platform': 'linux-aarch64', 'environments': identities},
                         indent=2, sort_keys=True) + '\n'
    if not IDENTITIES.is_file() or IDENTITIES.read_text() != payload:
        temporary = IDENTITIES.with_name(IDENTITIES.name + '.tmp')
        temporary.write_text(payload); temporary.replace(IDENTITIES)
    return identities


def command(key, argv, writable=()):
    item, activation = environment(key)
    scratch = Path(os.environ.get('CRC_SPARK_PREFLIGHT_ROOT', ''))
    root = Path(os.environ.get('CRC_PHAGE_WORK_ROOT', str(Path.home() / 'scratch/crc-phage'))).resolve()
    if not scratch.is_absolute() or not scratch.resolve().is_relative_to(root / 'spark_preflight'):
        raise RuntimeError('Missing safe preflight scratch; use bin/spark-run')
    scratch.mkdir(parents=True, exist_ok=True)
    activation.update(HOME=str(scratch), TMPDIR=str(scratch), TMP=str(scratch), TEMP=str(scratch),
                      XDG_CACHE_HOME=str(scratch / 'cache'), PYTHONDONTWRITEBYTECODE='1')
    # /usr/bin/env receives argv elements directly, never shell-expanded strings.
    guard = ['/usr/bin/python3', str(REPO / 'bin/spark_fs_guard.py'), '--writable', str(scratch)]
    for path in writable:
        guard += ['--writable', str(Path(path).resolve(strict=True))]
    return guard + ['--', '/usr/bin/env', '-u', 'PYTHONPATH', '-u', 'PYTHONHOME',
            *[k + '=' + v for k, v in activation.items()], *map(str, argv)]


def validate_output_inputs(out, inputs):
    from spark_fs_guard import reject_ipc_storage
    inputs=list(inputs)
    reject_ipc_storage([out,*inputs])
    out = Path(out).expanduser().resolve()
    protected = [REPO / name for name in ('bin', 'modules', 'workflows', 'conf', 'assets', 'envs', 'baseline_evidence')]
    protected += [Path(value).expanduser().resolve() for value in inputs if value]
    for source in protected:
        if source == out or out in source.parents or source in out.parents:
            raise ValueError('Output overlaps protected input/code: ' + str(source))
    return True


def task_temp_alias(temp, alias_root=None):
    """Short persistent cache alias; data still resides inside the original task."""
    temp=Path(temp).resolve(strict=True)
    from spark_fs_guard import reject_ipc_storage
    reject_ipc_storage([temp])
    if not temp.is_dir():raise ValueError('Task temporary path must be a directory')
    root=Path(alias_root) if alias_root is not None else Path.home()/'scratch/crc-phage/t'
    root=root.absolute()
    try:root.mkdir(mode=0o700)
    except FileExistsError:pass
    info=root.lstat()
    if root.is_symlink() or not root.is_dir() or info.st_uid!=os.getuid() or info.st_mode & 0o077:
        raise ValueError('Private temporary alias directory has unknown ownership/permissions; refusing to change it')
    alias=root/hashlib.sha256(os.fsencode(str(temp))).hexdigest()[:12]
    if len(os.fsencode(alias))>70:raise ValueError('Short temporary alias still exceeds the AF_UNIX-safe path budget')
    created=False
    try:alias.symlink_to(temp,target_is_directory=True);created=True
    except FileExistsError:pass
    info=alias.lstat()
    if not alias.is_symlink() or info.st_uid!=os.getuid() or os.readlink(alias)!=str(temp) or alias.resolve(strict=True)!=temp:
        raise ValueError('Temporary alias conflicts with an existing object; nothing overwritten')
    return {'path':str(alias),'target':str(temp),'created':created,'uid':info.st_uid,
            'device':info.st_dev,'inode':info.st_ino,'scope':'persistent project work/cache alias; no additional Landlock grant'}


def image_environment(image):
    mapping = {}
    containers = read(REPO / 'assets/containers.json')
    for old, new in [('python', 'python'), ('genomad', 'core-genomad'), ('checkv', 'core-checkv')]:
        mapping[containers[old]] = new
    for key, item in read(REPO / 'assets/reads_tools.json').items():
        if isinstance(item, dict) and item.get('container'):
            mapping[item['container']] = key
    assoc = read(REPO / 'assets/association_tools.json')
    mapping[assoc.get('container', assoc.get('image_id'))] = 'association'
    # Preserve old image references as source provenance; no image is launched.
    for name, key in [('core_genomes_tools.json', 'core-checkm2'),
                      ('core_callers_tools.json', 'core-callers'),
                      ('core_functions_tools.json', 'core-functions')]:
        path = REPO / 'assets' / name
        if not path.exists():
            continue
        def walk(value, hint=key):
            if isinstance(value, dict):
                for k, v in value.items():
                    child = {'fastani': 'core-fastani', 'blast': 'core-blast',
                             'checkm2': 'core-checkm2'}.get(k.lower(), hint)
                    if k in ('container', 'image_id') and isinstance(v, str):
                        mapping[v] = hint
                    else:
                        walk(v, child)
            elif isinstance(value, list):
                for v in value:
                    walk(v, hint)
        walk(read(path))
    # BLAST is independently pinned in the workflow configuration.
    for filename in ['conf/viral.config', 'modules/local/core_votu.nf']:
        for value in re.findall(r"['\"]([^'\"]*blast[^'\"]*@sha256:[0-9a-f]+)['\"]", (REPO / filename).read_text()):
            mapping[value] = 'core-blast'
    if image not in mapping:
        raise RuntimeError('No reviewed native environment for legacy image ' + str(image))
    return mapping[image]


def preflight_image(image):
    allocation()
    key = image_environment(image)
    item, _ = environment(key)
    return {'runtime': 'native-locked-prefix', 'legacy_container_reference': image,
            'environment': key, 'prefix': item['prefix'], 'lock_sha256': item['lock_sha256'],
            'platform': 'linux-aarch64', 'container_executed': False,
            'readonly_mounts_enforced': False}


def nextflow_options():
    allocation()
    options=['-c', str(REPO / 'conf/spark_native.config')]
    if os.environ.get('CRC_SPARK_RESOURCE_PROFILE')=='exclusive20':options+=['-c',str(REPO/'conf/spark_exclusive.config')]
    return options


def native_record(info, keys):
    from spark_fs_guard import IPC_POLICY
    evidence = {'runtime': 'native-locked-prefix', 'platform': 'linux-aarch64',
                'allocation': allocation(), 'site_manifest_sha256': sha(SITE),
                'dependency_manifest_sha256': sha(DEPENDENCIES),
                'environments': {}, 'readonly_mounts_enforced': False,
                'filesystem_write_protection': 'Landlock ABI>=3; task and approved output writes plus recorded POSIX IPC; inputs/databases denied',
                'container_isolation': False, 'network_or_metadata_chmod_isolation': False,
                'guard_sha256': sha(REPO / 'bin/spark_fs_guard.py'),
                'task_shell_sha256': sha(REPO / 'bin/spark-task-shell'),
                'runtime_adapter_sha256': sha(Path(__file__)),
                'posix_ipc_policy': IPC_POLICY,
                'validation': 'input_database_checksums_and_Landlock_data_write_guard'}
    for key in sorted(set(keys)):
        item, _ = environment(key)
        evidence['environments'][key] = item
    save(Path(info) / 'native_runtime_manifest.json', evidence)
    return evidence


def write_reads_software_manifest(path):
    """Give R1 preprocessing IDs the actual ARM lock, not the old image identity."""
    data = read(REPO / 'assets/reads_tools.json')
    for key, old in data.items():
        if not isinstance(old, dict) or not old.get('container'):
            continue
        item, _ = environment(key)
        old.update(legacy_container_reference=old['container'], container=None,
                   container_executed=False, runtime='native-locked-prefix',
                   platform='linux-aarch64', native_environment=key,
                   native_prefix=item['prefix'], lock_sha256=item['lock_sha256'],
                   native_manifest_sha256=sha(item['manifest']),
                   native_runtime_fingerprint=item['runtime_fingerprint'],
                   native_runtime_identity=item['runtime_identity'],
                   write_guard_sha256=sha(REPO / 'bin/spark_fs_guard.py'))
    path = Path(path)
    payload = json.dumps(data, indent=2, sort_keys=True) + '\n'
    # Preserve mtime when unchanged so resume sees the same staged path input.
    if not path.is_file() or path.read_text() != payload:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + '.tmp')
        temporary.write_text(payload); temporary.replace(path)
    return path


def reads_database_manifest(info):
    """Use server receipts only; unneeded databases may remain unconfigured."""
    root = Path(os.environ.get('CRC_PHAGE_DATABASE_ROOT', str(Path.home() / 'databases')))
    server = root / 'reads_database_manifest.spark.json'
    if server.is_file():
        return server
    path = Path(info) / 'native_databases_not_configured.json'
    payload = json.dumps({'schema_version': 1, 'runtime': 'native-locked-prefix',
                          'status': 'not_configured', 'expected_manifest': str(server)},
                         indent=2, sort_keys=True) + '\n'
    if not path.is_file() or path.read_text() != payload:
        path.write_text(payload)
    # All database-consuming branches retain db_verify(status=ready) gates.
    return path


def prepare_a0_runtime(db_root, info):
    """The same A0 database/model checks using recorded native executables."""
    import run_pipeline as a0
    allocation()
    evidence = native_record(info, ['python', 'core-genomad', 'core-checkv'])
    catalog = read(REPO / 'assets/required_assets.json')
    manifest = read(db_root / 'database_manifest.json')
    databases = a0.verify_databases(db_root, catalog, manifest)
    from spark_fs_guard import reject_ipc_storage
    reject_ipc_storage([db_root,*databases.values()])
    db = databases['checkv']
    source = db / 'genome_db/checkv_reps.faa'
    index = db / 'genome_db/checkv_reps.dmnd'
    identity = 'native-lock:core-checkv:' + environment('core-checkv')[0]['lock_sha256']
    version = a0.capture(command('core-checkv', ['diamond', 'version']))
    receipt = db_root / 'checkv_generated_index.spark.json'
    if index.exists():
        if not receipt.is_file():
            # Reused indices require an explicit cross-platform validation receipt.
            raise RuntimeError('Native CheckV index provenance missing: ' + str(receipt))
        a0.verify_index_receipt(index, receipt, version, identity, source)
        expected_binary = environment('core-checkv')[0]['native_source']['binary_sha256']
        if read(receipt).get('diamond_binary_sha256') != expected_binary:
            raise RuntimeError('Native CheckV index DIAMOND binary identity differs')
    else:
        tmp = db / 'genome_db/.checkv_reps.spark-building'
        if Path(str(tmp) + '.dmnd').exists():
            raise RuntimeError('Unfinished CheckV build exists; inspect before retrying')
        a0.stream(command('core-checkv', ['diamond', 'makedb', '--in', source,
                                        '--db', tmp, '--threads', '6'], writable=[db / 'genome_db']))
        built = Path(str(tmp) + '.dmnd')
        if not built.is_file():
            raise RuntimeError('Native DIAMOND did not produce its index')
        save(receipt, {'sha256': sha(built), 'bytes': built.stat().st_size,
                      'diamond_version': version, 'container': identity,
                      'diamond_binary_sha256': environment('core-checkv')[0]['native_source']['binary_sha256'],
                      'runtime': 'native-locked-prefix', 'source_sha256': sha(source)})
        built.replace(index)
    a0.save(info / 'checkv_generated_index.json', read(receipt))
    a0.save(info / 'database_manifest.json', manifest)
    versions = {'nextflow': a0.capture(['nextflow', '-version']),
                'java': a0.capture(['java', '-version']), 'pixi': a0.capture(['pixi', '--version']),
                'genomad': a0.capture(command('core-genomad', ['genomad', '--version'])),
                'checkv': a0.capture(command('core-checkv', ['python', '-c',
                            "import importlib.metadata;print(importlib.metadata.version('checkv'))"])),
                'runtime': 'native-locked-prefix', 'resource_budget': evidence['allocation']}
    for name in ('genomad', 'checkv'):
        if not re.search(r'(?<![0-9.])' + re.escape(catalog['tools'][name]['version']) + r'(?![0-9.])', versions[name]):
            raise RuntimeError('Native version mismatch: ' + name)
    a0.save(info / 'software_versions.json', versions)
    code = """import hashlib,json
from genomad._paths import GenomadData
p=GenomadData.data_dir
required=%r
for name in required:
    if not (p/name).is_file(): raise RuntimeError('Missing geNomad model: '+name)
files=[]
for f in sorted(p.rglob('*')):
    if f.is_file():
        h=hashlib.sha256()
        with f.open('rb') as handle:
            for block in iter(lambda:handle.read(1048576),b''):h.update(block)
        files.append(dict(path=str(f.relative_to(p)),bytes=f.stat().st_size,sha256=h.hexdigest()))
print(json.dumps(files))
""" % (a0.MODEL_FILES,)
    models = json.loads(a0.capture(command('core-genomad', ['python', '-c', code])))
    if not models:
        raise RuntimeError('Native geNomad model inventory is empty')
    a0.save(info / 'genomad_model_manifest.json', {'runtime': 'native-locked-prefix',
            'required_assets': a0.MODEL_FILES, 'files': models,
            'environment_lock_sha256': environment('core-genomad')[0]['lock_sha256']})
    return databases


def verify_viral_assets(dbroot):
    from run_pipeline import verify_databases, verify_index_receipt
    from core_functions import verify_databases as verify_functions
    paths = verify_databases(dbroot, read(REPO / 'assets/required_assets.json'),
                            read(dbroot / 'database_manifest.json'))
    from spark_fs_guard import reject_ipc_storage
    reject_ipc_storage([dbroot,*paths.values()])
    observed = subprocess.check_output(command('core-checkv', ['diamond', 'version']), text=True).strip()
    identity = 'native-lock:core-checkv:' + environment('core-checkv')[0]['lock_sha256']
    verify_index_receipt(paths['checkv'] / 'genome_db/checkv_reps.dmnd',
                         dbroot / 'checkv_generated_index.spark.json', observed, identity,
                         paths['checkv'] / 'genome_db/checkv_reps.faa')
    if read(dbroot / 'checkv_generated_index.spark.json').get('diamond_binary_sha256') != environment('core-checkv')[0]['native_source']['binary_sha256']:
        raise RuntimeError('Native CheckV index DIAMOND binary identity differs')
    functions = dbroot / 'core/functions'
    return {'genomad_db': str(paths['genomad']), 'checkv_db': str(paths['checkv']),
            'functions_db': str(functions)}, {
            'a1_manifest_sha256': sha(dbroot / 'database_manifest.json'),
            'checkv_index_manifest_sha256': sha(dbroot / 'checkv_generated_index.spark.json'),
            'functions_consumed_manifest_sha256': verify_functions(functions)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    q = sub.add_parser('shell'); q.add_argument('environment')
    q = sub.add_parser('exec'); q.add_argument('environment'); q.add_argument('command', nargs=argparse.REMAINDER)
    sub.add_parser('allocation')
    args = p.parse_args()
    allocation()
    if args.action == 'allocation':
        print(json.dumps(allocation(), indent=2))
    elif args.action == 'shell':
        _, activation = environment(args.environment)
        print('unset PYTHONPATH PYTHONHOME')
        for k, v in activation.items():
            print('export ' + k + '=' + shlex.quote(v))
    else:
        argv = args.command[1:] if args.command[:1] == ['--'] else args.command
        if not argv:
            p.error('exec requires a command')
        cmd = command(args.environment, argv)
        os.execv(cmd[0], cmd)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)
