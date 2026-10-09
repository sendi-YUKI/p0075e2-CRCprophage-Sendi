"""Metadata-only readiness for new opt-in core capabilities; no tool execution."""
import json
from pathlib import Path
from host_phylogeny import config
from mainline_common import sha
import spark_runtime

REPO = Path(__file__).resolve().parents[1]


def inspect_extensions(args, require_frozen=False):
    enabled_kofam = bool(getattr(args, 'enable_kofam', False))
    host_config = getattr(args, 'host_config', None)
    for name in ('host_config', 'host_mask_bundle', 'kofam_db','host_recombination_config','whole_genome_search_config'):
        value = getattr(args, name, None)
        if value and any(c in str(value) for c in ("'", '\n', '\r', '\0')):
            raise ValueError('New module paths must not contain quotes/newlines: ' + name)
    result = dict(kofam='not_enabled', host_phylogeny='not_enabled', analysis_started=False, environments={})
    from catalog_policy import inspect_policy
    result['catalog_research'] = inspect_policy(args,require_frozen)
    if require_frozen and getattr(args,'catalog_policy',None) and not spark_runtime.enabled():
        raise ValueError('Research catalog requires Spark native runtime')
    keys = []
    planned_recombination=False
    for option,envkey in [('host_recombination_config','recombination'),('whole_genome_search_config','core-blast')]:
        cp=getattr(args,option,None)
        if not cp:continue
        cp=Path(cp).resolve()
        if not cp.is_file() or cp.stat().st_size>1024*1024:raise ValueError('Missing/oversized configuration '+option)
        if option=='host_recombination_config':
            from host_recombination import config as read_config
        else:
            from whole_genome_evidence import config as read_config
            if not getattr(args,'catalog_policy',None):raise ValueError('Whole-genome search requires research catalogue policy')
        value=read_config(cp,execute=require_frozen)
        if require_frozen and value['status']!='frozen':raise ValueError('Production extension configuration must be frozen: '+option)
        if option=='host_recombination_config':planned_recombination=value['status']=='frozen' and value.get('quality_policy',{}).get('status')=='frozen' and bool(value['applicability_review'].get('evidence')) and all(value['applicability_review'].get(k) for k in ('reviewer','date','decision_ref'))
        result[option]=dict(status=value['status'],config_sha256=sha(cp),scientific_calibration=False)
        if option=='whole_genome_search_config':result[option]['callability_policy_status']=(value.get('absence_rule') or {}).get('status','not_configured')
        keys.append(envkey)
    association_config=getattr(args,'host_association_config',None)
    if association_config:
        from host_source import inspect
        for n in ('host_association_config','host_association_bundle'):
            if any(x in str(getattr(args,n,'')) for x in ("'",'\n','\r','\0')):raise ValueError('Unsafe host-source path')
        result['host_source']=inspect(association_config,getattr(args,'host_association_bundle',None),production=require_frozen,planned_recombination=planned_recombination)
        keys += ['association','python']
        if require_frozen and (not host_config or not getattr(args,'catalog_policy',None) or not spark_runtime.enabled()):
            raise ValueError('Host-source branch requires frozen host/catalog configurations and Spark native runtime')
    if enabled_kofam:
        keys += ['kofam']
        db = Path(args.kofam_db).expanduser().resolve()
        if not db.is_relative_to(Path('/srv/CRC-PHIRE/databases/kofam')):
            raise ValueError('KOfam must remain in the Spark database root')
        manifest = db / 'database_manifest.json'
        if not manifest.is_file() or manifest.stat().st_size > 16 * 1024 * 1024:
            raise ValueError('Missing/oversized KOfam manifest')
        value = json.loads(manifest.read_text())
        registry = json.loads((REPO / 'assets/spark_databases.json').read_text())['databases']['kofam']
        if value.get('status') != 'ready' or registry['status'] != 'ready' or sha(manifest) != registry['manifest_sha256'] or str(db) != registry['path']:
            raise ValueError('KOfam fixed database identity is not registered ready')
        if not (db / 'ko_list').is_file() or not (db / 'profiles').is_dir():
            raise ValueError('KOfam profiles/threshold list missing')
        result.update(kofam='metadata_ready', kofam_database=dict(path=str(db), version=value['version'], manifest_sha256=sha(manifest),
                      full_payload_verification='performed by the dedicated module once in a future authorized analysis'))
    if host_config:
        keys += ['host-core', 'host-tree']
        p = Path(host_config).expanduser().resolve()
        if not p.is_file() or p.stat().st_size > 1024 * 1024:
            raise ValueError('Missing/oversized host configuration')
        c = config(p, require_ready=require_frozen)
        if require_frozen and c['status'] != 'frozen':
            raise ValueError('Production host configuration must be explicitly frozen; synthetic is fixture-only')
        bundle = Path(args.host_mask_bundle).expanduser().resolve()
        mp = bundle / 'manifest.json'
        if not mp.is_file() or mp.stat().st_size > 1024 * 1024:
            raise ValueError('Missing/oversized independent mobile evidence manifest')
        mask = json.loads(mp.read_text())
        if mask.get('schema_version') != 1 or mask.get('purpose') != 'independent_mobile_intervals':
            raise ValueError('Unknown mobile evidence bundle')
        for row in mask['files']:
            f = (bundle / row['path']).resolve()
            if not f.is_relative_to(bundle) or not f.is_file() or f.stat().st_size > 1024 * 1024 or sha(f) != row['sha256']:
                raise ValueError('Unsafe/changed/oversized mobile coordinate metadata')
        if require_frozen and not mask['files'] and not c['mask']['use_current_a1_prophage_master']:
            raise ValueError('No independent mobile evidence source selected')
        result.update(host_phylogeny='configuration_' + c['status'], host_config=dict(path=str(p), sha256=sha(p),
            declared_genome_count=len(c['genome_ids']), mask_strategy=c['mask']['strategy'], mask_bundle_sha256=sha(mp),
            kinship='optional_host_source_consumer', scientific_choices_frozen=c['status'] == 'frozen'))
    for key in sorted(set(keys)):
        _, item, prefix, lock, _ = spark_runtime.environment_definition(key)
        result['environments'][key] = dict(prefix=str(prefix), manifest_sha256=item['manifest_sha256'], lock_sha256=sha(lock))
    if require_frozen and (host_config or enabled_kofam) and not spark_runtime.enabled():
        raise ValueError('New mainline capabilities require the Spark native launcher')
    return result
