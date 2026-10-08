#!/usr/bin/env python3
"""A3: member-level PHROGs, genome-wide DefenseFinder and conservative cargo."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict


def read_table(path, required=()):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: missing columns {required}')
        return list(reader)


def write_table(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify_databases(root):
    root = Path(root).resolve()
    manifest = json.loads((root / 'consumed_files_manifest.json').read_text())
    files = manifest['files']
    if not files or len({r['path'] for r in files}) != len(files):
        raise ValueError('Empty/duplicate database inventory')
    for row in files:
        p = (root / row['path']).resolve()
        if not p.is_relative_to(root) or not p.is_file() or p.stat().st_size != row['bytes'] or digest(p) != row['sha256']:
            raise ValueError('Consumed database asset changed: ' + row['path'])
    return digest(root / 'consumed_files_manifest.json')


def stable_id(prefix, parts):
    return prefix + hashlib.sha256(json.dumps(parts, separators=(',', ':')).encode()).hexdigest()[:32]


def fasta(path):
    seqs = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in seqs:
                raise ValueError('Duplicate protein ID: ' + name)
            seqs[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('Sequence before FASTA header')
            seqs[name] += line.strip().upper()
    return seqs


def overlap(a, b):
    if a['genome_id'] != b['genome_id'] or a['contig_id'] != b['contig_id']:
        return 0
    return max(0, min(int(a['end0']), int(b['end0'])) - max(int(a['start0']), int(b['start0'])))


def gene_relation(gene, candidate):
    n = overlap(gene, candidate)
    return ('fully_inside' if n == int(gene['end0']) - int(gene['start0']) else 'boundary_spanning') if n else 'outside'


def system_relation(required_ids, genes, candidate):
    if not required_ids or any(g not in genes for g in required_ids):
        return 'unresolved'
    relations = [gene_relation(genes[g], candidate) for g in required_ids]
    if all(x == 'fully_inside' for x in relations):
        return 'fully_inside'
    if all(x == 'outside' for x in relations):
        return 'outside'
    return 'boundary_spanning'


def membership(genes, candidates):
    return [dict(gene_id=g['gene_id'], candidate_id=c['candidate_id'], genome_id=g['genome_id'],
                 contig_id=g['contig_id'], overlap_bp=overlap(g, c), relation=gene_relation(g, c))
            for c in sorted(candidates, key=lambda x: x['candidate_id'])
            for g in sorted(genes.values(), key=lambda x: x['gene_id']) if overlap(g, c)]


HIT_FIELDS = ['hit_id', 'gene_id', 'method', 'database_version', 'phrog', 'annotation', 'category',
              'sequence_evalue', 'domain_i_evalue', 'bitscore', 'domain_bitscore', 'hmm_start1', 'hmm_end1',
              'hmm_length', 'protein_start1', 'protein_end1', 'protein_length', 'model_coverage',
              'protein_coverage', 'accepted', 'is_best', 'evidence_group']
SYSTEM_FIELDS = ['system_id', 'genome_id', 'contig_id', 'raw_system_id', 'activity', 'type', 'subtype',
                 'model_version', 'model_fqn', 'required_gene_ids', 'all_gene_ids', 'required_rule',
                 'required_resolution', 'min_mandatory_genes_required', 'min_genes_required']
SYSTEM_HIT_FIELDS = ['system_id', 'gene_id', 'gene_name', 'gene_status', 'activity', 'model_fqn',
                     'i_eval', 'score', 'profile_coverage', 'sequence_coverage', 'raw_system_id']


def phrogs_search(proteins, db, out, cpus, evalue, coverage):
    import pyhmmer
    hmm_path = db / 'phrogs_pharokka_1.8.0/pharokka_v1.8.0_databases/all_phrogs.h3m'
    ann_path = hmm_path.parent / 'phrog_annot_v4.tsv'
    annotations = {r['phrog']: r for r in read_table(ann_path, ['phrog', 'annot', 'category'])}
    if not {str(i) for i in range(1, 38881)}.issubset(annotations):
        raise ValueError('Complete PHROGs v4 annotation must cover models 1..38880')
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = pyhmmer.easel.DigitalSequenceBlock(alphabet, [pyhmmer.easel.TextSequence(name=k.encode(), sequence=v.rstrip('*')).digitize(alphabet) for k, v in sorted(proteins.items())])
    rows = []
    model_count = 0
    with pyhmmer.plan7.HMMFile(hmm_path) as hmms:
        for hits in pyhmmer.hmmsearch(hmms, sequences, cpus=cpus, E=evalue, domE=evalue):
            model_count += 1
            phrog = str(hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name).replace('phrog_', '')
            ann = annotations.get(phrog)
            if ann is None:
                raise ValueError('Unknown PHROG model ' + phrog)
            for hit in hits:
                gene = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                for dom in hit.domains:
                    ali = dom.alignment
                    mcov = (ali.hmm_to - ali.hmm_from + 1) / hits.query.M
                    pcov = (ali.target_to - ali.target_from + 1) / len(proteins[gene].rstrip('*'))
                    row = dict(hit_id=stable_id('fh_', [gene, phrog, ali.hmm_from, ali.hmm_to, ali.target_from, ali.target_to]),
                               gene_id=gene, method='PHROGs_PyHMMER', database_version='PHROGs_v4_Pharokka1.8.0',
                               phrog=phrog, annotation=ann['annot'], category=ann['category'], sequence_evalue=hit.evalue,
                               domain_i_evalue=dom.i_evalue, bitscore=hit.score, domain_bitscore=dom.score,
                               hmm_start1=ali.hmm_from, hmm_end1=ali.hmm_to, hmm_length=hits.query.M,
                               protein_start1=ali.target_from, protein_end1=ali.target_to, protein_length=len(proteins[gene].rstrip('*')),
                               model_coverage=mcov, protein_coverage=pcov,
                               accepted=str(hit.evalue <= evalue and dom.i_evalue <= evalue and mcov >= coverage).lower(),
                               is_best='false', evidence_group='PHROGs')
                    rows.append(row)
    if model_count != 38880:
        raise ValueError('Incomplete PHROGs model search: ' + str(model_count))
    best = {}
    for row in sorted(rows, key=lambda r: (r['sequence_evalue'], -r['bitscore'], r['domain_i_evalue'], r['hit_id'])):
        if row['accepted'] == 'true' and row['gene_id'] not in best:
            best[row['gene_id']] = row
            row['is_best'] = 'true'
    write_table(out / 'function_hits.tsv', sorted(rows, key=lambda r: (r['gene_id'], r['hit_id'])), HIT_FIELDS)
    return best, model_count


def model_requirements(model_fqn, model_root, gene_rows):
    parts = model_fqn.split('/')
    definitions = model_root / parts[0] / 'definitions'
    path = definitions.joinpath(*parts[1:]).with_suffix('.xml')
    if not path.exists():
        matches = list(definitions.rglob(parts[-1] + '.xml'))
        if len(matches) != 1:
            return [], 'unresolved_model', '', ''
        path = matches[0]
    root = ET.parse(path).getroot()
    roles, aliases = {}, {}
    for gene in root.findall('gene'):
        name, role = gene.attrib['name'], gene.attrib.get('presence', '')
        roles[name] = role
        aliases[name] = name
        for alternate in gene.findall('./exchangeables/gene'):
            aliases[alternate.attrib['name']] = name
    mandatory, accessory = {}, {}
    mm = int(root.get('min_mandatory_genes_required', sum(r == 'mandatory' for r in roles.values())))
    mg = int(root.get('min_genes_required', mm))
    unresolved = False
    for row in gene_rows:
        ref = row.get('hit_gene_ref') or aliases.get(row.get('gene_name', ''))
        if ref not in roles:
            ref = aliases.get(row.get('gene_name', ''))
        role = roles.get(ref)
        if role not in {'mandatory', 'accessory'} or row.get('gene_status', role).lower() != role:
            unresolved = True
            continue
        group = mandatory if role == 'mandatory' else accessory
        group.setdefault(ref, set()).add(row['hit_id'])
    ids = set().union(set(), *mandatory.values())
    if unresolved or len(mandatory) < mm:
        return sorted(ids), 'unresolved_model_role_or_mandatory_quorum', str(mm), str(mg)
    needed_accessory_roles = max(0, mg - len(mandatory))
    if needed_accessory_roles:
        if len(accessory) != needed_accessory_roles:
            # More accessory roles than needed has no unique required set. Do not
            # arbitrarily pick the ones inside the prophage to manufacture carriage.
            return sorted(ids), 'unresolved_accessory_quorum_choice', str(mm), str(mg)
        ids.update(set().union(set(), *accessory.values()))
        rule = 'model_mandatory_plus_uniquely_required_accessory_quorum'
    else:
        rule = 'model_observed_mandatory_roles_exchangeables_resolved'
    return sorted(ids), rule if ids else 'unresolved_empty_required_set', str(mm), str(mg)



def defense_search(genes, proteins, db, out, cpus):
    systems, system_hits = [], []
    by_contig = defaultdict(list)
    for g in genes.values():
        by_contig[g['contig_id']].append(g)
    from concurrent.futures import ThreadPoolExecutor
    items=sorted(by_contig.items())
    concurrency=max(1,min(int(cpus),len(items)))
    workers_per_contig=max(1,int(cpus)//concurrency)
    def execute_contig(item):
        contig,glist=item
        raw = out / 'raw/defensefinder' / contig
        raw.mkdir(parents=True, exist_ok=True)
        inp = raw / 'ordered.faa'
        inp.write_text(''.join(f">{g['gene_id']}\n{proteins[g['gene_id']]}\n" for g in sorted(glist, key=lambda x: int(x['gene_order']))))
        cmd = ['defense-finder', 'run', str(inp), '--out-dir', str(raw / 'results'), '--workers', str(workers_per_contig),
               '--db-type', 'ordered_replicon', '--preserve-raw', '--models-dir', str(db / 'defense_models_3.1.0'),
               '--antidefensefinder', '--skip-model-version-check']
        (raw / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
        with (raw / 'run.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    with ThreadPoolExecutor(max_workers=concurrency) as pool:
        list(pool.map(execute_contig,items))
    (out/'defense_execution.json').write_text(json.dumps({'contigs':len(items),'parallel_contigs':concurrency,'workers_per_contig':workers_per_contig,'cpu_budget':int(cpus),'replicon_boundaries':'preserved_separate_invocations'},indent=2)+'\n')
    for contig,glist in items:
        raw=out/'raw/defensefinder'/contig
        sf = list((raw / 'results').glob('*_defense_finder_systems.tsv'))
        gf = list((raw / 'results').glob('*_defense_finder_genes.tsv'))
        if len(sf) != 1 or len(gf) != 1:
            raise ValueError('DefenseFinder output missing/ambiguous')
        grows = read_table(gf[0], ['hit_id', 'sys_id', 'gene_name', 'hit_status', 'model_fqn'])
        for row in grows:
            row.update(gene_status=row['hit_status'], i_eval=row['hit_i_eval'], score=row['hit_score'], profile_coverage=row['hit_profile_cov'], sequence_coverage=row['hit_seq_cov'])
        for system in read_table(sf[0], ['sys_id', 'activity', 'type', 'subtype']):
            sh = [r for r in grows if r['sys_id'] == system['sys_id']]
            if not sh or any(r['hit_id'] not in genes or genes[r['hit_id']]['contig_id'] != contig for r in sh):
                raise ValueError('Unmapped system gene or contig mismatch')
            ids, rule, mm, mg = model_requirements(sh[0]['model_fqn'], db / 'defense_models_3.1.0', sh)
            sid = stable_id('sys_', [glist[0]['genome_id'], contig, system['activity'], sh[0]['model_fqn'], sorted(r['hit_id'] for r in sh), '3.1.0'])
            systems.append(dict(system_id=sid, genome_id=glist[0]['genome_id'], contig_id=contig,
                                raw_system_id=system['sys_id'], activity=system['activity'], type=system['type'], subtype=system['subtype'],
                                model_version='3.1.0', model_fqn=sh[0]['model_fqn'], required_gene_ids=','.join(ids),
                                all_gene_ids=','.join(sorted({r['hit_id'] for r in sh})), required_rule=rule,
                                required_resolution='resolved' if ids and not rule.startswith('unresolved') else 'unresolved', min_mandatory_genes_required=mm, min_genes_required=mg))
            for row in sh:
                system_hits.append(dict(row, system_id=sid, gene_id=row['hit_id'], raw_system_id=system['sys_id']))
    write_table(out / 'systems.tsv', sorted(systems, key=lambda r: r['system_id']), SYSTEM_FIELDS)
    write_table(out / 'system_hits.tsv', sorted(system_hits, key=lambda r: (r['system_id'], r['gene_id'])), SYSTEM_HIT_FIELDS)
    return systems


# --- element classification (2026-10-01) -------------------------------------
# Separates phage satellites, ICEs and IS clusters from genuine prophages before
# vOTU clustering. Without this layer both the carriage denominator and the
# defense / anti-defense cargo statistics are contaminated: satellites
# (PICI / SaPI / cfPICI / P4-like) carry real viral markers, and Bacteroides
# genomes are rich in CTnDOT-family ICEs whose load tracks host lineage.
# Spec and limitations: docs/ELEMENT_CLASS_SPEC_20261001.md
ELEMENT_FIELDS = ['candidate_id', 'genome_id', 'contig_id', 'start0', 'end0', 'length_bp',
                  'element_class', 'element_class_evidence', 'n_genes_inside', 'n_integrase',
                  'n_head_packaging', 'n_tail', 'n_connector', 'n_lysis', 'n_transposase',
                  'p4like_morphology', 'satellite_markers', 'n_satellite_markers',
                  'transposase_fraction', 'conjscan_overlap_fraction', 'conjscan_overlap_genes',
                  'conjscan_systems',
                  'preliminary_class', 'classification_status', 'phrogs_search_state', 'conjscan_search_state', 'marker_search_state', 'conjscan_method', 'rule_version']
# Project starting points, NOT published thresholds. Freeze after manual review.
SATELLITE_MAX_BP = 30000          # was 25000; see positive-control panel 2026-10-02
# P4-like morphology: integrase, no capsid, no tail, modest gene count. Shared
# with defective prophages and IMEs, so this is a review flag, never a class.
P4LIKE_MIN_GENES = 8
P4LIKE_MAX_GENES = 20
IS_TRANSPOSASE_FRACTION = 0.40
ICE_REQUIRED_OVERLAP = 0.50
ICE_AMBIGUOUS_LOWER = 0.15      # one incidental gene is not evidence of an ICE
ICE_AMBIGUOUS_MIN_GENES = 2
MIN_GENES_FOR_STRUCTURE = 10   # below this a missing module may just be a short call
ELEMENT_RULE_VERSION = 'structural_v6_states_20261002'
TRANSPOSASE_PAT = re.compile(r'transposase|\bIS[0-9]|insertion sequence|integrase/transposase', re.I)


def _conj_models(models_root):
    """Parse CONJScan XML definitions. Returns {model_fqn: spec}."""
    out = {}
    root = Path(models_root) / 'CONJScan' / 'definitions'
    if not root.exists():
        return out
    for xml in sorted(root.rglob('*.xml')):
        try:
            tree = ET.parse(xml).getroot()
        except ET.ParseError as exc:
            raise ValueError('Invalid CONJScan XML: '+str(xml)) from exc
        mandatory, accessory, neutral, forbidden, loners, alias = [], [], [], [], set(), {}
        for gene in tree.findall('gene'):
            name = gene.get('name')
            if not name:
                continue
            fam = name
            alias[name] = fam
            for ex in gene.findall('./exchangeables/gene'):
                if ex.get('name'):
                    alias[ex.get('name')] = fam
            if gene.get('loner') in ('1','true'):
                loners.add(fam)
            presence = gene.get('presence')
            if presence not in ('mandatory','accessory','neutral','forbidden'):
                raise ValueError('Unknown model component presence: '+str(presence))
            {'mandatory':mandatory,'accessory':accessory,'neutral':neutral,'forbidden':forbidden}[presence].append(fam)
        fqn = 'CONJScan/%s/%s' % (xml.parent.name, xml.stem)
        out[fqn] = dict(
            mandatory=set(mandatory), accessory=set(accessory), neutral=set(neutral), forbidden=set(forbidden), loners=loners, alias=alias,
            max_space=int(tree.get('inter_gene_max_space', 30)),
            min_mandatory=int(tree.get('min_mandatory_genes_required', len(mandatory))),
            min_genes=int(tree.get('min_genes_required', len(mandatory))),
            model_version=tree.get('vers', ''))
    return out


def conjscan_search(genes, proteins, db, out, cpus, models_dir=None):
    """Detect conjugative / mobilisable systems with hmmsearch directly.

    Returns [{contig_id, sys_id, model_fqn, gene_ids}]. On any missing
    dependency it writes state not_assessed and returns [], so an absent
    CONJScan is never read downstream as absence of ICE.
    """
    import shutil
    raw_root = out / 'raw/conjscan'
    raw_root.mkdir(parents=True, exist_ok=True)
    models_root = Path(models_dir) if models_dir else None
    hmmsearch = shutil.which('hmmsearch')
    profiles_dir = (models_root / 'CONJScan' / 'profiles') if models_root else None
    models = _conj_models(models_root) if models_root else {}
    if not hmmsearch or not profiles_dir or not profiles_dir.is_dir() or not models:
        (raw_root / 'skipped.json').write_text(json.dumps(
            {'state': 'not_assessed',
             'reason': 'hmmsearch, CONJScan profiles or model definitions unavailable',
             'hmmsearch': hmmsearch, 'models_dir': str(models_root),
             'n_models': len(models)}, indent=2) + '\n')
        return []

    combined = raw_root / 'conjscan_profiles.hmm'
    # Two things to get right here.
    # 1. hmmsearch reports the HMM's internal NAME while the XML definitions
    #    refer to the profile FILE stem (NAME "MOBV" vs gene name "T4SS_MOBV").
    #    Without this map every model silently matches nothing.
    # 2. The CONJScan set mixes HMM format versions: 124 profiles are
    #    HMMER3/b [3.0] and one is HMMER3/f [3.3.2]. Each is valid alone, but
    #    HMMER cannot read one file that switches format partway through; it
    #    aborts with "bad file format", writes a partial table and exits 1.
    #    hmmconvert normalises them so the concatenation is uniform.
    hmm_name_to_stem, fmts = {}, set()
    hmmconvert = shutil.which('hmmconvert')
    profiles = sorted(profiles_dir.glob('*.hmm'))
    with combined.open('w') as fh:
        for hmm in profiles:
            text = hmm.read_text()
            first = text.split('\n', 1)[0].strip()
            fmts.add(first.split('[')[0].strip())
            for line in text.splitlines():
                if line.startswith('NAME'):
                    hmm_name_to_stem[line.split(None, 1)[1].strip()] = hmm.stem
                    break
            if hmmconvert:
                conv = subprocess.run([hmmconvert, str(hmm)], capture_output=True, text=True)
                if conv.returncode == 0 and conv.stdout.strip():
                    text = conv.stdout
            fh.write(text if text.endswith('\n') else text + '\n')
    (raw_root / 'profile_name_map.json').write_text(
        json.dumps(hmm_name_to_stem, indent=2, sort_keys=True) + '\n')
    if not hmmconvert and len(fmts) > 1:
        (raw_root / 'skipped.json').write_text(json.dumps(
            {'state': 'not_assessed',
             'reason': 'profile set mixes HMM format versions and hmmconvert is unavailable',
             'formats': sorted(fmts)}, indent=2) + '\n')
        return []

    found = []
    by_contig = defaultdict(list)
    for g in genes.values():
        by_contig[g['contig_id']].append(g)
    for contig, glist in sorted(by_contig.items()):
        glist = sorted(glist, key=lambda x: int(x['gene_order']))
        raw = raw_root / contig
        raw.mkdir(parents=True, exist_ok=True)
        inp = raw / 'ordered.faa'
        inp.write_text(''.join(f">{g['gene_id']}\n{proteins[g['gene_id']]}\n" for g in glist))
        tbl = raw / 'hmmsearch.tbl'
        cmd = [hmmsearch, '--cut_ga', '--cpu', str(cpus), '--noali',
               '--tblout', str(tbl), str(combined), str(inp)]
        (raw / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
        with (raw / 'run.log').open('w') as log:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
        if rc != 0 or not tbl.exists():
            n_partial = 0
            if tbl.exists():
                n_partial = sum(1 for ln in tbl.read_text().splitlines()
                                if ln.strip() and not ln.startswith('#'))
            (raw / 'state.json').write_text(json.dumps(
                {'state': 'failed', 'returncode': rc,
                 'partial_hits_discarded': n_partial,
                 'note': ('hmmsearch failed; any hits it had already written are discarded '
                          'because a partial profile sweep cannot support a negative call')},
                indent=2) + '\n')
            continue
        (raw / 'state.json').write_text(json.dumps({'state':'completed','returncode':0})+'\n')
        # best profile per gene, by bitscore
        best_hit = {}
        for line in tbl.read_text().splitlines():
            if line.startswith('#') or not line.strip():
                continue
            f = line.split()
            if len(f) < 6:
                continue
            gene_id, profile, score = f[0], f[2], float(f[5])
            if gene_id not in best_hit or score > best_hit[gene_id][1]:
                best_hit[gene_id] = (profile, score)
        if not best_hit:
            continue
        order = {g['gene_id']: int(g['gene_order']) for g in glist}
        hits = sorted(((order[g], g, v[0], v[1]) for g, v in best_hit.items() if g in order))
        for fqn, spec in sorted(models.items()):
            fams = {}
            for pos, gid, prof, score in hits:
                fam = spec['alias'].get(hmm_name_to_stem.get(prof, prof))
                if fam:
                    fams[gid] = (pos, fam)
            if not fams:
                continue
            # cluster on gene order, gap <= inter_gene_max_space
            ordered = sorted(fams.items(), key=lambda kv: kv[1][0])
            clusters, cur = [], []
            for gid, (pos, fam) in ordered:
                if cur and pos - cur[-1][1][0] > spec['max_space'] + 1:
                    clusters.append(cur)
                    cur = []
                cur.append((gid, (pos, fam)))
            if cur:
                clusters.append(cur)
            loner_hits = [(gid, v) for gid, v in ordered if v[1] in spec['loners']]
            for idx, cl in enumerate(clusters):
                members = dict(cl)
                for gid, v in loner_hits:
                    members.setdefault(gid, v)
                present = {v[1] for v in members.values()}
                n_mand = len(present & spec['mandatory'])
                if (present & spec['forbidden'] or n_mand < spec['min_mandatory'] or
                        len(present & (spec['mandatory'] | spec['accessory'])) < spec['min_genes']):
                    continue
                gene_ids = sorted(g for g,v in members.items() if v[1] in (spec['mandatory'] | spec['accessory']))
                found.append(dict(contig_id=contig,
                                  sys_id='%s_%s_%d' % (contig, Path(fqn).name, idx),
                                  model_fqn=fqn, gene_ids=gene_ids,
                                  mandatory_found=sorted(present & spec['mandatory']),
                                  families_found=sorted(present)))
    all_fams = set()
    for spec in models.values():
        all_fams |= set(spec['alias'])
    unmapped = sorted(set(hmm_name_to_stem.values()) - all_fams)
    (raw_root / 'summary.json').write_text(json.dumps(
        {'state': 'partial' if any(json.loads(p.read_text()).get('state')=='failed' for p in raw_root.glob('*/state.json')) else 'completed', 'method': 'exploratory_direct_hmmsearch_cut_ga',
         'official_macsyfinder_equivalent': False,
         'model_files_sha256': {str(p.relative_to(models_root)):digest(p) for p in sorted((models_root/'CONJScan').rglob('*')) if p.is_file()},
         'n_profiles': len(hmm_name_to_stem),
         'profiles_not_referenced_by_any_model': unmapped,
         'reason_not_macsyfinder': ('MacSyFinder 2.1.4 rejects CONJScan 2.1.0 model XML '
                                    '(expects vers 2.0)'),
         'n_models': len(models), 'n_systems': len(found),
         'caveat': ('per-cluster rule check; MacSyFinder global best-solution scoring is not '
                    'reproduced, so a region may satisfy more than one model')},
        indent=2, ensure_ascii=False) + '\n')
    return found

# P4-like / capsid-less satellite markers. Accession -> short name, for the
# evidence string. See module docstring for why HTH_3 is not in this set.
SATELLITE_MARKER_PFAM = {
    'PF07455': 'Psu',
    'PF04606': 'Ogr_Delta',
    'PF03288': 'Pox_D5_alpha',
    'PF10554': 'Phage_ASH',
    'PF05930': 'Phage_AlpA',
}
SATELLITE_MARKER_MIN = 2        # one lone hit is common in ordinary prophages


def resolve_pfam_db(args):
    """Locate a Pfam-A HMM file, or return None so the search records not_assessed.

    Order: an explicit --pfam-db, then the VIBRANT-bundled Pfam-A v32 under the
    database root, which is what exists on this host. Returning None rather than
    guessing is deliberate; a wrong path would look like "no markers found".
    """
    explicit = getattr(args, 'pfam_db', None)
    if explicit and Path(explicit).exists():
        return explicit
    root = Path(getattr(args, 'database_root', '') or '')
    for pat in ('vibrant/*/databases/Pfam-A*.HMM', 'pfam/Pfam-A*.hmm',
                'pfam/Pfam-A*.HMM'):
        for cand in sorted(root.glob(pat)):
            if cand.is_file():
                return str(cand)
    return None


def satellite_marker_search(proteins, pfam_db, out, cpus):
    """Return {gene_id: sorted[marker names]} for the P4-like marker families.

    On any missing dependency this writes state not_assessed and returns {},
    so an absent Pfam is never read downstream as absence of markers.
    """
    import shutil
    raw = out / 'raw/satellite_markers'
    raw.mkdir(parents=True, exist_ok=True)
    hmmfetch = shutil.which('hmmfetch')
    hmmsearch = shutil.which('hmmsearch')
    db = Path(pfam_db) if pfam_db else None
    if not (hmmfetch and hmmsearch and db and db.exists()):
        (raw / 'skipped.json').write_text(json.dumps(
            {'state': 'not_assessed',
             'reason': 'hmmfetch, hmmsearch or the Pfam database is unavailable',
             'hmmfetch': hmmfetch, 'hmmsearch': hmmsearch,
             'pfam_db': str(db)}, indent=2) + '\n')
        return {}

    # hmmfetch needs an .ssi index; build it beside our own copy rather than
    # writing into the shared database directory.
    local_db = raw / 'pfam_subset.hmm'
    accs = sorted(SATELLITE_MARKER_PFAM)
    keyfile = raw / 'accessions.txt'
    keyfile.write_text('\n'.join(accs) + '\n')
    idx = db.with_suffix(db.suffix + '.ssi')
    fetch_src = str(db)
    if not idx.exists():
        # No index and the database is not ours to index: fall back to a linear
        # scan, which is slower but needs no write access and no index.
        text, keep, acc_of = [], False, None
        cur = []
        for line in db.read_text(errors='replace').splitlines(True):
            cur.append(line)
            if line.startswith('ACC '):
                acc_of = line.split()[1].split('.')[0]
            if line.startswith('//'):
                if acc_of in SATELLITE_MARKER_PFAM:
                    text.extend(cur)
                cur, acc_of = [], None
        local_db.write_text(''.join(text))
    else:
        r = subprocess.run([hmmfetch, '-f', fetch_src, str(keyfile)],
                           capture_output=True, text=True)
        if r.returncode != 0 or not r.stdout.strip():
            (raw / 'skipped.json').write_text(json.dumps(
                {'state': 'failed', 'reason': 'hmmfetch failed',
                 'stderr': r.stderr[-2000:]}, indent=2) + '\n')
            return {}
        local_db.write_text(r.stdout)

    found = set()
    for line in local_db.read_text(errors='replace').splitlines():
        if line.startswith('ACC '):
            found.add(line.split()[1].split('.')[0])
    missing = [a for a in accs if a not in found]
    if len(found) == 0:
        (raw / 'skipped.json').write_text(json.dumps(
            {'state': 'not_assessed',
             'reason': 'no marker profile could be extracted from the Pfam database',
             'requested': accs}, indent=2) + '\n')
        return {}

    faa = raw / 'proteins.faa'
    with faa.open('w') as fh:
        for gid, seq in sorted(proteins.items()):
            fh.write('>%s\n%s\n' % (gid, seq))
    tbl = raw / 'markers.tbl'
    cmd = [hmmsearch, '--cut_ga', '--cpu', str(max(1, int(cpus))),
           '--tblout', str(tbl), str(local_db), str(faa)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    # --cut_ga needs GA lines; Pfam has them, but fall back rather than return
    # a silent empty result if this build of the profiles lacks them.
    if r.returncode != 0 or not tbl.exists():
        cmd = [hmmsearch, '-E', '1e-5', '--cpu', str(max(1, int(cpus))),
               '--tblout', str(tbl), str(local_db), str(faa)]
        r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not tbl.exists():
        (raw / 'skipped.json').write_text(json.dumps(
            {'state': 'failed', 'reason': 'hmmsearch failed',
             'returncode': r.returncode, 'stderr': r.stderr[-2000:]}, indent=2) + '\n')
        return {}

    hits = {}
    for line in tbl.read_text(errors='replace').splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        gene, acc = f[0], f[3].split('.')[0]
        name = SATELLITE_MARKER_PFAM.get(acc)
        if name:
            hits.setdefault(gene, set()).add(name)
    out_hits = {g: sorted(v) for g, v in hits.items()}
    (raw / 'summary.json').write_text(json.dumps(
        {'state': 'assessed', 'n_profiles_used': len(found),
         'profiles_missing': missing, 'n_proteins': len(proteins),
         'n_genes_with_marker': len(out_hits),
         'markers': SATELLITE_MARKER_PFAM,
         'min_markers_for_path_b': SATELLITE_MARKER_MIN}, indent=2) + '\n')
    return out_hits


def classify_elements(candidates, members, genes, best, conj_systems, out,
                      marker_hits=None, assessment=None):
    """Assign element_class to every prophage candidate.

    Priority: ICE_candidate > satellite_candidate > IS_cluster > prophage.
    island_ambiguous is reserved for partial ICE overlap and forces manual review.
    """
    conj_by_contig = defaultdict(list)
    for sysrec in conj_systems or []:
        conj_by_contig[sysrec['contig_id']].append(sysrec)
    by_cand = defaultdict(list)
    for m in members:
        if m.get('relation') == 'fully_inside':
            by_cand[m['candidate_id']].append(m['gene_id'])
    rows = []
    for c in candidates:
        inside = by_cand.get(c['candidate_id'], [])
        cats = defaultdict(int)
        n_tp = 0
        for gid in inside:
            h = best.get(gid)
            if not h:
                continue
            cats[h.get('category', 'unknown function')] += 1
            if TRANSPOSASE_PAT.search(str(h.get('annotation', ''))):
                n_tp += 1
        n_genes = len(inside)
        length = int(c['end0']) - int(c['start0'])
        tp_frac = (n_tp / n_genes) if n_genes else 0.0
        # ICE overlap: fraction of a CONJScan system's genes falling inside this candidate
        best_ov, best_n, ov_sys = 0.0, 0, []
        for sysrec in conj_by_contig.get(c['contig_id'], []):
            gids = sysrec['gene_ids']
            if not gids:
                continue
            shared = len(set(gids) & set(inside))
            frac = shared / len(gids)
            if frac > 0:
                ov_sys.append('%s(%.2f,n=%d)' % (sysrec['model_fqn'] or sysrec['sys_id'],
                                                 frac, shared))
            if frac > best_ov:
                best_ov, best_n = frac, shared
        n_int = cats.get('integration and excision', 0)
        n_head = cats.get('head and packaging', 0)
        n_tail = cats.get('tail', 0)
        n_conn = cats.get('connector', 0)
        n_lys = cats.get('lysis', 0)
        # Recorded for every candidate so the overlap with structurally_incomplete
        # is auditable; it never changes the class on its own.
        p4like_shape = (length < SATELLITE_MAX_BP and n_int > 0
                        and n_head == 0 and n_tail == 0
                        and P4LIKE_MIN_GENES <= n_genes <= P4LIKE_MAX_GENES)
        marks = set()
        for g in inside:
            marks.update((marker_hits or {}).get(g, ()))
        marker_ok = len(marks) >= SATELLITE_MARKER_MIN
        if best_ov >= ICE_REQUIRED_OVERLAP:
            klass, why = 'ICE_candidate', 'exploratory_conjscan_quorum_gene_overlap>=%.2f' % ICE_REQUIRED_OVERLAP
        elif best_ov >= ICE_AMBIGUOUS_LOWER and best_n >= ICE_AMBIGUOUS_MIN_GENES:
            klass = 'island_ambiguous'
            why = ('partial_conjscan_overlap=%.2f_genes=%d;manual_review'
                   % (best_ov, best_n))
        elif (length < SATELLITE_MAX_BP and n_int > 0 and n_head > 0 and n_tail == 0):
            klass, why = ('satellite_candidate',
                          'pathA:short(<%dbp)+integrase+capsid+no_tail' % SATELLITE_MAX_BP)
        elif p4like_shape and marker_ok:
            # Path B. P4-like satellites encode no capsid, so path A cannot see
            # them; family markers carry the evidence instead of structure.
            klass, why = ('satellite_candidate',
                          'pathB:capsid-less+integrase+markers(%s)' % ','.join(sorted(marks)))
        elif n_genes and tp_frac >= IS_TRANSPOSASE_FRACTION and n_head == 0 and n_tail == 0:
            klass, why = 'IS_cluster', 'transposase_fraction=%.2f_no_capsid_no_tail' % tp_frac
        elif n_genes >= MIN_GENES_FOR_STRUCTURE and n_head == 0 and n_tail == 0 and n_conn == 0:
            # No structural module at all. With an integrase this is a defective
            # prophage or an integrative mobilizable element; without one it is an
            # unresolved region. Either way it cannot support a carriage call, so
            # it is singled out rather than counted as a prophage. This is the
            # operational form of the protocol's "残缺元件单列" rule (section 8.3).
            klass = 'structurally_incomplete'
            why = ('integrase_present_but_no_capsid_tail_connector'
                   if n_int else 'no_structural_module_and_no_integrase')
            if p4like_shape:
                why += (';p4like_morphology_markers=%d_below_min_%d'
                        % (len(marks), SATELLITE_MARKER_MIN))
        else:
            klass = 'prophage'
            why = ('short_with_capsid_and_tail' if length < SATELLITE_MAX_BP and n_head and n_tail
                   else 'default')
        preliminary = klass
        states = (assessment or {}).get(c['contig_id'], {})
        ph_state = states.get('phrogs', 'not_assessed')
        conj_state = states.get('conjscan', 'not_assessed')
        marker_state = states.get('markers', 'not_assessed')
        complete = all(x in ('completed_with_hits','completed_zero_hits') for x in (ph_state,conj_state,marker_state))
        if not complete or not n_genes:
            klass = 'not_assessed'
            why = 'incomplete_search_or_gene_coverage; preliminary='+preliminary
        elif klass == 'prophage' and not (n_head and (n_tail or n_conn)):
            klass, why = 'unclassified', 'insufficient_positive_structural_evidence'
        rows.append(dict(candidate_id=c['candidate_id'], genome_id=c.get('genome_id', ''),
                         contig_id=c['contig_id'], start0=c['start0'], end0=c['end0'],
                         length_bp=length, element_class=klass, element_class_evidence=why,
                         n_genes_inside=n_genes, n_integrase=n_int, n_head_packaging=n_head,
                         n_tail=n_tail, n_connector=n_conn, n_lysis=n_lys, n_transposase=n_tp,
                         p4like_morphology='yes' if p4like_shape else 'no',
                         satellite_markers=','.join(sorted(marks)),
                         n_satellite_markers=len(marks),
                         transposase_fraction='%.3f' % tp_frac,
                         conjscan_overlap_fraction='%.3f' % best_ov,
                         conjscan_overlap_genes=best_n,
                         conjscan_systems=';'.join(sorted(ov_sys)),
                         preliminary_class=preliminary, classification_status='assessed_provisional' if complete and n_genes else 'not_assessed',
                         phrogs_search_state=ph_state, conjscan_search_state=conj_state, marker_search_state=marker_state,
                         conjscan_method=states.get('conjscan_method','not_assessed'), rule_version=ELEMENT_RULE_VERSION))
    write_table(out / 'element_class.tsv',
                sorted(rows, key=lambda r: r['candidate_id']), ELEMENT_FIELDS)
    counts = defaultdict(int)
    for r in rows:
        counts[r['element_class']] += 1
    (out / 'element_class_summary.json').write_text(json.dumps(
        {'rule_version': ELEMENT_RULE_VERSION,
         'thresholds': {'satellite_max_bp': SATELLITE_MAX_BP,
                         'p4like_gene_range': [P4LIKE_MIN_GENES, P4LIKE_MAX_GENES],
                         'satellite_marker_min': SATELLITE_MARKER_MIN,
                         'satellite_markers': SATELLITE_MARKER_PFAM,
                        'is_transposase_fraction': IS_TRANSPOSASE_FRACTION,
                        'ice_required_overlap': ICE_REQUIRED_OVERLAP,
                        'ice_ambiguous_lower': ICE_AMBIGUOUS_LOWER,
                        'ice_ambiguous_min_genes': ICE_AMBIGUOUS_MIN_GENES,
                        'min_genes_for_structure': MIN_GENES_FOR_STRUCTURE},
         'conjscan_available': bool(assessment) and all(v.get('conjscan') in ('completed_with_hits','completed_zero_hits') for v in assessment.values()),
         'search_assessment': assessment or {},
         'official_macsyfinder_equivalent': False,
         'counts': dict(counts),
         'caveat': ('satellite and IS rules are structural heuristics, not model-based; '
                    'thresholds are project starting points pending manual calibration')},
        indent=2, ensure_ascii=False) + '\n')
    return rows


def annotate_members(genes, candidates, members, best, systems, out, phrogs_state='completed', defense_state='completed'):
    contexts, cargo, system_members, summary = [], [], [], []
    # Explicit narrow provisional whitelist excludes viral nucleotide replication enzymes.
    metabolic = {'photosystem ii protein d1', 'photosystem ii protein d2', 'transaldolase',
                 'transketolase', 'phosphoadenosine phosphosulfate reductase', 'cysteine synthase',
                 'phosphomannose isomerase'}
    hallmark = {'head and packaging', 'connector', 'tail', 'lysis'}
    for c in candidates:
        ms = [m for m in members if m['candidate_id'] == c['candidate_id']]
        inside = {m['gene_id'] for m in ms if m['relation'] == 'fully_inside'}
        nearby = [genes[g] for g in inside if g in best and best[g]['category'] in hallmark]
        for m in ms:
            g, h = genes[m['gene_id']], best.get(m['gene_id'])
            distance = min(int(g['start0']) - int(c['start0']), int(c['end0']) - int(g['end0']))
            neighbor_n = sum(abs(int(g['gene_order']) - int(n['gene_order'])) <= 5 and n['gene_id'] != g['gene_id'] for n in nearby)
            known = phrogs_state == 'completed'
            contexts.append(dict(candidate_id=c['candidate_id'], gene_id=g['gene_id'], relation=m['relation'],
                                 annotation=h['annotation'] if h else '', phrog=h['phrog'] if h else '',
                                 category=h['category'] if h else '', status='present' if h else ('absent_assessable' if known and m['relation'] == 'fully_inside' and g['partial'] == '00' else 'unknown'),
                                 status_scope='PHROGs_hit_under_fixed_thresholds', boundary_distance_bp=distance,
                                 nearby_viral_hallmark_genes_5=neighbor_n))
            if h:
                ismet = h['annotation'].strip().lower() in metabolic
                contam = c.get('contamination', c.get('checkv_contamination', 'NA'))
                try:
                    lowcontam = float(contam) <= 5
                except (ValueError, TypeError):
                    lowcontam = False
                supported = ismet and neighbor_n >= 2 and distance >= 1000 and lowcontam and m['relation'] == 'fully_inside' and g['partial'] == '00'
                mixed = h['category'] == 'moron, auxiliary metabolic gene and host takeover'
                cargo_class = 'viral_context_supported_candidate_AMG' if supported else ('putative_metabolic_gene' if ismet else ('other_cargo' if mixed else 'PHROG_function'))
                cargo.append(dict(candidate_id=c['candidate_id'], gene_id=g['gene_id'], cargo_class=cargo_class,
                                  function=h['annotation'], phrog=h['phrog'], category=h['category'],
                                  domain='PHROGs_profile_not_independent_Pfam', model_coverage=h['model_coverage'], protein_coverage=h['protein_coverage'],
                                  evalue=h['sequence_evalue'], relation=m['relation'], nearby_viral_hallmark_genes_5=neighbor_n,
                                  boundary_distance_bp=distance, checkv_contamination=contam,
                                  viral_identity_confidence=c.get('viral_identity_confidence', 'unknown'),
                                  boundary_confidence=c.get('boundary_confidence', 'unknown'), manual_review='pending',
                                  rule_version='cargo_provisional_v1_not_scientifically_frozen'))
        for s in systems:
            if s['genome_id'] != c['genome_id']:
                continue
            relation = system_relation(s['required_gene_ids'].split(',') if s['required_gene_ids'] else [], genes, c) if s.get('required_resolution') == 'resolved' else 'unresolved'
            all_relation = system_relation(s['all_gene_ids'].split(',') if s['all_gene_ids'] else [], genes, c)
            system_members.append(dict(system_id=s['system_id'], candidate_id=c['candidate_id'], activity=s['activity'],
                                       type=s['type'], subtype=s['subtype'], relation=relation, required_gene_relation=relation, all_detected_gene_relation=all_relation,
                                       required_gene_ids=s['required_gene_ids'], required_rule=s['required_rule']))
        summary.append(dict(candidate_id=c['candidate_id'], genome_id=c['genome_id'], n_genes_overlapping=len(ms),
                            n_genes_fully_inside=len(inside), n_phrog_best_inside=sum(g in best for g in inside),
                            phrogs_status=phrogs_state, defensefinder_status=defense_state,
                            assessability='assessable_observed_sequence' if inside and phrogs_state == defense_state == 'completed' else 'unknown',
                            biological_absence_claim='not_made'))
    write_table(out/'member_functions.tsv', contexts, ['candidate_id','gene_id','relation','annotation','phrog','category','status','status_scope','boundary_distance_bp','nearby_viral_hallmark_genes_5'])
    write_table(out/'cargo_candidates.tsv', cargo, ['candidate_id','gene_id','cargo_class','function','phrog','category','domain','model_coverage','protein_coverage','evalue','relation','nearby_viral_hallmark_genes_5','boundary_distance_bp','checkv_contamination','viral_identity_confidence','boundary_confidence','manual_review','rule_version'])
    write_table(out/'system_prophage_membership.tsv', system_members, ['system_id','candidate_id','activity','type','subtype','relation','required_gene_relation','all_detected_gene_relation','required_gene_ids','required_rule'])
    write_table(out/'candidate_function_status.tsv', summary, ['candidate_id','genome_id','n_genes_overlapping','n_genes_fully_inside','n_phrog_best_inside','phrogs_status','defensefinder_status','assessability','biological_absence_claim'])


def run(args):
    out = Path(args.outdir).resolve()
    database_root = Path(args.database_root).resolve()
    if out == database_root or out.is_relative_to(database_root) or database_root.is_relative_to(out):
        raise ValueError('Output directory overlaps database directory; refusing any writes')
    out.mkdir(parents=True, exist_ok=True)
    status = out/'function_status.json'
    if status.exists() and json.loads(status.read_text()).get('status') == 'completed':
        raise ValueError('Use a new task directory; Nextflow handles cache reuse')
    status.write_text(json.dumps(dict(status='running', stage='input_validation'))+'\n')
    try:
        database_inventory_sha256 = verify_databases(Path(args.database_root))
        gl = read_table(Path(args.genes_dir)/'gene_table.tsv', ['gene_id','genome_id','contig_id','start0','end0','gene_order','partial'])
        genes = {r['gene_id']: r for r in gl}
        if len(genes) != len(gl) or not gl or len({r['genome_id'] for r in gl}) != 1:
            raise ValueError('Exactly one nonempty genome; unique gene IDs required')
        proteins = fasta(Path(args.genes_dir)/'genes.faa')
        if set(proteins) != set(genes):
            raise ValueError('Gene/protein foreign key mismatch')
        for k, seq in proteins.items():
            if hashlib.sha256(seq.rstrip('*').encode()).hexdigest() != genes[k]['protein_sha256']:
                raise ValueError('Protein hash mismatch: ' + k)
        candidates = [c for c in read_table(args.master, ['candidate_id','genome_id','contig_id','start0','end0']) if c['genome_id'] == gl[0]['genome_id']]
        members = membership(genes, candidates)
        write_table(out/'gene_prophage_membership.tsv', members, ['gene_id','candidate_id','genome_id','contig_id','overlap_bp','relation'])
        selected = {m['gene_id']: proteins[m['gene_id']] for m in members}
        if selected:
            best, models = phrogs_search(selected, Path(args.database_root), out, args.cpus, args.phrog_evalue, args.phrog_model_coverage)
        else:
            best, models = {}, 0
            write_table(out/'function_hits.tsv', [], HIT_FIELDS)
        systems = defense_search(genes, proteins, Path(args.database_root), out, args.cpus)
        candidate_contigs={c['contig_id'] for c in candidates}
        conj_genes={k:g for k,g in genes.items() if g['contig_id'] in candidate_contigs}
        conj = conjscan_search(conj_genes, proteins, Path(args.database_root), out, args.cpus,
                               models_dir=getattr(args, 'conjscan_models_dir', None))
        marker_hits = satellite_marker_search(
            proteins, resolve_pfam_db(args), out, args.cpus)
        classify_elements(candidates, members, genes, best, conj, out,
                          marker_hits=marker_hits, assessment=search_assessment(out, candidates, genes, best, conj, marker_hits))
        annotate_members(genes, candidates, members, best, systems, out)
        manifest = dict(status='completed', genome_id=gl[0]['genome_id'], input_gene_sha256=digest(Path(args.genes_dir)/'gene_table.tsv'),
                        input_protein_sha256=digest(Path(args.genes_dir)/'genes.faa'), candidate_master_sha256=digest(args.master),
                        software={p:importlib.metadata.version(p) for p in ['mdmparis-defense-finder','macsyfinder','pyhmmer']},
                        parameters=dict(cpus=args.cpus,phrog_evalue=args.phrog_evalue,phrog_model_coverage=args.phrog_model_coverage),
                        database_manifest_sha256=digest(Path(args.database_root)/'database_manifest.json'),
                        consumed_files_manifest_sha256=database_inventory_sha256,
                        phrog_models_searched=models, n_genome_genes=len(genes), n_candidate_genes=len(selected), n_candidates=len(candidates),
                        n_phrog_best=len(best), n_defense_systems=sum(s['activity']=='Defense' for s in systems),
                        n_antidefense_systems=sum(s['activity']=='Antidefense' for s in systems),
                        outputs={p.name:digest(p) for p in out.glob('*.tsv')})
        status.write_text(json.dumps(manifest,indent=2)+'\n')
    except Exception as exc:
        status.write_text(json.dumps(dict(status='failed', error=repr(exc)),indent=2)+'\n')
        raise


def summarize(args):
    members = read_table(args.votu_members)
    observations, states, systems = [], {}, []
    for directory in args.functions_dirs:
        root = Path(directory)
        receipt = json.loads((root/'function_status.json').read_text())
        if receipt.get('status') != 'completed':
            raise ValueError('Refusing incomplete function output: '+str(root))
        for row in read_table(root/'candidate_function_status.tsv'):
            states[row['candidate_id']] = row
        observations.extend(read_table(root/'member_functions.tsv'))
        systems.extend(read_table(root/'system_prophage_membership.tsv'))
    sets = defaultdict(set)
    for row in observations:
        if row['phrog'] and row['relation'] == 'fully_inside':
            sets[row['candidate_id']].add('PHROG:'+row['phrog'])
    for row in systems:
        if row['relation'] == 'fully_inside':
            sets[row['candidate_id']].add(row['activity']+':'+row['subtype'])
    universe = sorted(set.union(set(), *sets.values()))
    long, summary, associations = [], [], []
    member_roles = {r['candidate_id']: r.get('member_role','legacy') for r in members}
    clusters = defaultdict(set)
    for row in members:
        cid = row.get('candidate_id', row.get('prophage_id'))
        vid = row.get('votu_id', row.get('vOTU_id'))
        if not cid:
            raise ValueError('vOTU membership requires candidate_id/prophage_id')
        member_status = row.get('membership_status', 'assigned')
        if not vid and member_status not in {'excluded_host_qc', 'excluded_element_class', 'ambiguous_multiple_representatives', 'unresolved_no_anchor'}:
            raise ValueError('Missing vOTU ID without explicit excluded/ambiguous membership status')
        associations.append(dict(candidate_id=cid, votu_id=vid or '', membership_status=member_status, function_assessment='completed' if cid in states else 'unknown', cluster_summary_included=str(bool(vid)).lower(), matching_representative_ids=row.get('matching_representative_ids', '')))
        if vid:
            clusters[vid].add(cid)
    for vid, cids in sorted(clusters.items()):
        for function in universe:
            counts = defaultdict(int)
            for cid in sorted(cids):
                state = 'present' if function in sets[cid] else ('absent_assessable' if states.get(cid,{}).get('assessability')=='assessable_observed_sequence' and member_roles[cid]!='fragment' else 'unknown')
                counts[state] += 1
                long.append(dict(votu_id=vid,candidate_id=cid,function_id=function,status=state,status_scope='observed_member_sequence_only'))
            denom = counts['present']+counts['absent_assessable']
            summary.append(dict(votu_id=vid,function_id=function,n_members=len(cids),n_present=counts['present'],n_absent_assessable=counts['absent_assessable'],n_unknown=counts['unknown'],
                                proportion_assessable=counts['present']/denom if denom else 'NA',proportion_all_members=counts['present']/len(cids)))
    write_table(Path(args.outdir)/'function_association_status.tsv', associations, ['candidate_id','votu_id','membership_status','function_assessment','cluster_summary_included','matching_representative_ids'])
    write_table(Path(args.outdir)/'votu_member_functions.tsv',long,['votu_id','candidate_id','function_id','status','status_scope'])
    write_table(Path(args.outdir)/'votu_function_summary.tsv',summary,['votu_id','function_id','n_members','n_present','n_absent_assessable','n_unknown','proportion_assessable','proportion_all_members'])
    Path(args.outdir,'function_summary_status.json').write_text(json.dumps(dict(status='completed',n_votus=len(clusters),n_features=len(universe),n_associations=len(associations),n_excluded_or_ambiguous=sum(not a['votu_id'] for a in associations),feature_universe='observed_unique_PHROG_and_system_models',input_votu_sha256=digest(args.votu_members)),indent=2)+'\n')


def search_assessment(out, candidates, genes, best, conj, marker_hits):
    """Per-contig execution evidence; never infer successful search from an empty list."""
    out=Path(out)
    def load(p):return json.loads(p.read_text()) if p.is_file() else {}
    root=out/'raw/conjscan'; cs=load(root/'summary.json'); skipped=load(root/'skipped.json')
    ms=load(out/'raw/satellite_markers/summary.json'); mskip=load(out/'raw/satellite_markers/skipped.json')
    result={}
    for contig in {c['contig_id'] for c in candidates}:
        ids={g for g,r in genes.items() if r['contig_id']==contig}
        state=load(root/contig/'state.json').get('state', skipped.get('state','not_assessed'))
        if state=='completed':state='completed_with_hits' if any(r['contig_id']==contig for r in conj) else 'completed_zero_hits'
        marker='completed_with_hits' if any(g in marker_hits for g in ids) else 'completed_zero_hits'
        if ms.get('state')!='assessed' or ms.get('profiles_missing'):marker=mskip.get('state','not_assessed')
        result[contig]={'phrogs':'completed_with_hits' if ids & set(best) else 'completed_zero_hits',
                        'conjscan':state,'markers':marker,'conjscan_method':cs.get('method','not_assessed')}
    return result


def collect_elements(args):
    master=read_table(args.master,['candidate_id','genome_id','contig_id','start0','end0'])
    by={r['candidate_id']:r for r in master}
    if len(by)!=len(master):raise ValueError('Duplicate master candidate')
    combined=[];seen_genomes=set()
    for d in args.functions_dirs:
        root=Path(d);receipt=json.loads((root/'function_status.json').read_text())
        if receipt.get('status')!='completed':raise ValueError('Refusing incomplete function output: '+str(root))
        if receipt.get('candidate_master_sha256')!=digest(args.master):raise ValueError('Stale candidate master receipt')
        gid=receipt['genome_id']
        if gid in seen_genomes:raise ValueError('Duplicate genome function directory')
        seen_genomes.add(gid)
        p=root/'element_class.tsv'
        if receipt.get('outputs',{}).get(p.name)!=digest(p):raise ValueError('Changed classification output')
        rs=read_table(p,ELEMENT_FIELDS)
        if {r['candidate_id'] for r in rs}!={r['candidate_id'] for r in master if r['genome_id']==gid}:raise ValueError('Classification coverage differs')
        for r in rs:
            if r['rule_version']!=ELEMENT_RULE_VERSION:raise ValueError('Stale element rule version')
            if any(str(r[k])!=str(by[r['candidate_id']][k]) for k in ('genome_id','contig_id','start0','end0')):raise ValueError('Classification candidate identity mismatch')
        combined.extend(rs)
    if len(combined)!=len(by) or {r['candidate_id'] for r in combined}!=set(by):raise ValueError('Missing or duplicated classification candidates')
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    write_table(out/'element_class.tsv',sorted(combined,key=lambda r:r['candidate_id']),ELEMENT_FIELDS)
    (out/'element_collection.json').write_text(json.dumps({'status':'completed','rule_version':ELEMENT_RULE_VERSION,
        'candidate_master_sha256':digest(args.master),'candidates':len(combined),'genomes':len(seen_genomes),
        'element_class_sha256':digest(out/'element_class.tsv')},indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    r = sub.add_parser('run')
    r.add_argument('--genes-dir',required=True)
    r.add_argument('--master',required=True)
    r.add_argument('--database-root',required=True)
    r.add_argument('--outdir',required=True)
    r.add_argument('--cpus',type=int,default=6,choices=range(1,21))
    r.add_argument('--conjscan-models-dir')
    r.add_argument('--pfam-db')
    r.add_argument('--phrog-evalue',type=float,default=1e-5)
    r.add_argument('--phrog-model-coverage',type=float,default=0.4)
    s = sub.add_parser('summarize')
    s.add_argument('--votu-members',required=True)
    s.add_argument('--functions-dirs',nargs='+',required=True)
    s.add_argument('--outdir',required=True)
    c = sub.add_parser('collect-elements')
    c.add_argument('--master',required=True)
    c.add_argument('--functions-dirs',nargs='+',required=True)
    c.add_argument('--outdir',required=True)
    args = p.parse_args()
    if args.command=='run':
        if not 0<args.phrog_evalue<=1 or not 0<args.phrog_model_coverage<=1:
            p.error('E-value/coverage outside supported range')
        run(args)
    elif args.command=='collect-elements':
        collect_elements(args)
    else:
        summarize(args)


if __name__=='__main__':
    main()
