#!/usr/bin/env python3
"""Opt-in KO homology evidence on canonical genes; no AMG/activity inference."""
import argparse
import csv
from collections import defaultdict
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys

from mainline_common import fasta, new_output, read_table, run, save, sha, unique, write_table

HIT = ['genome_id', 'gene_id', 'ko', 'threshold', 'score', 'evalue', 'score_type',
       'threshold_status', 'passes_model_threshold', 'passes_kofamscan_rule', 'definition', 'database_version']
STATUS = ['genome_id', 'gene_id', 'status', 'reported_hit_count', 'accepted_ko_count', 'database_version']
LINK = ['candidate_id', 'votu_id', 'membership_status', 'genome_id', 'gene_id', 'source_host_genome',
        'relation', 'gene_status', 'ko', 'ko_evidence_status', 'database_version']


def parse_hits(path, genes, kos, version):
    result = []
    seen = set()
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        cells = next(csv.reader([line], delimiter='\t'))
        if len(cells) != 7:
            raise ValueError('Unexpected KofamScan detail-tsv columns')
        mark, gid, ko, threshold, score, evalue, definition = cells
        mark = mark.strip()
        if mark not in ('', '*') or gid not in genes or ko not in kos:
            raise ValueError('Unknown gene/KO or malformed acceptance mark')
        if (gid, ko) in seen:
            raise ValueError('Duplicate gene/KO record')
        seen.add((gid, ko))
        s, e = float(score), float(evalue)
        if not math.isfinite(s) or not math.isfinite(e) or e < 0:
            raise ValueError('Non-finite score/E-value')
        original = kos[ko]['threshold']
        known = original != '-'
        if known and (threshold == '-' or abs(float(threshold) - float(original)) > 0.11):
            raise ValueError('Threshold does not match locked KO list')
        if not known and threshold not in ('', '-'):
            raise ValueError('Unexpected threshold for thresholdless model')
        if mark == '*' and not known:
            raise ValueError('Thresholdless KO cannot be a confident assignment')
        result.append(dict(genome_id=genes[gid]['genome_id'], gene_id=gid, ko=ko,
            threshold=original, score=score, evalue=evalue, score_type=kos[ko]['score_type'],
            threshold_status='available' if known else 'no_threshold',
            passes_model_threshold=str(mark == '*').lower() if known else 'not_assessed', passes_kofamscan_rule=str(mark == '*').lower(),
            definition=definition, database_version=version))
    return result


def validate_database(root):
    root = Path(root).resolve()
    manifest = json.loads((root / 'database_manifest.json').read_text())
    if manifest.get('status') != 'ready' or manifest.get('database') != 'KOfam':
        raise ValueError('KOfam database is not ready')
    records = manifest['files']
    if len({r['path'] for r in records}) != len(records) or not records:
        raise ValueError('Empty/duplicate database inventory')
    for row in records:
        p = (root / row['path']).resolve()
        if not p.is_relative_to(root) or not p.is_file() or p.stat().st_size != row['bytes'] or sha(p) != row['sha256']:
            raise ValueError('KOfam database asset changed: ' + row['path'])
    kos = unique(read_table(root / 'ko_list', ['knum', 'threshold', 'score_type']), 'knum')
    return manifest, kos


def annotate(a):
    out = new_output(a.outdir)
    genes_dir = Path(a.genes_dir)
    genes = unique(read_table(genes_dir / 'gene_table.tsv', ['gene_id', 'genome_id', 'protein_sha256']), 'gene_id')
    proteins = fasta(genes_dir / 'genes.faa')
    if set(genes) != set(proteins):
        raise ValueError('Protein/gene catalog mismatch')
    import hashlib
    for gid, seq in proteins.items():
        if hashlib.sha256(seq.encode()).hexdigest() != genes[gid]['protein_sha256']:
            raise ValueError('Canonical protein checksum mismatch: ' + gid)
    root = Path(a.database)
    receipt = json.loads(Path(a.database_receipt).read_text())
    if receipt.get('status') != 'verified' or receipt['manifest_sha256'] != sha(root / 'database_manifest.json'):
        raise ValueError('Missing/current KOfam database verification receipt')
    db = json.loads((root / 'database_manifest.json').read_text())
    kos = unique(read_table(root / 'ko_list', ['knum', 'threshold', 'score_type']), 'knum')
    status = dict(status='running', database_version=db['version'], database_manifest_sha256=sha(Path(a.database) / 'database_manifest.json'),
                  input_faa_sha256=sha(genes_dir / 'genes.faa'), scientific_calibration=False)
    save(out / 'kofam_status.json', status)
    try:
        # Three serial searches at most in the 6-CPU profile; coordinators do not multiply HMMER worker pools.
        cpus = int(a.cpus)
        if not 1 <= cpus <= 6:
            raise ValueError('KOfam CPU budget must be 1..6')
        workers = max(1, cpus // 2)
        shim = Path(a.shim).resolve()
        hmmer = shutil.which('hmmsearch')
        if not hmmer or not shim.is_file():
            raise ValueError('Missing reviewed HMMER/shim')
        os.environ['CRC_KOFAM_HMMSEARCH'] = hmmer
        config = out / 'kofam_config.yml'
        config.write_text('hmmsearch: ' + json.dumps(str(shim)) + '\n')
        status.update(parallel_searches=workers, hmmer_cpu_argument=0, hmmer_path=hmmer, shim_sha256=sha(shim))
        if proteins:
            run(['exec_annotation', '-c', config, '-p', Path(a.database) / 'profiles', '-k', Path(a.database) / 'ko_list',
                 '--cpu', str(workers), '--tmp-dir', out / 'raw', '--keep-tabular', '-f', 'detail-tsv',
                 '-o', out / 'raw_hits.tsv', genes_dir / 'genes.faa'], out, 'kofamscan', timeout=a.timeout)
        else:
            (out / 'raw_hits.tsv').write_text('')
        hits = parse_hits(out / 'raw_hits.tsv', genes, kos, db['version'])
        write_table(out / 'gene_ko_hits.tsv', hits, HIT)
        by_gene = defaultdict(list)
        for h in hits:
            by_gene[h['gene_id']].append(h)
        states = []
        for gid, gene in sorted(genes.items()):
            hh = by_gene[gid]
            accepted = {h['ko'] for h in hh if h['passes_kofamscan_rule'] == 'true'}
            state = 'assigned' if accepted else 'no_threshold' if any(h['threshold_status'] == 'no_threshold' for h in hh) else 'below_threshold' if hh else 'no_reported_hit'
            states.append(dict(genome_id=gene['genome_id'], gene_id=gid, status=state, reported_hit_count=len(hh),
                               accepted_ko_count=len(accepted), database_version=db['version']))
        write_table(out / 'gene_ko_status.tsv', states, STATUS)
        status.update(status='completed', genes=len(genes), reported_hits=len(hits), semantics='KO homology; not AMG, pathway completeness or activity')
        save(out / 'kofam_status.json', status)
    except Exception as e:
        status.update(status='failed', error=str(e))
        write_table(out / 'gene_ko_status.tsv', [dict(genome_id=g['genome_id'], gene_id=g['gene_id'], status='failed',
                    reported_hit_count='', accepted_ko_count='', database_version=db['version']) for g in genes.values()], STATUS)
        save(out / 'kofam_status.json', status)
        raise


def summarize(a):
    out = new_output(a.outdir)
    states, hits, relations, candidates = [], [], [], {}
    for directory in a.kofam_dirs:
        p = Path(directory)
        if json.loads((p / 'kofam_status.json').read_text())['status'] != 'completed':
            raise ValueError('KOfam upstream did not complete')
        states.extend(read_table(p / 'gene_ko_status.tsv', STATUS))
        hits.extend(read_table(p / 'gene_ko_hits.tsv', HIT))
    genemap = unique(states, 'gene_id')
    by_gene = defaultdict(list)
    for h in hits:
        if h['gene_id'] not in genemap:
            raise ValueError('Unknown KO gene')
        by_gene[h['gene_id']].append(h)
    for d in a.function_dirs:
        relations.extend(read_table(Path(d) / 'gene_prophage_membership.tsv', ['candidate_id', 'gene_id', 'relation']))
        for r in read_table(Path(d) / 'candidate_function_status.tsv', ['candidate_id']):
            if r['candidate_id'] in candidates:
                raise ValueError('Duplicate candidate')
            candidates[r['candidate_id']] = r
    members = unique(read_table(a.votu_members, ['candidate_id', 'votu_id']), 'candidate_id')
    if not set(members).issubset(candidates):
        raise ValueError('Unknown vOTU candidate')
    linked = []
    accepted_sets = defaultdict(set)
    relation_seen = set()
    for r in relations:
        gid, cid = r['gene_id'], r['candidate_id']
        if gid not in genemap or cid not in candidates or (cid, gid) in relation_seen:
            raise ValueError('Invalid/duplicate gene-candidate relation')
        relation_seen.add((cid, gid))
        s = genemap[gid]
        member = members.get(cid)
        if s['genome_id'] != candidates[cid]['genome_id'] or ('genome_id' in r and r['genome_id'] != s['genome_id']):
            raise ValueError('Gene/candidate source-host identity mismatch')
        if member and 'genome_id' in member and member['genome_id'] != s['genome_id']:
            raise ValueError('vOTU member source-host identity mismatch')
        if r['relation'] not in ('fully_inside', 'boundary_spanning'):
            raise ValueError('Unknown gene/candidate overlap relation')
        for h in by_gene[gid] or [None]:
            accepted = h is not None and h['passes_kofamscan_rule'] == 'true'
            linked.append(dict(candidate_id=cid, votu_id=member['votu_id'] if member else '',
                membership_status=member.get('membership_status', 'recorded') if member else 'not_assessed',
                genome_id=s['genome_id'], source_host_genome=s['genome_id'], gene_id=gid, relation=r['relation'],
                gene_status=s['status'], ko=h['ko'] if h else '',
                ko_evidence_status='accepted' if accepted else ('no_threshold' if h['threshold_status'] == 'no_threshold' else 'below_threshold') if h else 'no_reported_hit', database_version=s['database_version']))
            if accepted:
                accepted_sets[(cid, h['ko'], r['relation'])].add(gid)
    write_table(out / 'gene_ko_hits.tsv', hits, HIT)
    write_table(out / 'gene_ko_status.tsv', states, STATUS)
    write_table(out / 'member_ko_evidence.tsv', linked, LINK)
    counts = [dict(candidate_id=cid, votu_id=members[cid]['votu_id'] if cid in members else '', ko=ko, relation=rel,
                   distinct_gene_count=len(ids), gene_ids=';'.join(sorted(ids))) for (cid, ko, rel), ids in sorted(accepted_sets.items())]
    write_table(out / 'member_ko_counts.tsv', counts, ['candidate_id', 'votu_id', 'ko', 'relation', 'distinct_gene_count', 'gene_ids'])
    candidate_genes = defaultdict(set)
    for cid, gid in relation_seen:
        candidate_genes[cid].add(gid)
    write_table(out / 'candidate_ko_status.tsv', [dict(candidate_id=cid,
        genome_id=r['genome_id'], assessed_gene_count=len(candidate_genes[cid]),
        status='assessed' if candidate_genes[cid] else 'no_overlapping_canonical_gene',
        membership_status='recorded' if cid in members else 'not_assessed', biological_absence_claim='not_made')
        for cid, r in sorted(candidates.items())], ['candidate_id', 'genome_id', 'assessed_gene_count', 'status', 'membership_status', 'biological_absence_claim'])
    observed = defaultdict(set)
    for r in counts:
        if r['votu_id']:
            observed[(r['votu_id'], r['ko'], r['relation'])].add(r['candidate_id'])
    write_table(out / 'votu_ko_observations.tsv', [dict(votu_id=v, ko=k, relation=rel, observed_member_count=len(ids), member_ids=';'.join(sorted(ids)))
                for (v, k, rel), ids in sorted(observed.items())], ['votu_id', 'ko', 'relation', 'observed_member_count', 'member_ids'])
    save(out / 'kofam_summary.json', dict(status='completed', genes=len(states), candidates=len(candidates),
         semantics='Member-specific observed KO evidence; no representative propagation; no absence/AMG claim',
         membership_source_sha256=sha(a.votu_members), outputs={p.name:sha(p) for p in out.glob('*.tsv')}, real_validation='not_performed'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    q = sub.add_parser('verify-db')
    q.add_argument('--database', required=True)
    q.add_argument('--expected-manifest-sha256', required=True)
    q.add_argument('--outdir', required=True)
    q = sub.add_parser('run')
    for n in ['genes-dir', 'database', 'database-receipt', 'shim', 'outdir']:
        q.add_argument('--' + n, required=True)
    q.add_argument('--cpus', type=int, default=1)
    q.add_argument('--timeout', type=int)
    q = sub.add_parser('summarize')
    for n in ['kofam-dirs', 'function-dirs']:
        q.add_argument('--' + n, nargs='+', required=True)
    q.add_argument('--votu-members', required=True)
    q.add_argument('--outdir', required=True)
    a = p.parse_args()
    if a.action == 'verify-db':
        out = new_output(a.outdir)
        if sha(Path(a.database) / 'database_manifest.json') != a.expected_manifest_sha256:
            raise ValueError('KOfam manifest differs from registered database identity')
        db, kos = validate_database(a.database)
        save(out / 'database_verified.json', dict(status='verified', version=db['version'],
             manifest_sha256=sha(Path(a.database) / 'database_manifest.json'), ko_count=len(kos)))
    else:
        (annotate if a.action == 'run' else summarize)(a)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print('ERROR: ' + str(e), file=sys.stderr)
        raise SystemExit(1)
