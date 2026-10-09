#!/usr/bin/env python3
"""Auditable BLASTn HSP aggregation and deterministic representative-centred vOTUs."""
import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path

ALGORITHM = 'nonredundant_alignment_columns_reciprocal_v1'
BLAST_FIELDS = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq'


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fasta(path):
    records = {}; name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError('Duplicate FASTA ID: ' + name)
            records[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('Sequence before header')
            records[name] += line.strip().upper()
    if any(not seq for seq in records.values()):
        raise ValueError('Empty FASTA sequence')
    return records


def aggregate_hsps(hsps):
    """Greedily retain each nucleotide at most once on BOTH aligned sequences.

    Alignments are sorted by bitscore, identity, length and coordinates; full
    alignment strings resolve ties. An overlapping column is skipped if either
    non-gap nucleotide was already used. Uncovered gap columns count against
    identity. Thus repeats do not inflate ANI or AF, and input order is irrelevant.
    """
    used_q = set(); used_s = set(); identical = columns = accepted_hsps = 0
    def rank(h):
        return (-h['bitscore'], -h['pident'], -h['length'], h['qstart'], h['qend'], h['sstart'], h['send'], h['qseq'], h['sseq'])
    for h in sorted(hsps, key=rank):
        if h['evalue'] > 1e-5:
            continue
        qpos, spos = h['qstart'] - 1, h['sstart'] - 1
        qstep = 1 if h['qend'] >= h['qstart'] else -1
        sstep = 1 if h['send'] >= h['sstart'] else -1
        if len(h['qseq']) != len(h['sseq']) or len(h['qseq']) != h['length']:
            raise ValueError('Malformed BLAST alignment strings')
        nq = sum(c != '-' for c in h['qseq']); ns = sum(c != '-' for c in h['sseq'])
        if nq != abs(h['qend'] - h['qstart']) + 1 or ns != abs(h['send'] - h['sstart']) + 1:
            raise ValueError('BLAST coordinate/alignment length mismatch')
        kept = 0
        for qc, sc in zip(h['qseq'].upper(), h['sseq'].upper()):
            qi = qpos if qc != '-' else None; si = spos if sc != '-' else None
            if qi is None and si is None:
                raise ValueError('Double-gap BLAST column')
            if (qi is not None and not 0 <= qi < h['qlen']) or (si is not None and not 0 <= si < h['slen']):
                raise ValueError('BLAST coordinate out of bounds')
            if (qi is None or qi not in used_q) and (si is None or si not in used_s):
                if qi is not None:
                    used_q.add(qi)
                if si is not None:
                    used_s.add(si)
                columns += 1; kept += 1
                identical += int(qc == sc and qc in 'ACGT')
            if qi is not None:
                qpos += qstep
            if si is not None:
                spos += sstep
        accepted_hsps += int(kept > 0)
    if not hsps:
        return dict(ani_percent=0.0, af_shorter_percent=0.0, aligned_columns=0, identical_columns=0, query_covered_bp=0, subject_covered_bp=0, accepted_hsps=0)
    qlen, slen = hsps[0]['qlen'], hsps[0]['slen']
    if qlen == slen:
        af = 100 * min(len(used_q), len(used_s)) / qlen
    else:
        af = 100 * (len(used_q) / qlen if qlen < slen else len(used_s) / slen)
    return dict(ani_percent=100 * identical / columns if columns else 0.0, af_shorter_percent=af,
                aligned_columns=columns, identical_columns=identical, query_covered_bp=len(used_q),
                subject_covered_bp=len(used_s), accepted_hsps=accepted_hsps)


def parse_blast(path, sequences):
    grouped = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) != 16:
            raise ValueError('Expected BLAST outfmt 6 std qlen slen qseq sseq (16 fields)')
        h = dict(zip(BLAST_FIELDS.split(), parts))
        for key in ['length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen']:
            h[key] = int(h[key])
        for key in ['pident', 'evalue', 'bitscore']:
            h[key] = float(h[key])
        q, s = h['qseqid'], h['sseqid']
        if q not in sequences or s not in sequences or h['qlen'] != len(sequences[q]) or h['slen'] != len(sequences[s]):
            raise ValueError('BLAST ID/sequence length foreign key mismatch')
        for side, seq_id in [('q', q), ('s', s)]:
            start, end = h[side+'start'], h[side+'end']
            canonical = sequences[seq_id][min(start, end)-1:max(start, end)]
            if end < start:
                canonical = canonical.translate(str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN'))[::-1]
            if h[side+'seq'].replace('-', '').upper() != canonical:
                raise ValueError('BLAST alignment/source sequence mismatch')
        if q != s:
            grouped[q, s].append(h)
    return grouped


def pairwise(grouped, sequences, ani=95.0, af=85.0):
    rows = []; directions = []
    for a, b in itertools.combinations(sorted(sequences), 2):
        metrics = []
        for q, s in [(a, b), (b, a)]:
            m = aggregate_hsps(grouped.get((q, s), []))
            directions.append(dict(query_id=q, subject_id=s, **m))
            metrics.append(m)
        # Both directions must support the threshold. This conservative choice
        # is recorded separately from the upstream CheckV weighted-HSP route.
        ap = min(m['ani_percent'] for m in metrics)
        fp = min(m['af_shorter_percent'] for m in metrics)
        rows.append(dict(candidate_a=a, candidate_b=b, ani_percent=ap, af_shorter_percent=fp,
                         ani_threshold_percent=ani, af_threshold_percent=af,
                         passes=str(ap >= ani and fp >= af).lower(),
                         status='aligned_reciprocal' if all(m['aligned_columns'] for m in metrics) else 'no_reciprocal_alignment'))
    return rows, directions


def rank_candidate(row):
    quality = {'Complete': 0, 'High-quality': 1, 'Medium-quality': 2, 'Low-quality': 3, 'Not-determined': 4}
    boundary = {'experimentally_validated': 0, 'independently_supported': 1, 'caller_only_not_independently_validated': 2}
    return (quality.get(row.get('checkv_quality'), 5), boundary.get(row.get('boundary_confidence'), 3), -int(row['length']), row['candidate_id'])


def cluster(rows, pairs):
    adjacency = defaultdict(set)
    for pair in pairs:
        if pair['passes'] == 'true':
            a, b = pair['candidate_a'], pair['candidate_b']
            adjacency[a].add(b); adjacency[b].add(a)
    by_id = {r['candidate_id']: r for r in rows}
    bridges = set()
    for cid, row in by_id.items():
        longer = sorted(n for n in adjacency[cid] if n in by_id and int(by_id[n]['length']) > int(row['length']))
        if any(b not in adjacency[a] for a, b in itertools.combinations(longer, 2)):
            bridges.add(cid)
    # A short fragment cannot become a centroid that joins incompatible longer
    # genomes, even when its CheckV quality ranks higher. Defer such potential
    # bridges, then assess them against the independently established centroids.
    order = sorted(rows, key=lambda r: (r['candidate_id'] in bridges, -int(r['length']) if r['candidate_id'] in bridges else 0, rank_candidate(r)))
    reps = []
    for row in order:
        cid = row['candidate_id']
        if not any(rep in adjacency[cid] for rep in reps):
            reps.append(cid)
    members = []
    # Recheck against ALL final representatives. A later representative may
    # expose ambiguity that was invisible during the first greedy pass.
    for row in rows:
        cid = row['candidate_id']
        matches = [cid] if cid in reps else [r for r in reps if r in adjacency[cid]]
        if not matches:
            raise ValueError('Eligible member has no representative')
        status = 'representative' if cid in reps else 'member' if len(matches) == 1 else 'ambiguous_multiple_representatives'
        members.append(dict(candidate_id=cid, genome_id=row['genome_id'],
                            representative_id=matches[0] if len(matches) == 1 else '',
                            votu_id='votu_' + digest(matches[0].encode())[:20] if len(matches) == 1 else '',
                            membership_status=status, matching_representative_ids=';'.join(sorted(matches))))
    return reps, sorted(members, key=lambda r: r['candidate_id'])


def run(args):
    if not 50 <= args.ani_percent <= 100 or not 1 < args.af_shorter_percent <= 100:
        raise ValueError('Thresholds use percent units: expected ANI 50..100 and AF >1..100, e.g. 95 and 85; never 0.95/0.85')
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(args.master); seqs = fasta(args.fasta)
    if len({r['candidate_id'] for r in rows}) != len(rows) or {r['candidate_id'] for r in rows} != set(seqs):
        raise ValueError('Master/FASTA candidate ID foreign keys differ')
    for r in rows:
        if int(r['length']) != len(seqs[r['candidate_id']]) or digest(seqs[r['candidate_id']].encode()) != r['sequence_sha256']:
            raise ValueError('Master/FASTA sequence hash mismatch')
    qc_rows = [r for path in args.qc for r in read_tsv(path)]
    qc = {r['genome_id']: r for r in qc_rows}
    if len(qc) != len(qc_rows):
        raise ValueError('Duplicate genome QC foreign key')
    # Element classification gate (2026-10-01). Satellites, ICEs and IS clusters
    # carry real viral markers but are not prophages; leaving them in would
    # corrupt both the carriage denominator and the cargo statistics. A missing
    # table means "not assessed", never "all prophage", so nothing is silently
    # admitted. See docs/ELEMENT_CLASS_SPEC_20261001.md
    element = {}
    if getattr(args, 'element_class', None):
        ec_rows=read_tsv(args.element_class)
        for r in ec_rows:
            if r['candidate_id'] in element:raise ValueError('Duplicate element_class candidate')
            element[r['candidate_id']] = r.get('element_class', 'not_assessed')
        if set(element)-{r['candidate_id'] for r in rows}:raise ValueError('Unknown element_class candidate')
        if getattr(args,'require_element_evidence',False):
            by_id={r['candidate_id']:r for r in rows}
            versions={r.get('rule_version') for r in ec_rows}
            if ec_rows and (len(versions)!=1 or None in versions or '' in versions):raise ValueError('Missing/mixed element rule version')
            for r in ec_rows:
                if any(str(r.get(k))!=str(by_id[r['candidate_id']][k]) for k in ('genome_id','contig_id','start0','end0')):raise ValueError('Stale element candidate identity')
                if r.get('classification_status')!='assessed_provisional' and r.get('element_class')=='prophage':raise ValueError('Unassessed candidate cannot be prophage')
        missing = {r['candidate_id'] for r in rows} - set(element)
        if missing:
            raise ValueError('element_class table misses %d candidates; '
                             'a partial table cannot gate the catalog' % len(missing))
    if getattr(args,'require_element_evidence',False) and not getattr(args,'element_class',None):raise ValueError('Classification input required')
    eligible = []; filters = []; updated = []
    for row in rows:
        q = qc.get(row['genome_id'], {})
        eligibility = q.get('eligibility', 'unknown')
        keep = str(eligibility).lower() in ('true', 'yes', 'pass', 'eligible', '1')
        reason = 'host_qc_eligible' if keep else ('host_qc_not_assessed' if not q else 'host_qc_ineligible_or_unresolved')
        eclass = element.get(row['candidate_id'], 'not_assessed' if element else '')
        if keep and element and eclass != 'prophage':
            keep = False
            reason = 'element_class_' + eclass
        updated.append(dict(row, genome_qc_status=q.get('genome_qc_status', 'not_assessed'),
                            element_class=eclass,
                            main_catalog_eligible=str(keep).lower(), main_catalog_reason=reason))
        filters.append(dict(candidate_id=row['candidate_id'], genome_id=row['genome_id'],
                            element_class=eclass,
                            eligible=str(keep).lower(), reason=reason))
        if keep:
            eligible.append(row)
    grouped = parse_blast(args.blast, seqs)
    pairs, directions = pairwise(grouped, seqs, args.ani_percent, args.af_shorter_percent)
    if getattr(args, 'catalog_policy', None):
        from catalog_policy import config, bundle, build
        cfg = config(args.catalog_policy)
        if (args.ani_percent,args.af_shorter_percent)!=(cfg['ani_percent'],cfg['af_shorter_percent']):
            raise ValueError('Research policy and CLI thresholds differ')
        ev, bh = bundle(args.catalog_evidence)
        _, old_members = cluster(eligible, pairs)
        build(updated, seqs, pairs, directions, qc, cfg, ev, bh, out, old_members, master_fields=list(csv.DictReader(open(args.master), delimiter="	").fieldnames))
        return
    reps, members = cluster(eligible, pairs)
    policy = dict(algorithm=ALGORITHM, representative_rule='defer_short_bridges_then_quality_boundary_length_id_recheck_all_centroids_v1',
                  rank_evidence=[(r['candidate_id'], rank_candidate(r)) for r in sorted(rows, key=lambda r:r['candidate_id'])],
                  pairwise_results=pairs, membership_results=members, ani_percent=args.ani_percent, af_shorter_percent=args.af_shorter_percent,
                  host_qc_evidence=sorted(qc_rows, key=lambda r:r['genome_id']),
                  eligibility=sorted(filters, key=lambda r: r['candidate_id']),
                  sequences=[(r['candidate_id'], r['sequence_sha256'], r['boundary_version']) for r in sorted(rows, key=lambda r: r['candidate_id'])])
    policy = json.loads(json.dumps(policy, sort_keys=True))
    version = 'core_' + digest(json.dumps(policy, sort_keys=True, separators=(',', ':')).encode())[:20]
    by_member = {r['candidate_id']: r for r in members}
    excluded_reason = {r['candidate_id']:r['reason'] for r in filters}
    for row in rows:
        if row['candidate_id'] not in by_member:
            members.append(dict(candidate_id=row['candidate_id'], genome_id=row['genome_id'], representative_id='', votu_id='', membership_status='excluded_element_class' if excluded_reason[row['candidate_id']].startswith('element_class_') else 'excluded_host_qc', matching_representative_ids=''))
    members.sort(key=lambda r: r['candidate_id'])
    for collection in [members, updated, filters]:
        for row in collection:
            row['catalog_version'] = version
    write_tsv(out/'votu_pairwise.tsv', pairs, ['candidate_a', 'candidate_b', 'ani_percent', 'af_shorter_percent', 'ani_threshold_percent', 'af_threshold_percent', 'passes', 'status'])
    write_tsv(out/'votu_directional_ani.tsv', directions, ['query_id', 'subject_id', 'ani_percent', 'af_shorter_percent', 'aligned_columns', 'identical_columns', 'query_covered_bp', 'subject_covered_bp', 'accepted_hsps'])
    write_tsv(out/'vOTU_members.tsv', members, ['catalog_version', 'votu_id', 'candidate_id', 'genome_id', 'representative_id', 'membership_status', 'matching_representative_ids'])
    write_tsv(out/'vOTU_representatives.tsv', [r for r in members if r['membership_status'] == 'representative'], ['catalog_version', 'votu_id', 'candidate_id', 'genome_id', 'representative_id', 'membership_status', 'matching_representative_ids'])
    write_tsv(out/'catalog_filter_reasons.tsv', filters, ['catalog_version', 'candidate_id', 'genome_id', 'element_class', 'eligible', 'reason'])
    fields = list(rows[0]) if rows else list(csv.DictReader(open(args.master), delimiter='\t').fieldnames)
    write_tsv(out/'prophage_master.tsv', updated, fields + ['element_class', 'main_catalog_eligible', 'main_catalog_reason', 'catalog_version'])
    write_tsv(out/'catalog_crosswalk.tsv', [dict(old_candidate_id=r['candidate_id'], new_candidate_id=r['candidate_id'], old_boundary_version=r['boundary_version'], new_boundary_version=r['boundary_version'], relationship='unchanged_sequence_and_boundary', catalog_version=version) for r in rows], ['catalog_version', 'old_candidate_id', 'new_candidate_id', 'old_boundary_version', 'new_boundary_version', 'relationship'])
    with (out/'vOTU_representatives.fna').open('w') as handle:
        for rep in reps:
            handle.write('>' + rep + '\n' + seqs[rep] + '\n')
    manifest = dict(catalog_version=version, version_payload=policy, **policy, input_master_sha256=digest(Path(args.master).read_bytes()),
                    input_fasta_sha256=digest(Path(args.fasta).read_bytes()), blast_sha256=digest(Path(args.blast).read_bytes()),
                    helper_sha256=digest(Path(__file__).read_bytes()), candidates=len(rows), eligible_candidates=len(eligible), representatives=len(reps),
                    ambiguous=sum(r['membership_status'].startswith('ambiguous') for r in members),
                    scientific_status='provisional_engineering_policy_not_final_research_panel',
                    scope='all candidates retained; primary vOTUs restricted to measured host QC eligibility; no boundary refinement')
    (out/'catalog_manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master', required=True); parser.add_argument('--fasta', required=True)
    parser.add_argument('--blast', required=True); parser.add_argument('--qc', nargs='*', default=[])
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--ani-percent', type=float, default=95); parser.add_argument('--af-shorter-percent', type=float, default=85)
    parser.add_argument('--catalog-policy'); parser.add_argument('--catalog-evidence')
    parser.add_argument('--element-class', help='element_class.tsv from core_functions; non-prophage classes are excluded from the main catalog')
    parser.add_argument('--require-element-evidence',action='store_true')
    args = parser.parse_args(); run(args)


if __name__ == '__main__':
    main()
