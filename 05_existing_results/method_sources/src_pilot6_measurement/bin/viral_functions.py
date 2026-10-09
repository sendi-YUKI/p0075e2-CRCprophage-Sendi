#!/usr/bin/env python3
"""Run member-sequence PHROGs and ordered viral-contig DF/ADF evidence."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import core_functions as cf
import phageflow as pf
from viral_catalog import applicable, guard, table, write


def run(args):
    out = guard(args.outdir, [args.catalog, args.genes_dir, args.database_root])
    status_path = out / 'viral_function_status.json'
    pf.write_json(status_path, dict(status='running'))
    try:
        inventory = cf.verify_databases(args.database_root)
        lock = json.loads(Path(args.tools_manifest).read_text())
        if inventory != lock['databases']['consumed_files_manifest_sha256']:
            raise ValueError('Function model inventory differs from locked tool manifest')
        gs = Path(args.genes_dir)
        receipt = json.loads((gs / 'gene_calling_status.json').read_text())
        if receipt.get('status') != 'completed':
            raise ValueError('Unsuccessful gene caller cannot become a negative annotation')
        all_genes = table(gs / 'gene_table.tsv')
        all_proteins = cf.fasta(gs / 'genes.faa')
        if {g['gene_id'] for g in all_genes} != set(all_proteins) or len({g['gene_id'] for g in all_genes}) != len(all_genes):
            raise ValueError('Gene/protein foreign-key mismatch')
        for g in all_genes:
            if pf.sequence_hash(all_proteins[g['gene_id']].rstrip('*')) != g['protein_sha256']:
                raise ValueError('Changed viral protein')
        masters = table(Path(args.catalog) / 'viral_sequence_master.tsv')
        included = {r['viral_sequence_id'] for r in masters if applicable(r['virus_scope'], 'PHROGs') in {'applicable', 'provisional'}}
        genes = {g['gene_id']: g for g in all_genes if g['viral_sequence_id'] in included}
        proteins = {g: all_proteins[g] for g in genes}
        candidates = [dict(r, candidate_id=r['viral_sequence_id'], genome_id=r['viral_sequence_id'],
                           contig_id=r['viral_sequence_id'], start0=0, end0=int(r['length']),
                           viral_identity_confidence=r.get('viral_identity_confidence') or 'unknown') for r in masters if r['viral_sequence_id'] in included]
        if proteins:
            best, models = cf.phrogs_search(proteins, Path(args.database_root), out, args.cpus, 1e-5, 0.4)
            systems = cf.defense_search(genes, proteins, Path(args.database_root), out, args.cpus)
        else:
            best, models, systems = {}, 0, []
            cf.write_table(out / 'function_hits.tsv', [], cf.HIT_FIELDS)
            cf.write_table(out / 'systems.tsv', [], cf.SYSTEM_FIELDS)
            cf.write_table(out / 'system_hits.tsv', [], cf.SYSTEM_HIT_FIELDS)
        membership = cf.membership(genes, candidates)
        cf.annotate_members(genes, candidates, membership, best, systems, out)
        write(out / 'gene_sequence_membership.tsv', membership, ['gene_id', 'candidate_id', 'contig_id', 'overlap_bp', 'relation'])
        # Preserve raw core evidence, expose viral semantics and explicit context.
        by_id = {r['viral_sequence_id']: r for r in masters}
        for name in ['member_functions.tsv', 'cargo_candidates.tsv', 'system_prophage_membership.tsv', 'candidate_function_status.tsv', 'systems.tsv']:
            fields, records = pf.read_tsv(out / name)
            for r in records:
                sid = r.get('candidate_id', r.get('contig_id', ''))
                r.update(viral_sequence_id=sid, evidence_context='viral_contig_only_no_host_flanks',
                         tool_applicability=applicable(by_id[sid]['virus_scope'], 'PHROGs'))
            write(out / name, records, fields + ['viral_sequence_id', 'evidence_context', 'tool_applicability'])
        gcounts = defaultdict(int)
        for g in genes.values():
            gcounts[g['viral_sequence_id']] += 1
        states = [dict(viral_sequence_id=r['viral_sequence_id'], status='unsupported' if r['viral_sequence_id'] not in included
                       else 'completed' if gcounts[r['viral_sequence_id']] else 'completed_zero_genes',
                       applicability=applicable(r['virus_scope'], 'PHROGs'), gene_count=gcounts[r['viral_sequence_id']],
                       assessability='assessable_observed_sequence' if r['viral_sequence_id'] in included and gcounts[r['viral_sequence_id']] else 'unknown',
                       biological_absence_claim='not_made') for r in masters]
        write(out / 'viral_sequence_function_status.tsv', states, ['viral_sequence_id', 'status', 'applicability', 'gene_count', 'assessability', 'biological_absence_claim'])
        pf.write_json(status_path, dict(status='completed', database_inventory_sha256=inventory,
                      models_searched=models, n_genes=len(genes), n_systems=len(systems),
                      identity_confidence_policy='inherit_master; missing=unknown; never infer caller confidence',
                      system_context='one viral sequence per ordered_replicon; linear; no host/circular bridge',
                      applicability_policy='viral_tool_applicability_v1', real_data_validated=False,
                      scientific_calibration=False, output_files={p.name: pf.file_hash(p) for p in out.glob('*.tsv')}, inputs={p: pf.file_hash(Path(p)) for p in
                      [str(gs / 'gene_table.tsv'), str(gs / 'genes.faa'), str(Path(args.catalog) / 'viral_sequence_master.tsv')]}))
    except Exception as exc:
        pf.write_json(status_path, dict(status='failed', error=repr(exc)))
        raise


def summarize(args):
    source, catalog = Path(args.functions), Path(args.catalog)
    out = guard(args.outdir, [source, catalog])
    receipt = json.loads((source / 'viral_function_status.json').read_text())
    if receipt.get('status') != 'completed':
        raise ValueError('Failed functions do not become absent features')
    states = {r['viral_sequence_id']: r for r in table(source / 'viral_sequence_function_status.tsv')}
    features = defaultdict(set)
    for r in table(source / 'member_functions.tsv'):
        if r['phrog'] and r['relation'] == 'fully_inside':
            features[r['candidate_id']].add('PHROG:' + r['phrog'])
    for r in table(source / 'system_prophage_membership.tsv'):
        if r['required_gene_relation'] == 'fully_inside':
            features[r['candidate_id']].add(r['activity'] + ':' + r['subtype'])
    universe = sorted(set.union(set(), *features.values()))
    groups = defaultdict(list)
    associations = table(catalog / 'vOTU_members.tsv')
    for r in associations:
        if r['votu_id']:
            groups[r['votu_id']].append(r['viral_sequence_id'])
    long, summary = [], []
    for vid, sids in sorted(groups.items()):
        for feature in universe:
            counts = defaultdict(int)
            for sid in sids:
                state = 'present' if feature in features[sid] else 'absent_assessable' if states.get(sid, {}).get('assessability') == 'assessable_observed_sequence' else 'unknown'
                counts[state] += 1
                long.append(dict(votu_id=vid, viral_sequence_id=sid, function_id=feature, status=state))
            assessed = counts['present'] + counts['absent_assessable']
            summary.append(dict(votu_id=vid, function_id=feature, n_distinct_sequences=len(sids),
                                n_present=counts['present'], n_absent_assessable=counts['absent_assessable'], n_unknown=counts['unknown'],
                                proportion_assessable=counts['present'] / assessed if assessed else 'NA',
                                proportion_all_sequences=counts['present'] / len(sids)))
    write(out / 'votu_sequence_functions.tsv', long, ['votu_id', 'viral_sequence_id', 'function_id', 'status'])
    write(out / 'votu_function_summary.tsv', summary, ['votu_id', 'function_id', 'n_distinct_sequences', 'n_present', 'n_absent_assessable', 'n_unknown', 'proportion_assessable', 'proportion_all_sequences'])
    pf.write_json(out / 'summary_manifest.json', dict(status='completed', output_files={p.name: pf.file_hash(p) for p in out.glob('*.tsv')}, feature_universe='observed PHROGs/system models',
                  aggregation_unit='distinct source-oriented viral sequence; not biological samples or independent subjects',
                  representative_broadcast=False, source_occurrences_kept_in_catalog=True,
                  excluded_or_ambiguous=sum(not r['votu_id'] for r in associations)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    r = sub.add_parser('run'); r.add_argument('--catalog', required=True); r.add_argument('--genes-dir', required=True)
    r.add_argument('--database-root', required=True); r.add_argument('--tools-manifest', required=True); r.add_argument('--outdir', required=True); r.add_argument('--cpus', type=int, choices=range(1, 7), default=6)
    r = sub.add_parser('summarize'); r.add_argument('--catalog', required=True); r.add_argument('--functions', required=True); r.add_argument('--outdir', required=True)
    args = p.parse_args(); {'run': run, 'summarize': summarize}[args.command](args)


if __name__ == '__main__':
    main()
