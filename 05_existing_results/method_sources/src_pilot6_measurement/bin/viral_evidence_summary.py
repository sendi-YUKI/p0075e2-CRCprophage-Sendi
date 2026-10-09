#!/usr/bin/env python3
"""Strict evidence aggregation; never changes primary taxonomy, vOTUs or abundance."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

from viral_evidence import load_catalog, table, write, sha
from evidence_crosswalk import OCC_REQUIRED, validate_occurrences

TOOLS = ('bacphlip', 'vibrant', 'vcontact3')
SUMMARY_FIELDS = ['viral_sequence_id', 'tool', 'status', 'votu_id', 'membership_status',
                  'amg_entry_count', 'amg_unique_tool_genes', 'amg_exact_unique_viral_genes',
                  'amg_exact_unique_host_canonical_genes', 'tool_output_directory']


def verify_receipt(root, expected_tool=None, catalog_hash=None):
    root = Path(root).resolve()
    receipt = json.loads((root / 'status.json').read_text())
    if expected_tool and receipt.get('tool') != expected_tool:
        raise ValueError('Evidence tool identity mismatch')
    if receipt.get('status') not in ('completed', 'empty_input', 'not_assessed'):
        raise ValueError('Evidence did not complete successfully')
    if catalog_hash and receipt.get('catalog_manifest_sha256') != catalog_hash:
        raise ValueError('Evidence uses another catalog')
    if not receipt.get('output_files'):
        raise ValueError('Evidence receipt lacks output hashes')
    for relative, checksum in receipt['output_files'].items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file() or sha(path) != checksum:
            raise ValueError('Evidence output missing/changed: ' + relative)
    return receipt


def summarize(catalog, evidence_dirs, enabled, outdir):
    catalog, out = Path(catalog).resolve(), Path(outdir).resolve()
    if len(enabled) != len(set(enabled)) or set(enabled) - set(TOOLS):
        raise ValueError('Unknown/duplicate requested evidence tool')
    for item in [catalog, *map(Path, evidence_dirs)]:
        item = item.resolve()
        if out == item or out.is_relative_to(item) or item.is_relative_to(out):
            raise ValueError('Summary output overlaps input')
    if out.exists():
        raise ValueError('Summary requires new output directory')
    masters, seqs = load_catalog(catalog)
    occurrences = table(catalog / 'source_occurrences.tsv', OCC_REQUIRED)
    validate_occurrences(occurrences, seqs)
    members = table(catalog / 'vOTU_members.tsv', ('viral_sequence_id', 'votu_id', 'membership_status', 'catalog_version'))
    by_member = {r['viral_sequence_id']: r for r in members}
    if len(members) != len(by_member) or set(by_member) != set(seqs):
        raise ValueError('Catalog member foreign-key mismatch')
    manifest = json.loads((catalog / 'catalog_manifest.json').read_text())
    if len({r['catalog_version'] for r in members}) > 1 or any(r['catalog_version'] != manifest['catalog_version'] for r in members):
        raise ValueError('Mixed/stale catalog version')
    catalog_hash = sha(catalog / 'catalog_manifest.json')
    tool_inputs, receipts = {}, {}
    for directory in map(Path, evidence_dirs):
        receipt = verify_receipt(directory, catalog_hash=catalog_hash)
        tool = receipt['tool']
        if tool not in enabled or tool in tool_inputs:
            raise ValueError('Unexpected/duplicate evidence directory')
        tool_inputs[tool], receipts[tool] = directory, receipt
    if set(tool_inputs) != set(enabled):
        raise ValueError('Requested evidence branch missing')
    output = []
    for tool in TOOLS:
        if tool in enabled:
            states = table(tool_inputs[tool] / 'sequence_status.tsv', ('viral_sequence_id', 'status'))
            by_state = {r['viral_sequence_id']: r['status'] for r in states}
            if len(states) != len(by_state) or set(by_state) != set(seqs):
                raise ValueError('Evidence sequence status foreign-key mismatch')
            allowed = {'completed', 'completed_no_hit', 'not_assessed_scope', 'not_assessed_incomplete'}
            if set(by_state.values()) - allowed:
                raise ValueError('Invalid terminal per-sequence state')
        else:
            by_state = {sid: 'not_enabled' for sid in seqs}
        entries, cross = [], []
        if tool == 'vibrant' and tool in enabled:
            entries = table(tool_inputs[tool] / 'amg_entries.tsv', ('viral_sequence_id', 'tool_gene_id'))
            cross = table(tool_inputs[tool] / 'gene_crosswalk.tsv', ('viral_sequence_id', 'tool_gene_id', 'target_gene_id', 'target_gene_set_id', 'relationship', 'target_namespace'))
            if any(r['viral_sequence_id'] not in seqs for r in entries + cross):
                raise ValueError('AMG/crosswalk orphan sequence')
        entries_by_sequence, exact_by_sequence = defaultdict(list), defaultdict(list)
        for entry in entries:
            entries_by_sequence[entry['viral_sequence_id']].append(entry)
        for row in cross:
            if row['relationship'] == 'exact':
                exact_by_sequence[row['viral_sequence_id']].append(row)
        for master in masters:
            sid = master['viral_sequence_id']
            row = dict(viral_sequence_id=sid, tool=tool, status=by_state[sid],
                       votu_id=by_member[sid]['votu_id'], membership_status=by_member[sid]['membership_status'],
                       amg_entry_count='', amg_unique_tool_genes='', amg_exact_unique_viral_genes='',
                       amg_exact_unique_host_canonical_genes='',
                       tool_output_directory=f'evidence/{tool}/evidence' if tool in enabled else '')
            if tool == 'vibrant' and by_state[sid] in ('completed', 'completed_no_hit'):
                local_entries = entries_by_sequence[sid]
                genes = {r['tool_gene_id'] for r in local_entries}
                exact = [r for r in exact_by_sequence[sid] if r['tool_gene_id'] in genes]
                row.update(amg_entry_count=len(local_entries), amg_unique_tool_genes=len(genes),
                           amg_exact_unique_viral_genes=len({(r['target_gene_set_id'], r['target_gene_id']) for r in exact if r['target_namespace'] == 'viral'}),
                           amg_exact_unique_host_canonical_genes=len({(r['target_gene_set_id'], r['target_gene_id']) for r in exact if r['target_namespace'] == 'host_canonical'}))
            output.append(row)
    out.mkdir(parents=True)
    write(out / 'sequence_evidence.tsv', output, SUMMARY_FIELDS)
    # Preserve source-host, predicted-host and sample evidence as separate columns.
    host_fields = OCC_REQUIRED + ['candidate_type', 'integration_evidence', 'host_resolved_prophage', 'source_genome_eligibility', 'sample_host_evidence']
    host_rows = [{key: row.get(key, 'unknown') for key in host_fields} for row in occurrences]
    write(out / 'source_host_connections.tsv', host_rows, host_fields)
    (out / 'OUTPUTS.md').write_text(
        '# Optional viral evidence\n\n'
        'sequence_evidence.tsv: one catalog sequence x tool; join catalog/viral_sequence_master.tsv by viral_sequence_id. '
        'votu_id is inherited membership, never a vConTACT network cluster.\n\n'
        'Counts distinguish AMG annotation entries, unique tool ORFs, exact unique viral genes and exact unique source canonical genes. '
        'The two target namespaces must not be added together. Zero exact matches does not imply absence of genes. '
        'Blank counts mean not applicable/not assessed; inspect status.\n\n'
        'source_host_connections.tsv: one source occurrence, including original sample/library and source assembly/contig coordinates. '
        'Source host, predicted host and sample evidence remain separate. Do not sum abundance after a one-to-many source-host join. '
        'This table does not supply or infer clinical subject/cohort metadata.\n\n'
        'Tool details: ../evidence/TOOL/evidence/ (only when enabled). VIBRANT gene_crosswalk.tsv preserves exact/partial/ambiguous/unmapped; '
        'tool_genes.tsv coordinates refer to the catalog sequence; mapped coordinates refer to the selected target namespace. '
        'All coordinates are 0-based half-open. No annotation transfer or scientific calibration is claimed.\n\n'
        'not_enabled, not_assessed_scope, not_assessed_incomplete, completed and completed_no_hit are distinct. '
        'A failed or missing requested branch prevents a successful summary.\n')
    doc = dict(schema_version=1, status='completed', enabled_tools=enabled,
               catalog_manifest_sha256=catalog_hash, catalog_version=manifest['catalog_version'],
               sequence_count=len(seqs), source_occurrence_count=len(occurrences),
               source_host_semantics='source occurrence, not an independent predicted host or subject',
               count_semantics='unique exact matches only; zero exact mappings is not biological absence',
               annotation_transfer=False, scientifically_calibrated=False,
               evidence_receipt_sha256={tool: sha(path / 'status.json') for tool, path in tool_inputs.items()},
               output_files={p.name: sha(p) for p in out.iterdir() if p.is_file()})
    (out / 'summary_manifest.json').write_text(json.dumps(doc, indent=2) + '\n')
    return doc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', required=True)
    parser.add_argument('--evidence-dirs', nargs='*', default=[])
    parser.add_argument('--enabled-tools', default='')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    summarize(args.catalog, args.evidence_dirs, args.enabled_tools.split(',') if args.enabled_tools else [], args.outdir)
