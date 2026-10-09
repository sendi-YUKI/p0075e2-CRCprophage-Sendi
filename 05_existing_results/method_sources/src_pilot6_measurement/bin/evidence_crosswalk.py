#!/usr/bin/env python3
"""VIBRANT 1.2.1 ORFs to viral and source canonical genes; no annotation transfer."""
from collections import Counter, defaultdict
import hashlib
import re
from pathlib import Path

from viral_evidence import fasta, table, write, unique_file, sha

FIELDS = ['viral_sequence_id', 'tool_scaffold_id', 'tool_gene_id', 'tool_start0',
          'tool_end0', 'tool_strand', 'tool_cds_sha256', 'tool_protein_sha256',
          'tool_partial', 'target_namespace', 'source_occurrence_id', 'source_host_genome',
          'mapped_start0', 'mapped_end0', 'mapped_strand', 'target_gene_set_id',
          'target_gene_id', 'target_start0', 'target_end0', 'target_strand',
          'target_partial', 'target_translation_table', 'relationship', 'overlap_bp',
          'protein_identical', 'evidence', 'annotation_transfer']
TOOL_FIELDS = ['viral_sequence_id', 'tool_scaffold_id', 'tool_gene_id', 'start0',
               'end0', 'strand', 'cds_sha256', 'protein_sha256', 'partial', 'coordinate_basis']
GENE_REQUIRED = ['gene_id', 'start0', 'end0', 'strand', 'cds_sha256', 'protein_sha256']
OCC_REQUIRED = ['source_occurrence_id', 'viral_sequence_id', 'assembly_id', 'contig_id',
                'start0', 'end0', 'strand', 'sequence_sha256', 'slice_sha256',
                'slice_verification', 'source_host_genome', 'predicted_host',
                'biological_sample_id', 'library_id']


def digest(sequence):
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTRYMKSWBDHVN', 'TGCAYRKMSWVHDBN'))[::-1]


def headers(path):
    out = {}
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                header = line[1:].strip()
                uid = header.split()[0]
                if uid in out:
                    raise ValueError('Duplicate FASTA header: ' + uid)
                out[uid] = header
    return out


def parse_tool_genes(root, expected, seqs):
    """Exported FAA coordinates are 1-based inclusive on query, NOT fragment-local."""
    faa = unique_file(root, '*.phages_combined.faa')
    ffn = unique_file(root, '*.phages_combined.ffn')
    proteins, cds = fasta(faa), fasta(ffn)
    if set(proteins) != set(cds):
        raise ValueError('VIBRANT protein/CDS IDs disagree')
    mother = table(unique_file(root, 'VIBRANT_annotations_*.tsv'), ('protein', 'scaffold'))
    by_gene = {r['protein']: r['scaffold'] for r in mother}
    if len(by_gene) != len(mother) or not set(proteins) <= set(by_gene):
        raise ValueError('VIBRANT exported ORF/mother annotation join failed')
    protein_headers, cds_headers = headers(faa), headers(ffn)
    genes = []
    for gid in sorted(proteins):
        # Exact upstream export format, including tab-delimited coordinates and strand.
        match = re.fullmatch(r'\S+\s+\((\d+)\.\.(\d+)\)\s+(-?1)(?:\s+.*)?', protein_headers[gid])
        dna_match = re.fullmatch(r'\S+\s+\((\d+)\.\.(\d+)\)(?:\s+.*)?', cds_headers[gid])
        if not match or not dna_match or match.groups()[:2] != dna_match.groups():
            raise ValueError('Unrecognized/mismatched VIBRANT ORF coordinates: ' + gid)
        start, end = int(match[1]) - 1, int(match[2])
        strand = '+' if match[3] == '1' else '-'
        scaffold = by_gene[gid]
        parents = [uid for uid in expected if scaffold == uid or
                   re.fullmatch(re.escape(uid) + r'_fragment_\d+', scaffold)]
        if len(parents) != 1:
            raise ValueError('Ambiguous tool scaffold parent: ' + scaffold)
        uid = parents[0]
        if not 0 <= start < end <= len(seqs[uid]):
            raise ValueError('Tool ORF outside query')
        sliced = seqs[uid][start:end]
        if strand == '-':
            sliced = reverse_complement(sliced)
        if sliced != cds[gid]:
            raise ValueError('Tool ORF CDS differs from parent slice: ' + gid)
        protein = proteins[gid].removesuffix('*')
        if not protein or '*' in protein or not re.fullmatch(r'[A-Z]+', protein):
            raise ValueError('Invalid tool protein')
        genes.append(dict(viral_sequence_id=uid, tool_scaffold_id=scaffold, tool_gene_id=gid,
                          start0=start, end0=end, strand=strand, cds_sha256=digest(sliced),
                          protein_sha256=digest(protein), partial='unknown_upstream_export',
                          coordinate_basis='VIBRANT_1.2.1_query_1based_inclusive_to_0based_halfopen'))
    return genes


def validate_occurrences(rows, seqs):
    seen = set()
    for row in rows:
        sid, oid = row['viral_sequence_id'], row['source_occurrence_id']
        if not oid or oid in seen or sid not in seqs:
            raise ValueError('Duplicate/orphan occurrence')
        seen.add(oid)
        start, end = int(row['start0']), int(row['end0'])
        if not 0 <= start < end or end - start != len(seqs[sid]) or row['strand'] not in ('+', '-'):
            raise ValueError('Invalid occurrence coordinate contract')
        if row['sequence_sha256'] != digest(seqs[sid]) or row['slice_sha256'] != digest(seqs[sid]):
            raise ValueError('Occurrence/catalog sequence identity mismatch')
        if not row['assembly_id'] or not row['contig_id']:
            raise ValueError('Missing source genome/contig key')
        if row.get('host_resolved_prophage') == 'true':
            if (row['source_host_genome'].lower() in ('', 'unknown', 'not_assessed', 'not_applicable') or
                    row.get('candidate_type') != 'provirus_locus' or row['slice_verification'] != 'verified_this_run' or
                    row.get('integration_evidence', 'unknown').lower() in ('', 'unknown', 'none', 'absent', 'false', 'not_assessed', 'not_applicable')):
                raise ValueError('Unsupported host-resolved claim')


def validate_genes(rows, seqs=None):
    seen = set()
    for row in rows:
        gid = row['gene_id']
        if not gid or gid in seen:
            raise ValueError('Duplicate/empty target gene ID')
        seen.add(gid)
        start, end = int(row['start0']), int(row['end0'])
        if not 0 <= start < end or row['strand'] not in ('+', '-'):
            raise ValueError('Invalid target gene coordinates')
        if any(not re.fullmatch('[a-f0-9]{64}', row[k]) for k in ('cds_sha256', 'protein_sha256')):
            raise ValueError('Invalid target gene sequence hash')
        if seqs is not None:
            sid = row['viral_sequence_id']
            if sid not in seqs or end > len(seqs[sid]):
                raise ValueError('Orphan/out-of-bounds viral gene')
            seq = seqs[sid][start:end]
            if row['strand'] == '-':
                seq = reverse_complement(seq)
            if digest(seq) != row['cds_sha256']:
                raise ValueError('Viral gene CDS identity mismatch')


def crosswalk(tool_genes, viral_genes, occurrences, canonical, seqs):
    validate_genes(viral_genes, seqs)
    validate_genes(canonical)
    validate_occurrences(occurrences, seqs)
    viral, source, occs = defaultdict(list), defaultdict(list), defaultdict(list)
    for gene in viral_genes:
        viral[gene['viral_sequence_id']].append(gene)
    for gene in canonical:
        source[(gene['genome_id'], gene['contig_id'])].append(gene)
    for occ in occurrences:
        occs[occ['viral_sequence_id']].append(occ)
    rows = []
    for tool in tool_genes:
        sid, ts, te, strand = tool['viral_sequence_id'], int(tool['start0']), int(tool['end0']), tool['strand']
        contexts = [('viral', None, ts, te, strand, viral[sid], '')]
        for occ in occs[sid]:
            if occ['strand'] == '+':
                start, end, direction = int(occ['start0']) + ts, int(occ['start0']) + te, strand
            else:
                start, end = int(occ['end0']) - te, int(occ['end0']) - ts
                direction = '-' if strand == '+' else '+'
            reason = ''
            if occ['slice_verification'] != 'verified_this_run':
                reason = 'source_slice_not_verified'
            elif not canonical:
                reason = 'canonical_gene_set_not_supplied'
            targets = [] if reason else source[(occ['assembly_id'], occ['contig_id'])]
            contexts.append(('host_canonical', occ, start, end, direction, targets, reason))
        if not occs[sid]:
            contexts.append(('host_canonical', None, ts, te, strand, [], 'source_occurrence_not_supplied'))
        for namespace, occ, start, end, direction, targets, reason in contexts:
            candidates = [(g, min(end, int(g['end0'])) - max(start, int(g['start0']))) for g in targets
                          if g['strand'] == direction and min(end, int(g['end0'])) > max(start, int(g['start0']))]
            base = dict(viral_sequence_id=sid, tool_scaffold_id=tool['tool_scaffold_id'], tool_gene_id=tool['tool_gene_id'],
                        tool_start0=ts, tool_end0=te, tool_strand=strand, tool_cds_sha256=tool['cds_sha256'],
                        tool_protein_sha256=tool['protein_sha256'], tool_partial=tool['partial'], target_namespace=namespace,
                        source_occurrence_id=occ['source_occurrence_id'] if occ else '',
                        source_host_genome=occ['source_host_genome'] if occ else 'not_applicable',
                        mapped_start0=start, mapped_end0=end, mapped_strand=direction,
                        annotation_transfer='not_performed')
            if not candidates:
                rows.append(dict(base, target_gene_set_id='', target_gene_id='', target_start0='', target_end0='',
                                 target_strand='', target_partial='', target_translation_table='', relationship='unmapped',
                                 overlap_bp=0, protein_identical='unknown', evidence=reason or 'no_same_strand_coordinate_overlap'))
            for target, overlap in candidates:
                same_bounds = start == int(target['start0']) and end == int(target['end0'])
                same_cds = tool['cds_sha256'] == target['cds_sha256']
                if same_bounds and not same_cds:
                    raise ValueError('Exact coordinate target has conflicting CDS hash')
                identical = tool['protein_sha256'] == target['protein_sha256']
                relation = 'ambiguous' if len(candidates) > 1 else 'exact' if same_bounds and same_cds and identical else 'partial'
                rows.append(dict(base, target_gene_set_id=target.get('gene_set_id', 'v0.2_host_canonical'),
                                 target_gene_id=target['gene_id'], target_start0=target['start0'], target_end0=target['end0'],
                                 target_strand=target['strand'], target_partial=target.get('partial', 'unknown'),
                                 target_translation_table=target.get('translation_table', 'unknown'),
                                 relationship=relation, overlap_bp=overlap, protein_identical=str(identical).lower(),
                                 evidence='same_sequence_coordinates_CDS_and_protein' if relation == 'exact' else
                                 'same_source_same_strand_overlap_not_functional_equivalence'))
    reverse_degree = defaultdict(set)
    for row in rows:
        if row['target_gene_id']:
            reverse_degree[(row['target_namespace'], row['source_occurrence_id'], row['target_gene_set_id'], row['target_gene_id'])].add(row['tool_gene_id'])
    for row in rows:
        key = (row['target_namespace'], row['source_occurrence_id'], row['target_gene_set_id'], row['target_gene_id'])
        if row['target_gene_id'] and len(reverse_degree[key]) > 1:
            row.update(relationship='ambiguous', evidence='multiple_tool_ORFs_overlap_target_no_automatic_merge')
    return rows


def export(root, catalog, genes_dir, canonical_paths, tool_genes, regions, entries, seqs, out):
    if not genes_dir:
        raise ValueError('VIBRANT requires explicit independent viral gene directory')
    occurrences = table(catalog / 'source_occurrences.tsv', OCC_REQUIRED)
    viral = table(Path(genes_dir) / 'gene_table.tsv', GENE_REQUIRED + ['viral_sequence_id', 'gene_set_id'])
    canonical = [row for path in canonical_paths for row in table(path, GENE_REQUIRED + ['genome_id', 'contig_id'])]
    by_gene = {g['tool_gene_id']: g for g in tool_genes}
    if not {e['tool_gene_id'] for e in entries} <= set(by_gene):
        raise ValueError('AMG gene missing exported CDS/protein evidence')
    region_map = {r['tool_scaffold_id']: r for r in regions}
    for gene in tool_genes:
        region = region_map.get(gene['tool_scaffold_id'])
        if region is None or gene['viral_sequence_id'] != region['viral_sequence_id']:
            raise ValueError('Exported ORF has no matching viral region')
        # The published 1.2.1 fragment exporter slices [start-1:stop-1].
        allowance = 1 if gene['tool_scaffold_id'] != gene['viral_sequence_id'] else 0
        if int(gene['start0']) < int(region['start0']) or int(gene['end0']) > int(region['end0']) + allowance:
            raise ValueError('ORF is inconsistent with exported region bounds')
    rows = crosswalk(tool_genes, viral, occurrences, canonical, seqs)
    write(out / 'tool_genes.tsv', tool_genes, TOOL_FIELDS)
    write(out / 'gene_crosswalk.tsv', rows, FIELDS)
    amg = {e['tool_gene_id'] for e in entries}
    exact = [r for r in rows if r['relationship'] == 'exact' and r['tool_gene_id'] in amg]
    return dict(status='completed', tool_orf_count=len(tool_genes), amg_entry_count=len(entries),
                amg_unique_tool_genes=len(amg),
                amg_exact_unique_viral_genes=len({(r['target_gene_set_id'], r['target_gene_id']) for r in exact if r['target_namespace'] == 'viral'}),
                amg_exact_unique_host_canonical_genes=len({(r['target_gene_set_id'], r['target_gene_id']) for r in exact if r['target_namespace'] == 'host_canonical'}),
                relationship_rows=dict(Counter(r['relationship'] for r in rows)), annotation_transfer=False,
                count_semantics='unique genes per namespace; do not sum occurrence or annotation rows',
                viral_gene_table_sha256=sha(Path(genes_dir) / 'gene_table.tsv'),
                canonical_inputs=[dict(path=str(p), sha256=sha(p)) for p in canonical_paths])
