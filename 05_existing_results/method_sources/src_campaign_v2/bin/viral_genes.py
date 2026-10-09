#!/usr/bin/env python3
"""Independent Pyrodigal-gv gene set and coordinate-based canonical crosswalk."""
import argparse
from collections import defaultdict
import importlib.metadata
import json
from pathlib import Path

import phageflow as pf
from viral_catalog import guard, stable, table, write

GENE_FIELDS = ['gene_set_id', 'viral_sequence_id', 'gene_id', 'contig_id', 'start0', 'end0', 'strand',
               'gene_order', 'partial', 'translation_table', 'gene_calling_version', 'gene_calling_mode',
               'model_description', 'protein_sha256', 'cds_sha256', 'genome_id', 'original_contig_id']
CROSS_FIELDS = ['source_occurrence_id', 'viral_sequence_id', 'gene_set_id', 'viral_gene_id',
                'canonical_gene_set_id', 'canonical_gene_id', 'relationship', 'overlap_bp',
                'source_gene_start0', 'source_gene_end0', 'source_gene_strand',
                'protein_identical', 'annotation_transfer', 'evidence']


def call(args):
    import pyrodigal_gv
    version = importlib.metadata.version('pyrodigal-gv')
    base_version = importlib.metadata.version('pyrodigal')
    lock = json.loads(Path(args.tools_manifest).read_text())
    if lock['tools'] != dict(pyrodigal_gv=version, pyrodigal=base_version):
        raise ValueError('Gene caller version differs from locked manifest')
    distribution = importlib.metadata.distribution('pyrodigal-gv')
    for item in lock['model_package_inventory']:
        p = Path(distribution.locate_file(item['path']))
        if not p.is_file() or pf.file_hash(p) != item['sha256']:
            raise ValueError('Changed packaged viral model: ' + item['path'])
    model_hash = lock['model_inventory_sha256']
    source = Path(args.catalog)
    out = guard(args.outdir, [source])
    receipt = out / 'gene_calling_status.json'
    pf.write_json(receipt, dict(status='running'))
    try:
        seqs = pf.read_fasta(source / 'viral_sequences.fna', allow_empty=True)
        rows, proteins, cds, statuses = [], [], [], []
        finder = pyrodigal_gv.ViralGeneFinder(meta=True)
        with (out / 'genes.gff').open('w') as gff:
            gff.write('##gff-version 3\n')
            for sid, (_, sequence) in sorted(seqs.items()):
                predictions = finder.find_genes(sequence.encode())
                description = str(predictions.metagenomic_bin.description) if predictions.metagenomic_bin is not None else 'unknown'
                translation_table = predictions.training_info.translation_table
                gsid = stable('vgs_', [sid, version, base_version, model_hash, 'meta_all_models', description, translation_table])
                ordered = sorted(predictions, key=lambda g: (g.begin, g.end, g.strand))
                for order, gene in enumerate(ordered, 1):
                    start, end, strand = gene.begin - 1, gene.end, '+' if gene.strand == 1 else '-'
                    if not 0 <= start < end <= len(sequence):
                        raise ValueError('Gene coordinates out of bounds')
                    seq = sequence[start:end]
                    if strand == '-':
                        seq = pf.reverse_complement(seq)
                    normalized = ''.join(c if c in 'ACGT' else 'N' for c in seq)
                    if gene.sequence() != normalized:
                        raise ValueError('Pyrodigal gene/source slice mismatch')
                    protein = gene.translate(include_stop=False)
                    gid = stable('vg_', [gsid, start, end, strand, gene.translation_table, pf.sequence_hash(seq), pf.sequence_hash(protein)])
                    row = dict(gene_set_id=gsid, viral_sequence_id=sid, gene_id=gid, contig_id=sid,
                               start0=start, end0=end, strand=strand, gene_order=order,
                               partial=f'{int(gene.partial_begin)}{int(gene.partial_end)}',
                               translation_table=gene.translation_table, gene_calling_version=version,
                               gene_calling_mode='Pyrodigal-gv_meta_all_models', model_description=description,
                               protein_sha256=pf.sequence_hash(protein), cds_sha256=pf.sequence_hash(seq),
                               genome_id=sid, original_contig_id=sid)
                    rows.append(row); proteins.append((gid, protein)); cds.append((gid, seq))
                    gff.write(f'{sid}\tPyrodigal-gv:{version}\tCDS\t{start+1}\t{end}\t.\t{strand}\t0\tID={gid};gene_set_id={gsid};translation_table={gene.translation_table};partial={row["partial"]}\n')
                statuses.append(dict(viral_sequence_id=sid, gene_set_id=gsid, gene_count=len(ordered),
                                     status='completed' if ordered else 'completed_zero_genes',
                                     translation_table=translation_table, model_description=description,
                                     alternative_code_decision='selected_by_viral_meta_model_score',
                                     source_sha256=pf.sequence_hash(sequence)))
        write(out / 'gene_table.tsv', rows, GENE_FIELDS)
        pf.write_fasta(out / 'genes.faa', proteins); pf.write_fasta(out / 'genes.ffn', cds)
        write(out / 'gene_set_status.tsv', statuses, ['viral_sequence_id', 'gene_set_id', 'gene_count', 'status',
              'translation_table', 'model_description', 'alternative_code_decision', 'source_sha256'])
        pf.write_json(receipt, dict(status='completed', pyrodigal_gv=version, pyrodigal=base_version,
                      gene_calling_mode='meta=True; viral_only=False; model selected separately per sequence',
                      input_sha256=pf.file_hash(source / 'viral_sequences.fna'),
                      n_sequences=len(seqs), n_genes=len(rows), real_data_validated=False, model_inventory_sha256=model_hash,
                      genome_id_semantics='viral_sequence_id computational adapter; never source bacterial genome',
                      model_integrity='container/package lock supplied through workflow tools manifest',
                      output_files={p.name: pf.file_hash(p) for p in out.iterdir() if p.is_file() and p != receipt}))
    except Exception as exc:
        pf.write_json(receipt, dict(status='failed', error=repr(exc)))
        raise


def gene_crosswalk(genes, occurrences, canonical):
    by_seq, by_source = defaultdict(list), defaultdict(list)
    for g in genes:
        by_seq[g['viral_sequence_id']].append(g)
    for g in canonical:
        by_source[(g['genome_id'], g['contig_id'])].append(g)
    output = []
    for occ in occurrences:
        old = [g for g in by_source[(occ['assembly_id'], occ['contig_id'])]
               if int(g['end0']) > int(occ['start0']) and int(g['start0']) < int(occ['end0'])]
        new = by_seq[occ['viral_sequence_id']]
        pairs, new_degree, old_degree = [], defaultdict(int), defaultdict(int)
        coordinates = {}
        for ng in new:
            if occ['strand'] == '+':
                start, end, strand = int(occ['start0']) + int(ng['start0']), int(occ['start0']) + int(ng['end0']), ng['strand']
            elif occ['strand'] == '-':
                start, end = int(occ['end0']) - int(ng['end0']), int(occ['end0']) - int(ng['start0'])
                strand = '-' if ng['strand'] == '+' else '+'
            else:
                raise ValueError('Unsupported source orientation')
            coordinates[ng['gene_id']] = start, end, strand
            for og in old:
                overlap = max(0, min(end, int(og['end0'])) - max(start, int(og['start0'])))
                if strand == og['strand'] and overlap:
                    pairs.append((ng, og, overlap))
                    new_degree[ng['gene_id']] += 1; old_degree[og['gene_id']] += 1
        for ng, og, overlap in pairs:
            start, end, strand = coordinates[ng['gene_id']]
            nd, od = new_degree[ng['gene_id']], old_degree[og['gene_id']]
            relation = ('complex_split_merge' if nd > 1 and od > 1 else 'merge' if nd > 1 else 'split' if od > 1
                        else 'exact' if start == int(og['start0']) and end == int(og['end0']) and ng['cds_sha256'] == og['cds_sha256'] else 'partial')
            output.append(dict(source_occurrence_id=occ['source_occurrence_id'], viral_sequence_id=occ['viral_sequence_id'],
                               gene_set_id=ng['gene_set_id'], viral_gene_id=ng['gene_id'],
                               canonical_gene_set_id=og.get('gene_set_id', 'v0.2_host_canonical'), canonical_gene_id=og['gene_id'],
                               relationship=relation, overlap_bp=overlap, source_gene_start0=start, source_gene_end0=end,
                               source_gene_strand=strand, protein_identical=str(ng['protein_sha256'] == og['protein_sha256']).lower(),
                               annotation_transfer='not_performed', evidence='same_source_contig_same_strand_coordinate_overlap'))
        for ng in new:
            if not new_degree[ng['gene_id']]:
                start, end, strand = coordinates[ng['gene_id']]
                output.append(dict(source_occurrence_id=occ['source_occurrence_id'], viral_sequence_id=occ['viral_sequence_id'],
                                   gene_set_id=ng['gene_set_id'], viral_gene_id=ng['gene_id'], canonical_gene_set_id='', canonical_gene_id='',
                                   relationship='unmatched', overlap_bp=0, source_gene_start0=start, source_gene_end0=end,
                                   source_gene_strand=strand, protein_identical='unknown', annotation_transfer='not_performed',
                                   evidence='no_canonical_gene_overlap' if canonical else 'canonical_gene_set_not_supplied'))
        for og in old:
            if not old_degree[og['gene_id']]:
                output.append(dict(source_occurrence_id=occ['source_occurrence_id'], viral_sequence_id=occ['viral_sequence_id'],
                                   gene_set_id='', viral_gene_id='', canonical_gene_set_id=og.get('gene_set_id', 'v0.2_host_canonical'),
                                   canonical_gene_id=og['gene_id'], relationship='unmatched', overlap_bp=0,
                                   source_gene_start0=og['start0'], source_gene_end0=og['end0'], source_gene_strand=og['strand'],
                                   protein_identical='unknown', annotation_transfer='not_performed', evidence='canonical_gene_not_recovered'))
    return output


def crosswalk(args):
    out = guard(args.outdir, [args.genes_dir, args.catalog] + args.canonical_genes)
    genes = table(Path(args.genes_dir) / 'gene_table.tsv', GENE_FIELDS)
    occ = table(Path(args.catalog) / 'source_occurrences.tsv')
    canonical = [row for p in args.canonical_genes for row in table(p, ['genome_id', 'contig_id', 'gene_id', 'start0', 'end0', 'strand', 'cds_sha256', 'protein_sha256'])]
    if len({r['gene_id'] for r in canonical}) != len(canonical):
        raise ValueError('Duplicate canonical gene IDs')
    write(out / 'gene_set_crosswalk.tsv', gene_crosswalk(genes, occ, canonical), CROSS_FIELDS)
    pf.write_json(out / 'crosswalk_manifest.json', dict(status='completed', output_sha256=pf.file_hash(out / 'gene_set_crosswalk.tsv'), relationship_semantics='coordinate evidence; split/merge not verified functional equivalence',
                  annotation_transfer=False, canonical_gene_tables=[dict(path=p, sha256=pf.file_hash(Path(p))) for p in args.canonical_genes]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    r = sub.add_parser('call'); r.add_argument('--catalog', required=True); r.add_argument('--outdir', required=True); r.add_argument('--tools-manifest', required=True)
    r = sub.add_parser('crosswalk'); r.add_argument('--catalog', required=True); r.add_argument('--genes-dir', required=True)
    r.add_argument('--canonical-genes', nargs='*', default=[]); r.add_argument('--outdir', required=True)
    a = p.parse_args()
    {'call': call, 'crosswalk': crosswalk}[a.command](a)


if __name__ == '__main__':
    main()
