#!/usr/bin/env python3
"""Opt-in canonical-gene host core, alignment, independent masks and tree interfaces.

Concatenated coordinates are alignment coordinates, never physical chromosomes.
No CRC outcomes, function hits, reassigned ORFs, or kinship estimators are consumed.
"""
import argparse
from bisect import bisect_left
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import unquote

from mainline_common import fasta, new_output, read_table, run, safe_id, save, sha, unique, write_fasta, write_table

GENE = ['genome_id', 'gene_id', 'contig_id', 'original_contig_id', 'start0', 'end0', 'strand', 'cds_sha256', 'protein_sha256', 'partial']
MEMBER = ['group_id', 'panaroo_group', 'genome_id', 'gene_id', 'selected_core', 'reason', 'presence_fraction', 'copy_count']
BLOCK = ['group_id', 'genome_id', 'gene_id', 'group_start0', 'group_end0', 'concat_start0', 'concat_end0',
         'cds_start0', 'cds_end0', 'contig_id', 'source_start0', 'source_end0', 'strand']
PARTITION = ['group_id', 'panaroo_group', 'start0', 'end0', 'present_genomes', 'total_genomes']
MASK = ['evidence_id', 'evidence_type', 'source', 'genome_id', 'contig_id', 'start0', 'end0']


def config(path, require_ready=True):
    path = Path(path)
    if path.stat().st_size > 1024 * 1024:
        raise ValueError('Host configuration exceeds 1 MiB metadata limit')
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result: raise ValueError('Duplicate host configuration key: ' + key)
            result[key] = value
        return result
    def constant(value):
        raise ValueError('Nonfinite JSON value: ' + value)
    c = json.loads(path.read_text(), object_pairs_hook=pairs, parse_constant=constant)
    allowed = {'schema_version', 'status', 'scope', 'genome_ids', 'core_presence_fraction', 'panaroo', 'alignment', 'tree', 'mask', 'proposal_source'}
    if set(c) != allowed or c['schema_version'] != 1:
        raise ValueError('Unexpected host configuration keys/version')
    if c['status'] not in ('proposed', 'frozen', 'synthetic'):
        raise ValueError('Invalid host configuration status')
    if require_ready and c['status'] == 'proposed':
        raise ValueError('Host scientific choices are proposed; explicit study freeze required before execution')
    if c['scope'] != 'user_declared_within_species':
        raise ValueError('A within-species input declaration is required; pipeline does not select strains')
    ids = c['genome_ids']
    if not isinstance(ids, list) or len(ids) != len(set(ids)) or (require_ready and len(ids) < 2):
        raise ValueError('Declare at least two unique genome IDs')
    for gid in ids:
        safe_id(gid)
    if type(c['core_presence_fraction']) not in (float, int) or not 0 < c['core_presence_fraction'] <= 1:
        raise ValueError('Core presence must be a fraction in (0,1]')
    p = c['panaroo']
    if set(p) != {'clean_mode', 'identity', 'family_threshold', 'refind_mode', 'input_policy'}:
        raise ValueError('Unexpected Panaroo configuration')
    if p['clean_mode'] not in ('strict', 'moderate', 'sensitive') or p['refind_mode'] != 'off' or p['input_policy'] != 'complete_cds_only_with_exclusion_receipt':
        raise ValueError('Canonical-only route requires refind off and explicit input policy')
    if any(type(p[k]) not in (float, int) or not 0 < p[k] <= 1 for k in ('identity', 'family_threshold')):
        raise ValueError('Invalid Panaroo identity threshold')
    if c['alignment'] != 'mafft_auto_dna_all_sites':
        raise ValueError('Only full DNA alignments are wired; SNP-only ascertainment is not implemented')
    t = c['tree']
    if set(t) != {'model', 'seed', 'ufboot', 'partition_mode'} or not re.fullmatch(r'[A-Za-z0-9+_,.-]+', str(t['model'])) or t['partition_mode'] not in ('linked', 'unpartitioned'):
        raise ValueError('Invalid tree configuration/model')
    if type(t['seed']) is not int or t['seed'] <= 0 or type(t['ufboot']) is not int or (t['ufboot'] != 0 and t['ufboot'] < 1000):
        raise ValueError('Seed must be positive; UFBoot must be explicitly 0 or >=1000')
    m = c['mask']
    if set(m) != {'version', 'strategy', 'use_current_a1_prophage_master'} or type(m['use_current_a1_prophage_master']) is not bool:
        raise ValueError('Invalid mask declaration')
    if not isinstance(m['version'], str) or not m['version'] or m['strategy'] not in (None, 'whole_gene', 'alignment_columns_union'):
        raise ValueError('Invalid mask version/strategy')
    if require_ready and (not m['version'] or m['version'] == 'pending' or m['strategy'] not in ('whole_gene', 'alignment_columns_union')):
        raise ValueError('Declare versioned whole_gene or alignment_columns_union mask policy')
    return c


def adapt(a):
    c = config(a.config)
    if a.genome_id not in c['genome_ids']:
        raise ValueError('Genome is outside the declared host set')
    out = new_output(a.outdir)
    root = Path(a.genes_dir)
    genes = unique(read_table(root / 'gene_table.tsv', GENE), 'gene_id')
    cds, proteins = fasta(root / 'genes.ffn'), fasta(root / 'genes.faa')
    contigs = fasta(Path(a.prepared) / 'genome.fna')
    if set(genes) != set(cds) or set(genes) != set(proteins):
        raise ValueError('Canonical catalogs disagree')
    features = {}
    for line in (root / 'genes.gff').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) != 9 or f[2] != 'CDS' or f[7] != '0':
            raise ValueError('Incompatible canonical GFF')
        attrs = dict(v.split('=', 1) for v in f[8].split(';') if v)
        gid = unquote(attrs['ID'])
        if gid in features:
            raise ValueError('Duplicate GFF gene')
        features[gid] = f
    if set(features) != set(genes):
        raise ValueError('GFF/canonical gene IDs disagree')
    excluded, included = [], []
    reverse = str.maketrans('ACGTN', 'TGCAN')
    for gid, row in sorted(genes.items()):
        f = features[gid]
        start, end = int(row['start0']), int(row['end0'])
        if row['genome_id'] != a.genome_id or row['contig_id'] not in contigs or not 0 <= start < end <= len(contigs[row['contig_id']]):
            raise ValueError('Gene source identity/coordinates disagree')
        if (f[0], int(f[3]) - 1, int(f[4]), f[6]) != (row['contig_id'], start, end, row['strand']):
            raise ValueError('GFF coordinates differ from canonical gene table')
        sequence = contigs[row['contig_id']][start:end]
        if row['strand'] == '-':
            sequence = sequence.translate(reverse)[::-1]
        elif row['strand'] != '+':
            raise ValueError('Unknown gene strand')
        if sequence != cds[gid] or hashlib.sha256(cds[gid].encode()).hexdigest() != row['cds_sha256'] or hashlib.sha256(proteins[gid].encode()).hexdigest() != row['protein_sha256']:
            raise ValueError('Canonical sequence/source checksum mismatch')
        reasons = []
        if row['partial'] != '00': reasons.append('partial_cds')
        if len(sequence) < 34 or len(sequence) % 3: reasons.append('panaroo_length_constraint')
        if '*' in proteins[gid] or set(sequence) - set('ACGT'): reasons.append('ambiguous_or_stop')
        if reasons:
            excluded.append(dict(genome_id=a.genome_id, gene_id=gid, reason=';'.join(reasons)))
        else:
            included.append(row)
    if not included:
        raise ValueError('No compatible CDS; not a biological absence result')
    with (out / (safe_id(a.genome_id) + '.gff')).open('w') as f:
        f.write('##gff-version 3\n')
        for row in included:
            f.write('\t'.join(features[row['gene_id']]) + '\n')
        f.write('##FASTA\n')
        for name, seq in contigs.items():
            f.write('>' + name + '\n' + seq + '\n')
    write_table(out / 'canonical_genes.tsv', included, GENE)
    write_table(out / 'excluded_genes.tsv', excluded, ['genome_id', 'gene_id', 'reason'])
    write_fasta(out / 'canonical.ffn', {r['gene_id']: cds[r['gene_id']] for r in included})
    write_table(out / 'source_contigs.tsv', [dict(genome_id=a.genome_id, contig_id=k, length=len(v)) for k, v in contigs.items()], ['genome_id', 'contig_id', 'length'])
    save(out / 'adapter.json', dict(genome_id=a.genome_id, input_genes=len(genes), compatible_genes=len(included),
        excluded_genes=len(excluded), source_gff_sha256=sha(root / 'genes.gff'), source_fna_sha256=sha(Path(a.prepared) / 'genome.fna'),
        source_gene_table_sha256=sha(root / 'gene_table.tsv'), config_sha256=sha(a.config),
        attributes='Existing IDs/coordinates/attributes preserved; no invented gene/product annotation', new_orfs=False))


def cluster(a):
    c = config(a.config)
    out = new_output(a.outdir)
    genes, seqs, gffs, exclusions, source_contigs = [], {}, [], [], []
    gids = []
    for d in a.inputs:
        p = Path(d)
        receipt = json.loads((p / 'adapter.json').read_text())
        gid = receipt['genome_id']
        if receipt['config_sha256'] != sha(a.config): raise ValueError('Adapter used another configuration')
        gids.append(gid)
        genes.extend(read_table(p / 'canonical_genes.tsv', GENE))
        for name, seq in fasta(p / 'canonical.ffn').items():
            if name in seqs: raise ValueError('Duplicate canonical gene across genomes')
            seqs[name] = seq
        gffs.append(p / (gid + '.gff'))
        exclusions.extend(read_table(p / 'excluded_genes.tsv'))
        source_contigs.extend(read_table(p / 'source_contigs.tsv'))
    if sorted(gids) != sorted(c['genome_ids']): raise ValueError('Declared host set differs from actual inputs')
    qc = unique(read_table(a.qc, ['genome_id', 'eligibility']), 'genome_id')
    if any(g not in qc or qc[g]['eligibility'] != 'eligible' for g in gids):
        raise ValueError('Host set includes missing/ineligible genome QC')
    unique(genes, 'gene_id')
    p = c['panaroo']
    run(['panaroo', '-i', *sorted(gffs), '-o', out / 'panaroo', '--clean-mode', p['clean_mode'],
         '--refind-mode', 'off', '-c', str(p['identity']), '-f', str(p['family_threshold']), '-t', str(a.cpus)], out, 'panaroo', timeout=a.timeout)
    write_table(out / 'canonical_genes.tsv', genes, GENE)
    write_table(out / 'excluded_genes.tsv', exclusions, ['genome_id', 'gene_id', 'reason'])
    write_table(out / 'source_contigs.tsv', source_contigs, ['genome_id', 'contig_id', 'length'])
    write_fasta(out / 'canonical.ffn', seqs)
    shutil.copyfile(a.config, out / 'host_config.json')
    save(out / 'cluster_status.json', dict(status='completed', genome_ids=sorted(gids), config_sha256=sha(a.config),
         new_orfs=False, refind_mode='off', qc_sha256=sha(a.qc), exclusions=len(exclusions)))


def alignment(a):
    c = config(a.config)
    out = new_output(a.outdir)
    src = Path(a.clusters)
    genes = unique(read_table(src / 'canonical_genes.tsv', GENE), 'gene_id')
    cds = fasta(src / 'canonical.ffn')
    genomes = sorted(c['genome_ids'])
    with (src / 'panaroo/gene_presence_absence_roary.csv').open(newline='') as f:
        r = csv.DictReader(f)
        if not set(genomes).issubset(r.fieldnames or []): raise ValueError('Panaroo output genome IDs disagree')
        groups = list(r)
    members, partitions, blocks = [], [], []
    combined = {g: '' for g in genomes}
    seen = set()
    for row in sorted(groups, key=lambda r: r['Gene']):
        by_genome = {g: row[g].split(';') if row[g] else [] for g in genomes}
        group = 'core_' + hashlib.sha256(json.dumps(sorted(gid for ids in by_genome.values() for gid in ids), separators=(',', ':')).encode()).hexdigest()[:20]
        count = sum(bool(v) for v in by_genome.values())
        selected = count / len(genomes) >= c['core_presence_fraction'] and all(len(v) <= 1 for v in by_genome.values())
        reason = 'selected_single_copy' if selected else 'paralog_present' if any(len(v) > 1 for v in by_genome.values()) else 'below_presence'
        for genome, ids in by_genome.items():
            for gid in ids:
                if gid not in genes or genes[gid]['genome_id'] != genome or gid in seen:
                    raise ValueError('Panaroo member is unknown, reassigned, or duplicated: ' + gid)
                seen.add(gid)
                members.append(dict(group_id=group, panaroo_group=row['Gene'], genome_id=genome, gene_id=gid,
                    selected_core=str(selected).lower(), reason=reason, presence_fraction=count / len(genomes), copy_count=len(ids)))
        if not selected: continue
        folder = out / 'genes' / group
        folder.mkdir(parents=True)
        source = {g: cds[ids[0]] for g, ids in by_genome.items() if ids}
        write_fasta(folder / 'input.ffn', source)
        if len(source) > 1:
            run(['mafft', '--thread', str(a.cpus), '--auto', '--nuc', folder / 'input.ffn'], folder, 'mafft', stdout=folder / 'aligned.fna', timeout=a.timeout)
        else: shutil.copyfile(folder / 'input.ffn', folder / 'aligned.fna')
        aligned = fasta(folder / 'aligned.fna')
        lengths = {len(s) for s in aligned.values()}
        if set(aligned) != set(source) or len(lengths) != 1: raise ValueError('Alignment identities/lengths disagree')
        size = lengths.pop()
        start = len(combined[genomes[0]])
        partitions.append(dict(group_id=group, panaroo_group=row['Gene'], start0=start, end0=start + size, present_genomes=count, total_genomes=len(genomes)))
        for genome in genomes:
            seq = aligned.get(genome, '-' * size)
            combined[genome] += seq
            if genome not in source: continue
            if seq.replace('-', '') != source[genome]: raise ValueError('MAFFT changed source sequence')
            gid = by_genome[genome][0]
            gene = genes[gid]
            cds_pos, i = 0, 0
            while i < size:
                if seq[i] == '-': i += 1; continue
                j = i
                while j < size and seq[j] != '-': j += 1
                n = j - i
                lo = int(gene['start0']) + cds_pos if gene['strand'] == '+' else int(gene['end0']) - cds_pos - n
                blocks.append(dict(group_id=group, genome_id=genome, gene_id=gid, group_start0=i, group_end0=j,
                    concat_start0=start + i, concat_end0=start + j, cds_start0=cds_pos, cds_end0=cds_pos + n,
                    contig_id=gene['contig_id'], source_start0=lo, source_end0=lo + n, strand=gene['strand']))
                cds_pos += n
                i = j
    write_table(out / 'core_members.tsv', members, MEMBER)
    write_table(out / 'partitions.tsv', partitions, PARTITION)
    write_table(out / 'alignment_blocks.tsv', blocks, BLOCK)
    write_table(out / 'unclustered_genes.tsv', [dict(gene_id=g, reason='not_retained_by_panaroo') for g in sorted(set(genes) - seen)], ['gene_id', 'reason'])
    shutil.copyfile(src / 'excluded_genes.tsv', out / 'excluded_genes.tsv')
    shutil.copyfile(src / 'source_contigs.tsv', out / 'source_contigs.tsv')
    write_fasta(out / 'alignment.fna', combined if partitions else {})
    (out / 'partitions.nex').write_text('#nexus\nbegin sets;\n' + ''.join(f"charset {p['group_id']} = {p['start0'] + 1}-{p['end0']};\n" for p in partitions) + 'end;\n')
    write_table(out / 'tree_tips.tsv', [dict(tip_id=g, genome_id=g, source_host_genome=g) for g in genomes], ['tip_id', 'genome_id', 'source_host_genome'])
    save(out / 'alignment_status.json', dict(status='completed' if partitions else 'no_eligible_core', selected_groups=len(partitions),
         genome_ids=genomes, columns=len(next(iter(combined.values()))), config_sha256=sha(a.config),
         coordinates='Concatenated alignment; not physical chromosome; no Gubbins input claim'))


def tree(a):
    c = config(a.config)
    out = new_output(a.outdir)
    source = Path(a.alignment)
    seqs = fasta(source / 'alignment.fna')
    state = dict(stage=a.stage, config_sha256=sha(a.config), input_alignment_sha256=sha(source / 'alignment.fna'),
                 tree_spec=c['tree'], scientific_calibration=False)
    if len(seqs) < 4 or not seqs:
        state.update(status='insufficient_tree_input', reason='Need at least four genomes and nonempty alignment; no invented tree')
    else:
        t = c['tree']
        cmd = ['iqtree2', '-s', source / 'alignment.fna', '-st', 'DNA', '-m', t['model'], '-seed', str(t['seed']),
               '-nt', str(a.cpus), '-pre', out / 'host_tree']
        if t['partition_mode'] == 'linked': cmd += ['-p', source / 'partitions.nex']
        if t['ufboot']: cmd += ['-B', str(t['ufboot'])]
        run(cmd, out, 'iqtree2', timeout=a.timeout)
        treefile = out / 'host_tree.treefile'
        if not treefile.is_file(): raise ValueError('IQ-TREE tree output missing')
        # Simple IDs are enforced at configuration time; no quoted Newick tip labels are accepted.
        tips = re.findall(r'(?:\(|,)\s*([A-Za-z0-9_.-]+)\s*:', treefile.read_text())
        if sorted(tips) != sorted(seqs): raise ValueError('Tree tip IDs differ from declared genomes')
        state.update(status='completed', tree_sha256=sha(treefile))
    shutil.copyfile(source / 'tree_tips.tsv', out / 'tree_tips.tsv')
    save(out / 'tree_status.json', state)
    if a.stage == 'final':
        save(out / 'lineage_kinship_interface.json', dict(schema_version=1, final_tree_status=state['status'],
             tree='host_tree.treefile' if state['status'] == 'completed' else None, tree_tips='tree_tips.tsv',
             lineage_assignment='not_assessed; clade definition not frozen', kinship_status='method_pending',
             kinship_method=None, kinship_normalization=None, distance_matrix_is_not_kinship=True,
             alignment_sha256=state['input_alignment_sha256'], mask_manifest_sha256=sha(source / 'mask_manifest.json')))


def mask(a):
    c = config(a.config)
    out = new_output(a.outdir)
    source = Path(a.alignment)
    preliminary = json.loads((Path(a.preliminary_tree) / 'tree_status.json').read_text())
    if preliminary['input_alignment_sha256'] != sha(source / 'alignment.fna'):
        raise ValueError('Preliminary tree/alignment identity mismatch')
    seqs = fasta(source / 'alignment.fna')
    blocks = read_table(source / 'alignment_blocks.tsv', BLOCK)
    partitions = unique(read_table(source / 'partitions.tsv', PARTITION), 'group_id')
    intervals, sources = [], []
    if c['mask']['use_current_a1_prophage_master']:
        master = read_table(a.master, ['candidate_id', 'genome_id', 'contig_id', 'start0', 'end0'])
        sources.append(dict(kind='independent_A1_prophage_coordinates', sha256=sha(a.master)))
        intervals.extend(dict(evidence_id=r['candidate_id'], evidence_type='prophage', source='current_a1_prophage_master',
                         **{k: r[k] for k in ['genome_id', 'contig_id', 'start0', 'end0']}) for r in master if r['genome_id'] in c['genome_ids'])
    bundle = Path(a.mask_bundle)
    manifest = json.loads((bundle / 'manifest.json').read_text())
    if manifest.get('schema_version') != 1 or manifest.get('purpose') != 'independent_mobile_intervals':
        raise ValueError('Only declared independent mobile evidence is permitted')
    for item in manifest['files']:
        p = (bundle / item['path']).resolve()
        if not p.is_relative_to(bundle.resolve()) or sha(p) != item['sha256']:
            raise ValueError('Changed/unsafe mobile evidence')
        rows = read_table(p, MASK)
        if any(r['evidence_type'] not in ('prophage', 'mobile') for r in rows):
            raise ValueError('Function/outcome evidence is not a host mask input')
        intervals.extend(rows)
        sources.append(item)
    if not sources:
        raise ValueError('No assessed mobile evidence source; cannot claim a completed mask')
    index = defaultdict(list)
    contigs = {(r['genome_id'], r['contig_id']): int(r['length']) for r in read_table(source / 'source_contigs.tsv')}
    for row in intervals:
        start, end = int(row['start0']), int(row['end0'])
        if not 0 <= start < end <= contigs.get((row['genome_id'], row['contig_id']), 0) or row['genome_id'] not in c['genome_ids']:
            raise ValueError('Invalid mobile evidence coordinates/host')
        index[(row['genome_id'], row['contig_id'])].append(row)
    events, remove = [], set()
    strategy = c['mask']['strategy']
    for block in blocks:
        lo, hi = int(block['source_start0']), int(block['source_end0'])
        for interval in index[(block['genome_id'], block['contig_id'])]:
            x, y = max(lo, int(interval['start0'])), min(hi, int(interval['end0']))
            if x >= y: continue
            if strategy == 'whole_gene':
                p = partitions[block['group_id']]
                start, end = int(p['start0']), int(p['end0'])
            else:
                start = int(block['concat_start0']) + (x - lo if block['strand'] == '+' else hi - y)
                end = start + y - x
            remove.update(range(start, end))
            events.append(dict(group_id=block['group_id'], genome_id=block['genome_id'], gene_id=block['gene_id'],
                  evidence_id=interval['evidence_id'], evidence_type=interval['evidence_type'], source=interval['source'],
                  contig_id=block['contig_id'], source_start0=x, source_end0=y, strand=block['strand'],
                  removed_start0=start, removed_end0=end, strategy=strategy))
    lengths = {len(v) for v in seqs.values()}
    if len(lengths) > 1: raise ValueError('Unequal alignment lengths')
    size = next(iter(lengths), 0)
    keep = [i for i in range(size) if i not in remove]
    write_fasta(out / 'alignment.fna', {g: ''.join(s[i] for i in keep) for g, s in seqs.items()} if keep else {})
    # Run-length maps retain exact source columns without a genome-by-base giant table.
    maps = []
    for final, original in enumerate(keep):
        if maps and original == maps[-1]['original_end0']:
            maps[-1]['original_end0'] += 1; maps[-1]['final_end0'] += 1
        else: maps.append(dict(original_start0=original, original_end0=original + 1, final_start0=final, final_end0=final + 1))
    final_parts = []
    for group, p in partitions.items():
        start, end = bisect_left(keep, int(p['start0'])), bisect_left(keep, int(p['end0']))
        if end > start:
            final_parts.append(dict(p, start0=start, end0=end))
    write_table(out / 'partitions.tsv', final_parts, PARTITION)
    (out / 'partitions.nex').write_text('#nexus\nbegin sets;\n' + ''.join(f"charset {p['group_id']} = {p['start0'] + 1}-{p['end0']};\n" for p in final_parts) + 'end;\n')
    write_table(out / 'column_crosswalk.tsv', maps, ['original_start0', 'original_end0', 'final_start0', 'final_end0'])
    write_table(out / 'mask_events.tsv', events, ['group_id', 'genome_id', 'gene_id', 'evidence_id', 'evidence_type', 'source', 'contig_id',
              'source_start0', 'source_end0', 'strand', 'removed_start0', 'removed_end0', 'strategy'])
    write_table(out / 'mask_intervals.tsv', intervals, MASK)
    shutil.copyfile(source / 'tree_tips.tsv', out / 'tree_tips.tsv')
    save(out / 'mask_manifest.json', dict(schema_version=1, status='completed' if keep else 'no_core_after_mask',
        mask_version=c['mask']['version'], strategy=strategy, assessed_sources=sources, mask_bundle_sha256=sha(bundle / 'manifest.json'),
        evidence_rows=len(intervals), removed_columns=len(remove), original_columns=size, final_columns=len(keep),
        input_alignment_sha256=sha(source / 'alignment.fna'), preliminary_tree_status=preliminary['status'],
        config_sha256=sha(a.config), scientific_calibration=False, coordinates='alignment columns; source genome coordinates retained in event/block crosswalks'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    q = sub.add_parser('check-config'); q.add_argument('--config', required=True); q.add_argument('--allow-proposed', action='store_true')
    q = sub.add_parser('adapt')
    for name in ['genes-dir', 'prepared', 'genome-id', 'config', 'outdir']: q.add_argument('--' + name, required=True)
    q = sub.add_parser('cluster')
    q.add_argument('--inputs', nargs='+', required=True)
    for name in ['qc', 'config', 'outdir']: q.add_argument('--' + name, required=True)
    q.add_argument('--cpus', type=int, default=1); q.add_argument('--timeout', type=int)
    q = sub.add_parser('align')
    for name in ['clusters', 'config', 'outdir']: q.add_argument('--' + name, required=True)
    q.add_argument('--cpus', type=int, default=1); q.add_argument('--timeout', type=int)
    q = sub.add_parser('tree')
    for name in ['alignment', 'config', 'outdir']: q.add_argument('--' + name, required=True)
    q.add_argument('--stage', choices=['preliminary', 'final'], required=True)
    q.add_argument('--cpus', type=int, default=1); q.add_argument('--timeout', type=int)
    q = sub.add_parser('mask')
    for name in ['alignment', 'preliminary-tree', 'master', 'mask-bundle', 'config', 'outdir']: q.add_argument('--' + name, required=True)
    a = p.parse_args()
    if hasattr(a, 'cpus') and not 1 <= a.cpus <= 6: raise ValueError('CPU budget must be 1..6')
    if a.action == 'check-config':
        c = config(a.config, require_ready=not a.allow_proposed)
        print(json.dumps(dict(status=c['status'], analysis_started=False, genome_count=len(c['genome_ids']))))
    else: {'adapt': adapt, 'cluster': cluster, 'align': alignment, 'tree': tree, 'mask': mask}[a.action](a)


if __name__ == '__main__':
    try: main()
    except Exception as e:
        print('ERROR: ' + str(e), file=sys.stderr)
        raise SystemExit(1)
