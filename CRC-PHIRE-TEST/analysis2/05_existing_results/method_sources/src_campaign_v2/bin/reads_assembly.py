#!/usr/bin/env python3
"""R2/R3 assemblers, assembly provenance, annotation validation and reports."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

from reads_io import dump, sha, run, write_tsv, fq_counts, ident, table


def fasta(path):
    key, seq, seen = None, [], set()
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if key is not None:
                    if not seq:
                        raise ValueError('Empty FASTA record')
                    yield key, ''.join(seq).upper()
                header = line[1:].split()
                if not header or header[0] in seen:
                    raise ValueError('Empty/duplicate FASTA header')
                key, seq = header[0], []
                seen.add(key)
            elif key is None:
                raise ValueError('Sequence before FASTA header')
            else:
                seq.append(line)
        if key is not None:
            if not seq:
                raise ValueError('Empty FASTA record')
            yield key, ''.join(seq).upper()


def normalize(source, target, namespace):
    ident(namespace, 'assembly unit ID')
    mapping = []
    with open(target, 'w') as out:
        for old, seq in fasta(source):
            if re.search('[^ACGTN]', seq):
                raise ValueError('Unsupported assembly alphabet')
            digest = hashlib.sha256(seq.encode()).hexdigest()
            name = namespace + '__' + digest[:20] + '__' + hashlib.sha256(old.encode()).hexdigest()[:8]
            out.write('>' + name + '\n' + seq + '\n')
            mapping.append(dict(contig_id=name, original_contig_id=old, length=len(seq), sequence_sha256=digest,
                                strand='+', start0=0, end0=len(seq)))
    return mapping


def verify_gbff(fna, gbff):
    from Bio import SeqIO
    seqs, ids = dict(fasta(fna)), set()
    with Path(gbff).open() as handle:
        records = list(SeqIO.parse(handle, 'genbank'))
    if not seqs or not records:
        raise ValueError('A matching annotation requires nonempty assembly and GBFF')
    for record in records:
        if record.id in ids or record.id not in seqs or str(record.seq).upper() != seqs[record.id]:
            raise ValueError('GBFF contig identity/sequence does not exactly match assembly')
        ids.add(record.id)
    if ids != set(seqs):
        raise ValueError('GBFF does not cover exactly the assembly contigs')
    return dict(status='matched', records=len(records), fasta_sha256=sha(fna), gbff_sha256=sha(gbff))


def staged_file(value, unit_parent):
    original = Path(value)
    choices = [original, Path(unit_parent) / original.name, Path.cwd() / original.name]
    for path in choices:
        if path.is_file():
            return str(path.resolve())
    raise ValueError('Required input missing after staging: ' + str(value))


def clean_inputs(unit, unitfile):
    u = dict(unit)
    translated = {}
    for filename, digest in u['checksums'].items():
        p = staged_file(u.get('staged_input_paths', {}).get(filename, filename), Path(unitfile).parent)
        if sha(p) != digest:
            raise ValueError('Clean FASTQ checksum mismatch')
        translated[filename] = p
    for key in ['reads_1', 'reads_2']:
        if u.get(key):
            if u[key] not in translated:
                raise ValueError('Clean read path lacks checksum')
            u[key] = translated[u[key]]
    u['singletons'] = [translated[p] for p in u.get('singletons', [])]
    actual = fq_counts(u['reads_1'], u.get('reads_2'))
    if any(actual[k] != int(u['counts'][k]) for k in ['fragments', 'reads', 'bases']):
        raise ValueError('Clean FASTQ counts disagree with manifest')
    u['checksums'] = {translated[k]: v for k, v in u['checksums'].items()}
    return u


def assembly_command(u, work, threads=6, memory_gib=10):
    if not 1 <= threads <= 20:
        raise ValueError('Assembly threads must be 1..20')
    if not 1 <= memory_gib <= 100:
        raise ValueError('Assembly memory must be 1..100 GiB')
    singles = [p for p in u.get('singletons', []) if fq_counts(p)['reads']]
    primary = int(u['counts']['fragments']) > 0
    r1, r2 = u['reads_1'], u.get('reads_2')
    paired = bool(primary and r2)
    single_reads = ([r1] if primary and not r2 else []) + singles
    if not paired and not single_reads:
        return None
    if u['material_type'] == 'bacterial_isolate':
        cmd = ['spades.py', '-t', str(threads), '-m', str(memory_gib), '--phred-offset', '33', '-o', str(work)]
        if paired:
            cmd += ['--isolate', '-1', r1, '-2', r2]
            for p in singles:
                cmd += ['--pe1-s', p]
        else:
            cmd += ['--careful']
            for p in single_reads:
                cmd += ['--s1', p]
        return 'SPAdes', cmd, work / 'contigs.fasta', work / 'assembly_graph_with_scaffolds.gfa', paired, single_reads
    if u['material_type'] in {'bulk_metagenome', 'vlp'}:
        cmd = ['megahit', '-t', str(threads), '-m', str(memory_gib * 1024**3), '--min-contig-len', '500', '-o', str(work)]
        if paired:
            cmd += ['-1', r1, '-2', r2]
        if single_reads:
            if any(',' in str(p) for p in single_reads):
                raise ValueError('MEGAHIT input filename contains reserved comma')
            cmd += ['-r', ','.join(single_reads)]
        return 'MEGAHIT', cmd, work / 'final.contigs.fa', None, paired, single_reads
    raise ValueError('Unsupported assembly material')


def assemble(unitfile, out, threads=6, memory_gib=10):
    u = clean_inputs(json.load(open(unitfile)), unitfile)
    out = Path(out).resolve()
    if any(Path(p).resolve().is_relative_to(out) for p in u['checksums']):
        raise ValueError('Assembly output overlaps raw/clean input path')
    out.mkdir(parents=True, exist_ok=True)
    selection = assembly_command(u, out / 'assembler', threads, memory_gib)
    if selection is None:
        dump(out / 'assembly_unit.json', dict(u, status='not_assessed', reason='zero_clean_reads',
             assembly_fasta=None, assembly_qc_status='not_assessed'))
        return
    tool, cmd, raw, graph, paired, singles = selection
    run(cmd, out / 'commands.log')
    if not raw.is_file():
        raise ValueError('Assembler did not produce contigs')
    rawcopy = out / 'raw_contigs.fna'; shutil.copyfile(raw, rawcopy)
    normalized = out / 'assembly.fna'; mapping = normalize(rawcopy, normalized, u['unit_id'])
    write_tsv(out / 'contig_map.tsv', mapping, ['contig_id', 'original_contig_id', 'length', 'sequence_sha256', 'strand', 'start0', 'end0'])
    if not mapping:
        dump(out / 'assembly_unit.json', dict(u, status='zero_contigs', reason='assembler_returned_no_contigs',
             assembly_fasta=str(normalized), sha256=sha(normalized), contig_map=str(out / 'contig_map.tsv'), assembly_qc_status='not_assessed'))
        return
    run(['quast.py', normalized, '--threads', str(threads), '--min-contig', '0', '--no-plots', '--no-html', '-o', out / 'quast'], out / 'commands.log')
    index = out / 'readback_index'
    run(['bowtie2-build', '--threads', str(threads), normalized, index], out / 'commands.log')
    cmd = ['bowtie2', '--very-sensitive', '-p', str(threads), '-x', index, '-S', out / 'readback.sam']
    if paired:
        cmd += ['-1', u['reads_1'], '-2', u['reads_2']]
    if singles:
        cmd += ['-U', ','.join(singles)]
    run(cmd, out / 'readback.log')
    run(['samtools', 'sort', '-@', '1', '-m', '512M', '-o', out / 'readback.bam', out / 'readback.sam'], out / 'commands.log')
    with open(out / 'readback_flagstat.tsv', 'w') as handle:
        subprocess.run(['samtools', 'flagstat', '-O', 'tsv', out / 'readback.bam'], stdout=handle, check=True)
    if graph and graph.exists():
        shutil.copyfile(graph, out / graph.name)
    result = {k: u.get(k) for k in ['unit_id', 'biological_sample_id', 'library_id', 'subject_id', 'material_type', 'preprocessing_id', 'source_readsets', 'metadata']}
    result.update(dict(status='success', assembler=tool, assembly_id=u['unit_id'],
                       input_kind='bacterial_genome' if u['material_type'] == 'bacterial_isolate' else 'virome_assembly' if u['material_type'] == 'vlp' else 'metagenome_assembly',
                       assembly_fasta=str(normalized), raw_assembly_fasta=str(rawcopy), raw_sha256=sha(rawcopy),
                       sha256=sha(normalized), contig_map=str(out / 'contig_map.tsv'), contig_map_sha256=sha(out / 'contig_map.tsv'),
                       source_clean_checksums=u['checksums'], assembly_qc_status='not_assessed',
                       assembly_qc_reason='QUAST/basic_readback_only;biological_QC_is_separate',
                       readback_input_policy='all_used_primary_reads_and_retained_singletons',
                       assembler_paired_reads_used=paired, assembler_single_read_files=singles, fastq_quality_encoding='Phred33_contract',
                       graph=str(out / graph.name) if graph and graph.exists() else None,
                       graph_status='retained' if graph and graph.exists() else 'not_produced_by_configured_assembler', genbank=None))
    dump(out / 'assembly_unit.json', result)


def annotate(unitfile, out, db, threads=6):
    u = json.loads(Path(unitfile).read_text()); out = Path(out).resolve()
    source = Path(unitfile).resolve().parent
    if out == source or out.is_relative_to(source) or source.is_relative_to(out):
        raise ValueError('Annotation output must not overlap assembly input')
    if db and (out.is_relative_to(Path(db).resolve()) or Path(db).resolve().is_relative_to(out)):
        raise ValueError('Annotation output overlaps database')
    out.mkdir(parents=True, exist_ok=True)
    if u['material_type'] != 'bacterial_isolate':
        raise ValueError('Bakta only applies to bacterial isolate assemblies')
    # A self-contained bundle avoids reliance on an earlier Nextflow work mount.
    bundle = out / 'assembly_source'
    shutil.copytree(source, bundle)
    for field in ['assembly_fasta', 'raw_assembly_fasta', 'contig_map', 'graph']:
        if u.get(field):
            located = Path(staged_file(u[field], source))
            try:
                rel = located.relative_to(source)
            except ValueError:
                rel = Path(located.name)
                shutil.copy2(located, bundle / rel)
            u[field] = str(bundle / rel)
    if u['status'] != 'success':
        u.update(annotation_status='not_assessed', annotation_reason=u.get('reason', u['status']), genbank=None)
        dump(out / 'assembly_unit.json', u)
        return
    if sha(u['assembly_fasta']) != u['sha256']:
        raise ValueError('Assembly changed before Bakta annotation')
    provided = (u.get('metadata') or {}).get('genbank')
    if provided:
        provided = staged_file(provided, source)
        expected = u['metadata'].get('genbank_sha256')
        if not expected or sha(provided) != expected:
            raise ValueError('Provided GBFF checksum missing or changed since input registration')
        match = verify_gbff(u['assembly_fasta'], provided)
        gb = out / 'annotation.gbff'; shutil.copyfile(provided, gb)
    else:
        if not db:
            raise ValueError('PhiSpy requested without matching GBFF: full Bakta database required')
        run(['bakta', '--db', db, '--threads', str(threads), '--keep-contig-headers', '--output', out / 'bakta',
             '--prefix', 'annotation', '--locus-tag', 'CRCP', u['assembly_fasta']], out / 'bakta.log')
        gb = out / 'bakta/annotation.gbff'
        match = verify_gbff(u['assembly_fasta'], gb)
    u.update(genbank=str(gb), genbank_sha256=sha(gb), annotation_validation=match, annotation_status='completed')
    dump(out / 'assembly_unit.json', u)


def collect(unitfiles, out):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    units, genomes, contigs, statuses, ids = [], [], [], [], set()
    for filename in unitfiles:
        u = json.load(open(filename))
        ident(u['unit_id'], 'assembly unit ID')
        if u['unit_id'] in ids:
            raise ValueError('Duplicate assembly unit ID')
        ids.add(u['unit_id']); units.append(u)
        statuses.append(dict(unit_id=u['unit_id'], biological_sample_id=u.get('biological_sample_id') or 'unknown', library_id=u.get('library_id') or 'unknown', material_type=u['material_type'], status=u['status'], reason=u.get('reason', '')))
        if u['status'] != 'success':
            continue
        if u['material_type'] == 'bacterial_isolate':
            genomes.append(dict(genome_id=u['unit_id'], fasta=u['assembly_fasta'], genbank=u.get('genbank') or ''))
        else:
            contigs.append(u)
    dump(out / 'assembly_manifest.json', dict(schema_version=1, kind='assembly_manifest', units=units))
    dump(out / 'viral_assembly_manifest.json', dict(schema_version=1, kind='assembly_manifest', units=contigs,
         source_statuses=statuses, excluded_non_success_units=[r['unit_id'] for r in statuses if r['status'] != 'success']))
    write_tsv(out / 'genomes.tsv', genomes, ['genome_id', 'fasta', 'genbank'])
    write_tsv(out / 'sample_status.tsv', statuses, ['unit_id', 'biological_sample_id', 'library_id', 'material_type', 'status', 'reason'])


def crosswalk(genes, gbff, out):
    from Bio import SeqIO
    ann, contigs, seen = [], set(), set()
    for rec in SeqIO.parse(gbff, 'genbank'):
        contigs.add(rec.id)
        for n, feature in enumerate(rec.features):
            if feature.type != 'CDS':
                continue
            gid = feature.qualifiers.get('locus_tag', [f'{rec.id}:CDS{n}'])[0]
            if gid in seen:
                raise ValueError('Duplicate GBFF CDS identifier')
            seen.add(gid)
            parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
            strand = '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else 'unknown'
            ann.append(dict(contig_id=rec.id, gene_id=gid, parts=parts, strand=strand, length=sum(b-a for a,b in parts)))
    canonical = table(genes); edges, counts_a, counts_b = [], {}, {}
    if len({g['gene_id'] for g in canonical}) != len(canonical):
        raise ValueError('Duplicate canonical gene ID')
    for a in canonical:
        if not {'start0', 'end0', 'strand', 'gene_id', 'contig_id'}.issubset(a):
            raise ValueError('Canonical schema requires explicit 0-based start0/end0; no guessed start conversion')
        aliases = [c for c in [a.get('contig_id'), a.get('original_contig_id')] if c in contigs]
        if len(set(aliases)) > 1:
            raise ValueError('Ambiguous canonical/GBFF contig mapping')
        cid = aliases[0] if aliases else None
        start, end = int(a['start0']), int(a['end0'])
        if not 0 <= start < end:
            raise ValueError('Invalid canonical coordinates')
        for b in ann:
            overlap = sum(max(0, min(end, y) - max(start, x)) for x,y in b['parts'])
            if cid == b['contig_id'] and a['strand'] == b['strand'] and overlap:
                edges.append(dict(canonical_gene_id=a['gene_id'], annotation_gene_id=b['gene_id'], overlap_bp=overlap,
                                  canonical_length=end-start, annotation_length=b['length'],
                                  exact=b['parts'] == [(start,end)], annotation_transfer='not_performed'))
                counts_a[a['gene_id']] = counts_a.get(a['gene_id'],0)+1
                counts_b[b['gene_id']] = counts_b.get(b['gene_id'],0)+1
    for edge in edges:
        na, nb = counts_a[edge['canonical_gene_id']], counts_b[edge['annotation_gene_id']]
        edge['relationship'] = 'complex' if na > 1 and nb > 1 else 'split' if na > 1 else 'merge' if nb > 1 else 'exact' if edge['exact'] else 'partial'
    for a in canonical:
        if a['gene_id'] not in counts_a:
            edges.append(dict(canonical_gene_id=a['gene_id'], annotation_gene_id='', relationship='unmatched', annotation_transfer='not_performed'))
    for b in ann:
        if b['gene_id'] not in counts_b:
            edges.append(dict(canonical_gene_id='', annotation_gene_id=b['gene_id'], relationship='unmatched', annotation_transfer='not_performed'))
    write_tsv(out, edges, ['canonical_gene_id','annotation_gene_id','relationship','overlap_bp','canonical_length','annotation_length','annotation_transfer'])


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='cmd', required=True)
    for name in ['assemble', 'annotate']:
        a = sub.add_parser(name); a.add_argument('--unit', required=True); a.add_argument('--outdir', required=True)
        a.add_argument('--threads', type=int, default=6, choices=range(1,21));a.add_argument('--memory-gib',type=int,default=10,choices=range(1,101))
        if name == 'annotate':
            a.add_argument('--db')
    a = sub.add_parser('collect'); a.add_argument('--units', nargs='*', default=[]); a.add_argument('--outdir', required=True)
    a = sub.add_parser('crosswalk'); a.add_argument('--genes', required=True); a.add_argument('--gbff', required=True); a.add_argument('--out', required=True)
    args = p.parse_args()
    if args.cmd == 'assemble':
        assemble(args.unit, args.outdir, args.threads, args.memory_gib)
    elif args.cmd == 'annotate':
        annotate(args.unit, args.outdir, args.db, args.threads)
    elif args.cmd == 'collect':
        collect(args.units, args.outdir)
    else:
        crosswalk(args.genes, args.gbff, args.out)


if __name__ == '__main__':
    main()
