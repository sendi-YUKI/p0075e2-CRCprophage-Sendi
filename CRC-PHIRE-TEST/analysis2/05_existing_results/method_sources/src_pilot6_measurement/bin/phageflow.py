#!/usr/bin/env python3
"""Deterministic, dependency-free data contracts for the A0/A1 phage workflow.

geNomad schema: https://portal.nersc.gov/genomad/quickstart.html
CheckV schema: actual pinned CheckV 1.0.3 quality_summary.tsv (14 columns).
The separate complete_genomes.tsv is preserved as raw output; its optional
predictions are not invented as an absent quality_summary column.
Coordinates in exported tables are zero-based, half-open. The canonical sequence
uses the source contig's orientation; caller orientation is recorded separately.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

VERSION = "0.1.0"
BOUNDARY_VERSION = "genomad_raw_v1"
DNA = frozenset("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")
GENOMAD_FIELDS = ["seq_name", "length", "topology", "coordinates", "n_genes", "genetic_code", "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy"]
MAP_FIELDS = ["genome_id", "genome_namespace", "contig_id", "original_contig_id", "original_header", "length", "sequence_sha256"]
STATS_FIELDS = ["genome_id", "length", "contigs", "N50", "gc_percent", "gc_denominator_bases", "ambiguous_bases", "gc_percent_total_length", "input_sha256", "normalized_fasta_sha256", "genome_qc_status"]
HIT_FIELDS = ["candidate_id", "genome_id", "source_host_genome", "genome_namespace", "contig_id", "original_contig_id", "candidate_type", "start0", "end0", "strand", "caller_strand", "length", "contig_length", "touches_contig_start", "touches_contig_end", "boundary_version", "parent_candidate_id", "sequence_sha256", "caller_sequence_sha256", "caller", "caller_seq_name", "caller_topology", "caller_coordinates", "virus_score", "fdr", "n_genes", "n_hallmarks", "marker_enrichment", "genetic_code", "taxonomy", "viral_identity_confidence", "integration_evidence", "boundary_confidence", "genome_qc_status"]
COORD_FIELDS = ["candidate_id", "genome_id", "contig_id", "original_contig_id", "start0", "end0", "strand", "caller_strand", "length", "boundary_version", "parent_candidate_id", "sequence_sha256", "caller_sequence_sha256", "caller_seq_name", "caller_coordinates", "coordinate_system", "slice_verified"]
STATUS_FIELDS = ["genome_id", "sample_status", "genome_qc_status", "genomad_status", "candidate_count", "checkv_status", "error"]
CHECKV_RAW_FIELDS = ["contig_id", "contig_length", "provirus", "proviral_length", "gene_count", "viral_genes", "host_genes", "checkv_quality", "miuvig_quality", "completeness", "completeness_method", "contamination", "kmer_freq", "warnings"]
CHECKV_FIELDS = ["candidate_id", "genome_id"] + ["checkv_" + k if not k.startswith("checkv_") else k for k in CHECKV_RAW_FIELDS if k != "contig_id"]
MASTER_FIELDS = HIT_FIELDS + [k for k in CHECKV_FIELDS if k not in HIT_FIELDS] + ["completeness", "completeness_method", "contamination", "checkv_status"]


class ContractError(ValueError):
    """An input/output contract was violated; never convert this to zero hits."""


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sequence_hash(sequence: str) -> str:
    return sha256_bytes(sequence.upper().encode("ascii"))


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def validate_id(value: str) -> str:
    if not ID_PATTERN.fullmatch(value):
        raise ContractError(f"Invalid genome_id {value!r}; use 1-128 ASCII letters/digits/._- beginning with a letter or digit")
    return value


def read_tsv(path: Path, required=()) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ContractError(f"Missing required table: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields) or any(not f for f in fields):
            raise ContractError(f"Empty or duplicate table header: {path}")
        if set(required) - set(fields):
            raise ContractError(f"Missing columns in {path}: {sorted(set(required) - set(fields))}")
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ContractError(f"Wrong number of columns at {path}:{line}")
            rows.append(row)
    return fields, rows


def write_tsv(path: Path, fields: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def read_fasta(path: Path, allow_empty: bool = False, *, max_bases: int | None = None) -> dict[str, tuple[str, str]]:
    """Return id -> (full original header, uppercase DNA); gzip magic is detected."""
    if not path.is_file():
        raise ContractError(f"Missing FASTA: {path}")
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    records: dict[str, tuple[str, str]] = {}
    header = None
    parts: list[str] = []
    total_bases = 0

    def finish():
        if header is None:
            return
        seq_id = header.split()[0]
        if seq_id in records:
            raise ContractError(f"Duplicate contig ID {seq_id!r} in {path}")
        sequence = "".join(parts).upper()
        if not sequence:
            raise ContractError(f"Empty sequence for {seq_id!r} in {path}")
        invalid = set(sequence) - DNA
        if invalid:
            raise ContractError(f"Invalid DNA characters {sorted(invalid)!r} for {seq_id!r} in {path}")
        records[seq_id] = (header, sequence)

    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                if not header or any(ord(c) < 32 for c in header):
                    raise ContractError(f"Invalid FASTA header at {path}:{line_number}")
                parts = []
            else:
                if header is None:
                    raise ContractError(f"Sequence before first FASTA header at {path}:{line_number}")
                sequence_line = "".join(line.split())
                total_bases += len(sequence_line)
                if max_bases is not None and total_bases > max_bases:
                    raise ContractError(f"FASTA exceeds inspection base budget: {max_bases}")
                parts.append(sequence_line)
    finish()
    if not records and not allow_empty:
        raise ContractError(f"Empty FASTA: {path}")
    return records


def write_fasta(path: Path, records) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def unique_index(rows, key: str, context: str):
    result = {}
    for row in rows:
        value = row[key]
        if not value or value in result:
            raise ContractError(f"Empty or duplicate {key} {value!r} in {context}")
        result[value] = row
    return result


def load_samplesheet(path: Path) -> list[dict[str, str]]:
    """Read/resolve the production samplesheet contract without writing outputs."""
    source = Path(path).expanduser().resolve()
    _, rows = read_tsv(source, ["genome_id", "fasta"])
    if not rows:
        raise ContractError("The samplesheet has no genomes")
    seen_ids = set()
    seen_files = set()
    normalized = []
    for row in rows:
        genome_id = validate_id(row["genome_id"])
        if genome_id in seen_ids:
            raise ContractError(f"Duplicate genome_id {genome_id!r}")
        seen_ids.add(genome_id)
        if not row["fasta"] or row["fasta"] != row["fasta"].strip():
            raise ContractError(f"Empty or whitespace-padded FASTA path for {genome_id}")
        path = Path(row["fasta"]).expanduser()
        if not path.is_absolute():
            path = source.parent / path
        path = path.resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ContractError(f"Missing or empty FASTA for {genome_id}: {path}")
        if path in seen_files:
            raise ContractError(f"The same FASTA file is assigned to multiple genome IDs: {path}")
        seen_files.add(path)
        # Content validation belongs to the prepare module, avoiding a second
        # full input read here while still failing before geNomad for bad DNA.
        normalized.append({"genome_id": genome_id, "fasta": str(path)})
    return sorted(normalized, key=lambda r: r["genome_id"])


def validate_sheet(args) -> None:
    source = Path(args.input).expanduser().resolve()
    normalized = load_samplesheet(source)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "normalized_samples.tsv", ["genome_id", "fasta"], normalized)
    write_json(out / "samplesheet_manifest.json", {"samplesheet": str(source), "sha256": file_hash(source), "sample_count": len(normalized), "relative_paths_resolved_against": str(source.parent)})


def prepare(args) -> None:
    genome_id = validate_id(args.genome_id)
    source = Path(args.fasta).expanduser().resolve()
    records = read_fasta(source)
    namespace = "g" + sha256_bytes(genome_id.encode("utf-8"))[:20]
    mapping = []
    canonical = []
    for original_id, (header, sequence) in sorted(records.items()):
        contig_id = namespace + "_c" + sha256_bytes(original_id.encode("utf-8"))[:24]
        mapping.append(dict(genome_id=genome_id, genome_namespace=namespace, contig_id=contig_id, original_contig_id=original_id, original_header=header, length=len(sequence), sequence_sha256=sequence_hash(sequence)))
        canonical.append((contig_id, sequence))
    unique_index(mapping, "contig_id", "namespaced contigs")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    write_fasta(out / "genome.fna", canonical)
    write_tsv(out / "contig_map.tsv", MAP_FIELDS, mapping)
    lengths = sorted((len(seq) for _, seq in canonical), reverse=True)
    total_length = sum(lengths)
    cumulative = 0
    for n50 in lengths:
        cumulative += n50
        if cumulative * 2 >= total_length:
            break
    bases = Counter("".join(seq for _, seq in canonical))
    called = sum(bases[b] for b in "ACGT")
    gc = bases["G"] + bases["C"]
    acquired_hash = file_hash(source)
    normalized_hash = file_hash(out / "genome.fna")
    stats = dict(genome_id=genome_id, length=total_length, contigs=len(records), N50=n50, gc_percent=f"{100 * gc / called:.6f}" if called else "unknown", gc_denominator_bases=called, ambiguous_bases=total_length - called, gc_percent_total_length=f"{100 * gc / total_length:.6f}", input_sha256=acquired_hash, normalized_fasta_sha256=normalized_hash, genome_qc_status="not_assessed")
    write_tsv(out / "genome_stats.tsv", STATS_FIELDS, [stats])
    write_json(out / "input_manifest.json", {"genome_id": genome_id, "input_path": str(source), "input_bytes": source.stat().st_size, "input_sha256": acquired_hash, "normalized_fasta_sha256": normalized_hash, "normalization": "uppercase IUPAC DNA; source orientation; sorted original IDs; 80-column FASTA", "namespace_method": "g + sha256(genome_id)[0:20]; contig: namespace + _c + sha256(original_id)[0:24]", "gc_definition": "100*(G+C)/(A+C+G+T), ambiguous bases excluded", "genome_qc_status": "not_assessed", "phageflow_version": VERSION})


def copy_contract_files(source: Path, target: Path, names) -> None:
    for name in names:
        path = source / name
        destination = target / name
        if path.is_file() and path.resolve() != destination.resolve():
            shutil.copy2(path, destination)


def normalize(args) -> None:
    genome_id = validate_id(args.genome_id)
    source_fasta = Path(args.fasta)
    records = read_fasta(source_fasta)
    _, map_rows = read_tsv(Path(args.mapping), MAP_FIELDS)
    mapping = unique_index(map_rows, "contig_id", "contig_map.tsv")
    if set(mapping) != set(records):
        raise ContractError("Contig map and normalized genome IDs disagree")
    for contig, row in mapping.items():
        if row["genome_id"] != genome_id or int(row["length"]) != len(records[contig][1]) or row["sequence_sha256"] != sequence_hash(records[contig][1]):
            raise ContractError(f"Contig mapping integrity failure: {contig}")
    summaries = sorted(Path(args.genomad_dir).rglob("*_virus_summary.tsv"))
    if len(summaries) != 1:
        raise ContractError(f"Expected exactly one geNomad virus summary; found {len(summaries)} in {args.genomad_dir}")
    summary = summaries[0]
    _, caller_rows = read_tsv(summary, GENOMAD_FIELDS)
    caller_index = unique_index(caller_rows, "seq_name", str(summary))
    virus_fasta = summary.with_name(summary.name.replace("_virus_summary.tsv", "_virus.fna"))
    fasta_choices = [p for p in (virus_fasta, virus_fasta.with_suffix(".fna.gz")) if p.is_file()]
    if len(fasta_choices) != 1:
        raise ContractError(f"Expected one plain or gzipped geNomad virus FASTA; found {len(fasta_choices)}")
    virus_fasta = fasta_choices[0]
    viral_sequences = read_fasta(virus_fasta, allow_empty=True)
    if set(caller_index) != set(viral_sequences):
        raise ContractError("geNomad virus summary and virus FASTA IDs disagree; output may be incomplete")
    provirus_tables = sorted(Path(args.genomad_dir).rglob("*_provirus.tsv"))
    if len(provirus_tables) != 1:
        raise ContractError(f"Expected one authoritative geNomad provirus table; found {len(provirus_tables)}")
    _, provirus_rows = read_tsv(provirus_tables[0], ["seq_name", "source_seq", "start", "end", "length"])
    provirus_index = unique_index(provirus_rows, "seq_name", str(provirus_tables[0]))
    hits = []
    coordinates = []
    sequences = []
    for caller_name, raw in sorted(caller_index.items()):
        provirus_match = re.fullmatch(r"(.+)\|provirus_(\d+)_(\d+)", caller_name)
        is_provirus = raw["topology"] == "Provirus"
        if is_provirus != bool(provirus_match):
            raise ContractError(f"Inconsistent geNomad provirus topology/name for {caller_name}")
        if is_provirus:
            contig, start1, end1 = provirus_match.groups()
            authoritative = provirus_index.get(caller_name)
            if authoritative is None:
                raise ContractError(f"Retained provirus is absent from find-proviruses output: {caller_name}")
            if authoritative["source_seq"] != contig or int(authoritative["start"]) != int(start1) or int(authoritative["end"]) != int(end1) or int(authoritative["length"]) != int(raw["length"]):
                raise ContractError(f"Provirus source/coordinate fields disagree with summary for {caller_name}")
            if raw["coordinates"] != f"{int(start1)}-{int(end1)}":
                raise ContractError(f"geNomad name/coordinate disagreement for {caller_name}")
            start, end = int(start1) - 1, int(end1)
            candidate_type = "provirus_locus"
            integration = "caller_predicted_provirus_locus"
            boundary = "caller_only_not_independently_validated"
        else:
            contig = caller_name
            if raw["coordinates"] not in ("NA", "", "None"):
                raise ContractError(f"Unexpected coordinates for non-provirus {caller_name}")
            start, end = 0, len(records.get(contig, ("", ""))[1])
            candidate_type = "full_contig_virus" if raw["topology"] in ("DTR", "ITR", "No terminal repeats") else "uncertain"
            integration = "unknown"
            boundary = "full_contig_extent_integration_unknown" if candidate_type == "full_contig_virus" else "unknown"
        if contig not in records:
            raise ContractError(f"geNomad references an unknown source contig: {contig}")
        source_sequence = records[contig][1]
        if not 0 <= start < end <= len(source_sequence):
            raise ContractError(f"Invalid/circular-wrap coordinates for {caller_name}; v0.1 supports non-wrapping intervals only")
        canonical_sequence = source_sequence[start:end]
        caller_sequence = viral_sequences[caller_name][1]
        if int(raw["length"]) != end - start or len(caller_sequence) != end - start:
            raise ContractError(f"Length/coordinate disagreement for {caller_name}")
        if caller_sequence == canonical_sequence:
            caller_strand = "+"
        elif caller_sequence == reverse_complement(canonical_sequence):
            caller_strand = "-"
        else:
            raise ContractError(f"Sequence does not match source-genome slice on either strand for {caller_name}")
        m = mapping[contig]
        digest = sequence_hash(canonical_sequence)
        identity = [m["genome_namespace"], contig, start, end, BOUNDARY_VERSION, digest]
        candidate_id = "ph_" + m["genome_namespace"] + "_" + sha256_bytes(json.dumps(identity, separators=(",", ":")).encode("ascii"))[:32]
        hit = dict(candidate_id=candidate_id, genome_id=genome_id, source_host_genome=genome_id, genome_namespace=m["genome_namespace"], contig_id=contig, original_contig_id=m["original_contig_id"], candidate_type=candidate_type, start0=start, end0=end, strand="+", caller_strand=caller_strand, length=end-start, contig_length=len(source_sequence), touches_contig_start=str(start == 0).lower(), touches_contig_end=str(end == len(source_sequence)).lower(), boundary_version=BOUNDARY_VERSION, parent_candidate_id="", sequence_sha256=digest, caller_sequence_sha256=sequence_hash(caller_sequence), caller="geNomad", caller_seq_name=caller_name, caller_topology=raw["topology"], caller_coordinates=raw["coordinates"], viral_identity_confidence="geNomad_predicted_score_reported", integration_evidence=integration, boundary_confidence=boundary, genome_qc_status="not_assessed")
        hit.update({key: raw[key] for key in ("virus_score", "fdr", "n_genes", "n_hallmarks", "marker_enrichment", "genetic_code", "taxonomy")})
        hits.append(hit)
        coordinates.append({**hit, "coordinate_system": "0-based-half-open", "slice_verified": "true"})
        sequences.append((candidate_id, canonical_sequence))
    unique_index(hits, "candidate_id", "normalized candidates")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    hits.sort(key=lambda r: r["candidate_id"])
    write_tsv(out / "caller_hits.tsv", HIT_FIELDS, hits)
    write_tsv(out / "coordinate_mapping.tsv", COORD_FIELDS, sorted(coordinates, key=lambda r: r["candidate_id"]))
    write_fasta(out / "candidates.fna", sorted(sequences))
    (out / "candidate_count.txt").write_text(str(len(hits)) + "\n", encoding="ascii")
    status = dict(genome_id=genome_id, sample_status="candidates_pending_checkv" if hits else "zero_candidates", genome_qc_status="not_assessed", genomad_status="completed", candidate_count=len(hits), checkv_status="pending" if hits else "skipped_empty", error="")
    write_tsv(out / "sample_status.tsv", STATUS_FIELDS, [status])
    copy_contract_files(Path(args.mapping).parent, out, ["contig_map.tsv", "genome_stats.tsv", "input_manifest.json"])
    write_json(out / "normalization_manifest.json", {"genome_id": genome_id, "geNomad_summary_sha256": file_hash(summary), "geNomad_provirus_table_sha256": file_hash(provirus_tables[0]), "geNomad_virus_fasta_sha256": file_hash(virus_fasta), "source_normalized_fasta_sha256": file_hash(source_fasta), "candidate_fasta_sha256": file_hash(out / "candidates.fna"), "candidate_count": len(hits), "coordinate_system": "0-based-half-open", "orientation_policy": "source forward; caller orientation and caller sequence hash retained", "boundary_version": BOUNDARY_VERSION, "circular_origin_spanning": "explicitly_unsupported", "schema_source": "https://portal.nersc.gov/genomad/quickstart.html"})


def quality(args) -> None:
    candidates_path = Path(args.candidates)
    _, candidates = read_tsv(candidates_path, HIT_FIELDS)
    indexed = unique_index(candidates, "candidate_id", str(candidates_path))
    sequences = read_fasta(candidates_path.parent / "candidates.fna", allow_empty=True)
    if set(sequences) != set(indexed) or any(sequence_hash(seq) != indexed[candidate_id]["sequence_sha256"] for candidate_id, (_, seq) in sequences.items()):
        raise ContractError("Candidate FASTA and table IDs/hashes disagree before CheckV integration")
    _, statuses = read_tsv(candidates_path.parent / "sample_status.tsv", STATUS_FIELDS)
    if len(statuses) != 1 or statuses[0]["genomad_status"] != "completed" or statuses[0]["sample_status"] not in ("zero_candidates", "candidates_pending_checkv"):
        raise ContractError("Candidate status must represent one successfully completed geNomad sample")
    status = statuses[0]
    if int(status["candidate_count"]) != len(candidates) or any(r["genome_id"] != status["genome_id"] for r in candidates):
        raise ContractError("Candidate table and sample status disagree")
    quality_rows = []
    master = []
    if candidates:
        if not args.checkv:
            raise ContractError("CheckV output is required for nonempty candidates")
        quality_path = Path(args.checkv) / "quality_summary.tsv"
        _, raw_quality = read_tsv(quality_path, CHECKV_RAW_FIELDS)
        quality_index = unique_index(raw_quality, "contig_id", str(quality_path))
        if set(quality_index) != set(indexed):
            raise ContractError(f"CheckV/candidate IDs disagree; missing={sorted(set(indexed)-set(quality_index))}, unexpected={sorted(set(quality_index)-set(indexed))}")
        for candidate_id, candidate in sorted(indexed.items()):
            raw = quality_index[candidate_id]
            if int(raw["contig_length"]) != int(candidate["length"]):
                raise ContractError(f"CheckV input length differs for {candidate_id}")
            q = {"candidate_id": candidate_id, "genome_id": candidate["genome_id"]}
            q.update({key if key.startswith("checkv_") else "checkv_" + key: value for key, value in raw.items() if key != "contig_id"})
            quality_rows.append(q)
            master.append({**candidate, **q, "completeness": raw["completeness"] if raw["completeness"] not in ("", "NA") else "unknown", "completeness_method": raw["completeness_method"] if raw["completeness_method"] not in ("", "NA") else "unknown", "contamination": raw["contamination"] if raw["contamination"] not in ("", "NA") else "unknown", "checkv_status": "completed"})
        status.update(sample_status="completed_with_candidates", checkv_status="completed")
    else:
        # A pre-existing nonempty CheckV table on the empty branch is a stale
        # result, never silently reuse it as if it belonged to this run.
        if args.checkv and (Path(args.checkv) / "quality_summary.tsv").exists():
            _, stale = read_tsv(Path(args.checkv) / "quality_summary.tsv", CHECKV_RAW_FIELDS)
            if stale:
                raise ContractError("Nonempty CheckV output was supplied for zero candidates")
        status.update(sample_status="zero_candidates", checkv_status="skipped_empty")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    copy_contract_files(candidates_path.parent, out, ["caller_hits.tsv", "coordinate_mapping.tsv", "candidates.fna", "candidate_count.txt", "contig_map.tsv", "genome_stats.tsv", "input_manifest.json", "normalization_manifest.json"])
    write_tsv(out / "prophage_master.tsv", MASTER_FIELDS, master)
    write_tsv(out / "checkv_quality.tsv", CHECKV_FIELDS, quality_rows)
    write_tsv(out / "sample_status.tsv", STATUS_FIELDS, [status])
    write_json(out / "quality_manifest.json", {"genome_id": status["genome_id"], "candidate_count": len(candidates), "checkv_status": status["checkv_status"], "boundary_policy": "report_only_no_trimming_no_boundary_changes", "checkv_quality_summary_sha256": file_hash(Path(args.checkv) / "quality_summary.tsv") if candidates else None, "source_caller_hits_sha256": file_hash(candidates_path)})


def aggregate(args) -> None:
    inputs, out = Path(args.inputs), Path(args.outdir)
    if inputs.resolve() == out.resolve():
        raise ContractError("Aggregate output must be separate from its per-sample inputs")
    final_tables = sorted(p for p in inputs.rglob("prophage_master.tsv") if out.resolve() not in p.resolve().parents)
    if not final_tables:
        raise ContractError("No final per-sample prophage_master.tsv files found; this is not a zero-hit run")
    specs = {"prophage_master.tsv": MASTER_FIELDS, "caller_hits.tsv": HIT_FIELDS, "coordinate_mapping.tsv": COORD_FIELDS, "checkv_quality.tsv": CHECKV_FIELDS, "sample_status.tsv": STATUS_FIELDS, "contig_map.tsv": MAP_FIELDS, "genome_stats.tsv": STATS_FIELDS}
    combined = {name: [] for name in specs}
    fasta_records = []
    manifests = []
    sample_ids = set()
    for table in final_tables:
        sample_dir = table.parent
        _, statuses = read_tsv(sample_dir / "sample_status.tsv", STATUS_FIELDS)
        if len(statuses) != 1 or statuses[0]["sample_status"] not in ("completed_with_candidates", "zero_candidates"):
            raise ContractError(f"Incomplete or failed sample in aggregate: {sample_dir}")
        genome_id = statuses[0]["genome_id"]
        if genome_id in sample_ids:
            raise ContractError(f"Multiple final directories for sample {genome_id}")
        sample_ids.add(genome_id)
        sample_tables = {}
        for name, fields in specs.items():
            path = sample_dir / name
            if not path.is_file() and name in ("contig_map.tsv", "genome_stats.tsv"):
                matches = []
                for possible in inputs.rglob(name):
                    _, possible_rows = read_tsv(possible, fields)
                    if possible_rows and {r["genome_id"] for r in possible_rows} == {genome_id}:
                        matches.append((possible, possible_rows))
                if not matches:
                    raise ContractError(f"Missing {name} for sample {genome_id}")
                # Repeated published copies must have identical content.
                if len({file_hash(p) for p, _ in matches}) != 1:
                    raise ContractError(f"Conflicting {name} copies for {genome_id}")
                rows = matches[0][1]
            else:
                _, rows = read_tsv(path, fields)
            if any(r.get("genome_id") != genome_id for r in rows):
                raise ContractError(f"Mixed genome IDs in {path}")
            sample_tables[name] = rows
            combined[name].extend(rows)
        fasta = read_fasta(sample_dir / "candidates.fna", allow_empty=True)
        _, sample_master = read_tsv(table, MASTER_FIELDS)
        expected = unique_index(sample_master, "candidate_id", str(table))
        if set(fasta) != set(expected) or len(fasta) != int(statuses[0]["candidate_count"]):
            raise ContractError(f"Aggregate FASTA/table/status counts disagree for {genome_id}")
        for name in ("caller_hits.tsv", "coordinate_mapping.tsv", "checkv_quality.tsv"):
            index = unique_index(sample_tables[name], "candidate_id", str(sample_dir / name))
            if set(index) != set(expected):
                raise ContractError(f"Aggregate candidate IDs disagree between {name} and master for {genome_id}")
            if name != "checkv_quality.tsv":
                for candidate_id, row in index.items():
                    for key in ("contig_id", "start0", "end0", "strand", "length", "boundary_version", "sequence_sha256"):
                        if row[key] != expected[candidate_id][key]:
                            raise ContractError(f"Aggregate {key} differs between {name} and master for {candidate_id}")
            else:
                for candidate_id, row in index.items():
                    for key in CHECKV_FIELDS:
                        if row[key] != expected[candidate_id][key]:
                            raise ContractError(f"Aggregate CheckV field {key} differs from master for {candidate_id}")
        for candidate_id, (_, sequence) in fasta.items():
            if sequence_hash(sequence) != expected[candidate_id]["sequence_sha256"]:
                raise ContractError(f"Aggregate sequence hash mismatch for {candidate_id}")
            fasta_records.append((candidate_id, sequence))
        for manifest_name in ("input_manifest.json", "normalization_manifest.json", "quality_manifest.json"):
            path = sample_dir / manifest_name
            if not path.exists() and manifest_name == "input_manifest.json":
                choices = [p for p in inputs.rglob(manifest_name) if json.loads(p.read_text(encoding="utf-8"))["genome_id"] == genome_id]
                if not choices:
                    raise ContractError(f"Missing input manifest for {genome_id}")
                path = choices[0]
            if not path.exists():
                raise ContractError(f"Missing {manifest_name} for {genome_id}")
            manifest_content = json.loads(path.read_text(encoding="utf-8"))
            if manifest_content.get("genome_id") != genome_id:
                raise ContractError(f"Manifest genome identity disagrees for {path}")
            manifests.append({"genome_id": genome_id, "manifest_type": manifest_name, "content": manifest_content, "manifest_sha256": file_hash(path)})
    unique_index(combined["prophage_master.tsv"], "candidate_id", "aggregate master")
    out.mkdir(parents=True, exist_ok=True)
    for name, fields in specs.items():
        write_tsv(out / name, fields, sorted(combined[name], key=lambda r: (r["genome_id"], r.get("candidate_id", r.get("contig_id", "")))))
    write_fasta(out / "candidates.fna", sorted(fasta_records))
    types = Counter(r["candidate_type"] for r in combined["prophage_master.tsv"])
    count_rows = [dict(stage="genomes", entity_type="genome", count=len(sample_ids)), dict(stage="input_contigs", entity_type="contig", count=len(combined["contig_map.tsv"])), dict(stage="genomad_candidates", entity_type="candidate", count=len(fasta_records)), dict(stage="coordinate_verified", entity_type="candidate", count=len(combined["coordinate_mapping.tsv"])), dict(stage="checkv_assessed", entity_type="candidate", count=len(combined["checkv_quality.tsv"])), dict(stage="zero_candidate_genomes", entity_type="genome", count=sum(s["sample_status"] == "zero_candidates" for s in combined["sample_status.tsv"]))]
    count_rows.extend(dict(stage=k, entity_type="candidate", count=types[k]) for k in ("provirus_locus", "full_contig_virus", "uncertain"))
    write_tsv(out / "stage_counts.tsv", ["stage", "entity_type", "count"], count_rows)
    write_json(out / "data_manifests.json", manifests)
    write_json(out / "aggregate_manifest.json", {"phageflow_version": VERSION, "samples": sorted(sample_ids), "files": {p.name: file_hash(p) for p in sorted(out.iterdir()) if p.is_file() and p.suffix in (".tsv", ".fna")}})
    rows_html = "".join("<tr>" + "".join(f"<td>{html.escape(str(row[k]))}</td>" for k in ("genome_id", "sample_status", "candidate_count", "checkv_status", "genome_qc_status")) + "</tr>" for row in sorted(combined["sample_status.tsv"], key=lambda r: r["genome_id"]))
    links = "".join(f'<li><a href="{name}">{name}</a></li>' for name in [*specs, "stage_counts.tsv", "candidates.fna", "data_manifests.json", "aggregate_manifest.json"])
    document = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>CRC phage pipeline A0+A1 report</title><style>body{{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#162a35}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #ccd5da;padding:10px;text-align:left}}.note{{background:#eef5f8;padding:18px}}</style><h1>CRC phage pipeline v{VERSION}</h1><p>{len(sample_ids)} genomes; {len(fasta_records)} viral candidates, including {types['provirus_locus']} caller-predicted provirus loci.</p><div class="note">This A0+A1 report records geNomad predictions and CheckV quality estimates. Bacterial genome QC is <b>not_assessed</b>. A full-contig virus is not evidence of integration. CheckV does not change the original candidate boundaries. Zero candidates means the caller completed without retained viral calls; it does not establish biological absence. These outputs do not establish caller sensitivity or specificity.</div><h2>Sample status</h2><table><tr><th>Genome</th><th>Status</th><th>Candidates</th><th>CheckV</th><th>Genome QC</th></tr>{rows_html}</table><h2>Auditable outputs</h2><ul>{links}</ul><p>Internal coordinates: zero-based, half-open. Candidate FASTA follows the source contig orientation; caller orientation is retained in the mapping. Raw outputs and execution manifests are retained by the Nextflow workflow.</p></html>'''
    (out / "report.html").write_text(document + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--version", action="version", version=VERSION)
    sub = result.add_subparsers(dest="command", required=True)
    val = sub.add_parser("validate-sheet")
    val.add_argument("--input", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--genome-id", required=True)
    prep.add_argument("--fasta", required=True)
    norm = sub.add_parser("normalize")
    norm.add_argument("--genome-id", required=True)
    norm.add_argument("--fasta", required=True)
    norm.add_argument("--mapping", required=True)
    norm.add_argument("--genomad-dir", required=True)
    qual = sub.add_parser("quality")
    qual.add_argument("--candidates", required=True)
    qual.add_argument("--checkv")
    agg = sub.add_parser("aggregate")
    agg.add_argument("--inputs", required=True)
    for p, fn in ((val, validate_sheet), (prep, prepare), (norm, normalize), (qual, quality), (agg, aggregate)):
        p.add_argument("--outdir", required=True)
        p.set_defaults(handler=fn)
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        args.handler(args)
    except (ContractError, OSError, UnicodeError, ValueError, csv.Error) as exc:
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        write_json(out / "failure.json", {"status": "failed", "stage": args.command, "error": str(exc), "phageflow_version": VERSION})
        genome_id = getattr(args, "genome_id", None)
        if args.command == "quality":
            try:
                _, previous = read_tsv(Path(args.candidates).parent / "sample_status.tsv", STATUS_FIELDS)
                genome_id = previous[0]["genome_id"] if len(previous) == 1 else None
            except (OSError, ContractError):
                pass
        if genome_id:
            write_tsv(out / "sample_status.tsv", STATUS_FIELDS, [dict(genome_id=genome_id, sample_status="failed", genome_qc_status="not_assessed", genomad_status="completed" if args.command == "quality" else "failed" if args.command == "normalize" else "not_run", candidate_count="unknown", checkv_status="failed" if args.command == "quality" else "not_run", error=str(exc))])
        print(f"ERROR [{args.command}]: {exc}", file=sys.stderr)
        return 2
    # A successful explicit rerun must not leave an earlier failure marker.
    (Path(args.outdir) / "failure.json").unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
