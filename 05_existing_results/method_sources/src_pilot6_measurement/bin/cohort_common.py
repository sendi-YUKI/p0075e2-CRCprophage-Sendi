"""Shared strict B1 I/O; independent of biological software PATH."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import os
import importlib.util
import re
from pathlib import Path

class ContractError(ValueError):
    pass

def require(condition, message):
    if not condition:
        raise ContractError(message)

_IPC_REJECT = None


def _native_input(path):
    if os.environ.get('CRC_PHAGE_RUNTIME') == 'spark-native':
        global _IPC_REJECT
        if _IPC_REJECT is None:
            repo = os.environ.get('CRC_SPARK_REPO')
            helper = Path(repo) / 'bin/spark_fs_guard.py' if repo else Path(__file__).resolve().with_name('spark_fs_guard.py')
            spec = importlib.util.spec_from_file_location('_crc_cohort_fs_guard', helper)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            _IPC_REJECT = module.reject_ipc_storage
        _IPC_REJECT([path])


def read_json(path):
    _native_input(path)
    return json.loads(Path(path).read_text())

def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")

def sha256(path):
    _native_input(path)
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

def text_open(path):
    _native_input(path)
    return gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()

def fasta(path):
    records = {}
    name = None
    with text_open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                require(name and name not in records, f"Duplicate/empty FASTA ID: {name}")
                records[name] = ""
            else:
                require(name is not None, "FASTA sequence precedes header")
                require(re.fullmatch("[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line), "Invalid DNA FASTA")
                records[name] += line.upper()
    require(all(records.values()), "Empty FASTA sequence")
    return records

def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for pos in range(0, len(sequence), 80):
                handle.write(sequence[pos:pos + 80] + "\n")

def table(path):
    _native_input(path)
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"Missing TSV header: {path}")
        require(all(reader.fieldnames) and len(set(reader.fieldnames)) == len(reader.fieldnames),
                f"Duplicate/empty TSV header: {path}")
        rows = list(reader)
        require(all(None not in row and all(value is not None for value in row.values()) for row in rows),
                f"TSV row width differs from header: {path}")
        return rows

def write_table(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("NA" if row.get(key) is None else row.get(key)) for key in fields})

def safe_id(value):
    require(isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value), f"Unsafe identifier: {value}")
    return value

def child_path(base, value):
    path = Path(value)
    path = path.resolve() if path.is_absolute() else (Path(base) / path).resolve()
    _native_input(path)
    return path

def guard_output(outdir, inputs):
    out = Path(outdir).resolve()
    for source in inputs:
        source = Path(source).resolve()
        require(out != source and source not in out.parents, f"Output overlaps input/database: {out}")
    out.mkdir(parents=True, exist_ok=True)
    return out

def validate_profile(profile):
    expected = {
        "schema_version", "profile_id", "status", "aligner", "alignment_mode", "preset",
        "identity_min", "aligned_query_fraction_min", "mapq_min", "mapq_source", "base_quality_min",
        "report_k", "duplicate_policy", "secondary_policy", "supplementary_policy",
        "mate_overlap_policy", "end_exclusion", "pair_orientation", "insert_min", "insert_max",
        "ambiguity_policy", "breadth_min_depth", "mappability_read_length",
        "mappability_stride", "mappability_max_tiles_per_reference"
    }
    require(expected <= profile.keys(), f"Missing measurement fields: {sorted(expected - profile.keys())}")
    require(profile["schema_version"] == 1 and profile["status"] == "development", "Unsupported profile schema/status")
    require(profile["aligner"] == "bowtie2" and profile["alignment_mode"] == "end-to-end", "Unsupported mapper mode")
    require(profile["preset"] == "very-sensitive", "Unsupported Bowtie2 preset")
    for key in ["identity_min", "aligned_query_fraction_min"]:
        require(isinstance(profile[key], (int, float)) and 0 <= profile[key] <= 1, f"{key} must be fraction")
    for key in ["mapq_min", "base_quality_min", "insert_min", "insert_max", "report_k", "end_exclusion",
                "breadth_min_depth", "mappability_read_length", "mappability_stride", "mappability_max_tiles_per_reference"]:
        require(isinstance(profile[key], int) and profile[key] >= 0, f"Invalid integer: {key}")
    require(profile["mapq_source"] == "bowtie2_default_primary_exact_location_join", "Unsupported MAPQ source")
    require(profile["report_k"] >= 2 and profile["insert_max"] >= profile["insert_min"], "Invalid reporting/pair bounds")
    require(profile["breadth_min_depth"] >= 1 and profile["mappability_read_length"] > 0
            and profile["mappability_stride"] > 0 and profile["mappability_max_tiles_per_reference"] > 0, "Invalid coverage/tile settings")
    require(profile["pair_orientation"] in {"fr", "rf", "ff"}, "Invalid pair orientation")
    require(profile["duplicate_policy"] == "retain"
            and profile["secondary_policy"] == "audit_only"
            and profile["supplementary_policy"] == "exclude"
            and profile["mate_overlap_policy"] == "union_observed_bases"
            and profile["ambiguity_policy"] == "any_passing_alternative_location", "Unsupported measurement policies")
    return profile

def validate_clean_manifest(manifest, verify_files=True):
    require(manifest.get("schema_version") == 1 and manifest.get("preprocessing_id"), "Missing clean reads schema/preprocessing ID")
    require(isinstance(manifest.get("units"), list), "Missing clean reads units")
    seen = set()
    for unit in manifest["units"]:
        identifier = safe_id(unit.get("unit_id"))
        require(identifier not in seen, "Duplicate measurement unit")
        seen.add(identifier)
        for key in ["biological_sample_id", "library_id"]:
            require(unit.get(key), f"Missing {key}: {identifier}")
        require(unit.get("platform") == "ILLUMINA" and unit.get("molecule") == "DNA", "B1 requires Illumina DNA")
        require(unit.get("layout") in {"PE", "SE"}, "Invalid layout")
        require(unit.get("material_type") in {"bacterial_isolate", "bulk_metagenome", "vlp"}, "Invalid material")
        require(unit.get("reads_1") and (bool(unit.get("reads_2")) == (unit["layout"] == "PE")), "Reads/layout mismatch")
        counts = unit.get("counts", {})
        require(all(isinstance(counts.get(k), int) and counts[k] >= 0 for k in ["fragments", "reads", "bases"]), "Invalid clean-read counts")
        require(counts["reads"] == counts["fragments"] * (2 if unit["layout"] == "PE" else 1), "Read/fragment denominator mismatch")
        require(all(k in unit.get("eligibility", {}) for k in ["quantification", "host_adjustment", "activity"]), "Missing eligibility")
        if verify_files:
            for key in ["reads_1", "reads_2"]:
                if not unit.get(key):
                    continue
                path = Path(unit[key])
                require(path.is_absolute() and path.is_file(), f"Missing absolute clean reads: {path}")
                require(unit.get("checksums", {}).get(str(path)) == sha256(path), f"Clean reads checksum mismatch: {path}")
    return manifest
