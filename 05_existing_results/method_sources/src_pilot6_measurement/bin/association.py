#!/usr/bin/env python3
"""Validate B2 inputs, prepare strict normalized tables, and execute the pinned R model."""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import json
import math
import subprocess
from pathlib import Path

UNKNOWN = {"", "unknown", "na", "nan", "none", "null"}
STATUSES = {"assessed", "unknown", "not_assessed", "failed", "not_applicable", "ambiguous", "ineligible"}
ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path, required=()):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Missing or duplicate TSV headers: " + str(path))
        if not set(required).issubset(reader.fieldnames):
            raise ValueError("Missing columns in " + str(path) + ": " + str(set(required) - set(reader.fieldnames)))
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError("Malformed TSV row: " + str(path))
        return rows, reader.fieldnames


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_spec(path, schema_path=None):
    import jsonschema
    import yaml
    spec = yaml.safe_load(Path(path).read_text())
    schema = json.loads(Path(schema_path or ROOT / "assets/association_spec.schema.json").read_text())
    jsonschema.Draft202012Validator(schema).validate(spec)
    if spec["exposure"]["case"] == spec["exposure"]["control"]:
        raise ValueError("Exposure contrast must contain distinct case/control levels")
    if spec["input_kind"] == "genome" and not {"species_column", "target_species"}.issubset(spec):
        raise ValueError("Genome models require an explicit within-species stratum")
    if spec["input_kind"] == "cohort" and "material_type" not in spec:
        raise ValueError("Cohort models require a bulk/VLP stratum")
    if spec.get("material_type") == "vlp" and "host_adjusted" in spec["models"]:
        raise ValueError("VLP host-adjusted inference is outside this development contract")
    if spec["outcome"] in {"carriage", "count"} and spec["input_kind"] != "genome":
        raise ValueError("carriage/count require genome input")
    if spec["outcome"] in {"prevalence", "positive_abundance"} and spec["input_kind"] != "cohort":
        raise ValueError("prevalence/positive_abundance require cohort input")
    if "meta" in spec["models"] and "per_cohort" not in spec["models"]:
        raise ValueError("meta requires explicit per_cohort models")
    if "meta" in spec["models"] and not {"meta_method", "meta_ci_method"}.issubset(spec):
        raise ValueError("meta method and CI method must be explicit")
    if "host_adjusted" in spec["models"] and "host" not in spec:
        raise ValueError("host_adjusted requires a host field/transform/zero policy")
    if spec.get("host", {}).get("transform") == "log" and spec["host"]["zero_policy"] != "exclude":
        raise ValueError("log host requires exclusion of zero, no implicit pseudocount")
    if spec["repeated_subject"] == "select_first" and "selection_column" not in spec:
        raise ValueError("select_first requires a predeclared selection_column")
    if spec["repeated_subject"] == "random_intercept" and spec["independence_policy"] != "known_subject_required":
        raise ValueError("Random subject effect requires observed subject identities")
    covars = [item["name"] for item in spec["covariates"]]
    if len(covars) != len(set(covars)) or spec["exposure"]["name"] in covars:
        raise ValueError("Duplicate covariate or exposure used as its own covariate")
    if any(item["type"] == "categorical" and "reference" not in item for item in spec["covariates"]):
        raise ValueError("Categorical covariates require explicit reference levels")
    if spec["outcome"] == "count":
        if spec.get("count_offset") not in {"none", "genome_length_mb"}:
            raise ValueError("count requires explicit count_offset")
        if spec["count_offset"] == "genome_length_mb" and "length_column" not in spec:
            raise ValueError("Per-Mb count requires genome length column")
    if spec["outcome"] == "positive_abundance" and "abundance_units" not in spec:
        raise ValueError("positive abundance requires a named normalization and units")
    return spec


def prepare(long_path, metadata_path, spec_path, output, schema_path=None):
    output = Path(output)
    spec = load_spec(spec_path, schema_path)
    required = {"unit_id", spec["exposure"]["name"], spec["cohort_column"], spec["subject_column"]}
    required.update(item["name"] for item in spec["covariates"])
    required.add(spec["species_column"] if spec["input_kind"] == "genome" else "material_type")
    for key in ("selection_column", "length_column", "lineage_column"):
        if key in spec:
            required.add(spec[key])
    if "host" in spec:
        required.add(spec["host"]["column"])
    if spec["independence_policy"] == "explicit_independent_units":
        required.add("independence_confirmed")
    metadata, meta_fields = read_tsv(metadata_path, required)
    if any(row["unit_id"].lower() in UNKNOWN for row in metadata):
        raise ValueError("unit_id cannot be missing/unknown")
    ids = [row["unit_id"] for row in metadata]
    if len(ids) != len(set(ids)):
        raise ValueError("metadata unit_id must be unique; technical merging is upstream")
    collision = {"feature_id", "value", "status", "measurement_reason"} & set(meta_fields)
    if collision:
        raise ValueError("Metadata conflicts with reserved measurement columns: " + str(collision))
    if any(name.startswith("b2_") for name in meta_fields):
        raise ValueError("Metadata b2_ prefix is reserved")
    for row in metadata:
        for field in ("catalog_version", "catalog_scope", "measurement_profile", "detection_profile"):
            if field in row and row[field].lower() not in UNKNOWN and row[field] != spec[field]:
                raise ValueError("Metadata/spec version conflict for " + field)
    by_id = {row["unit_id"]: row for row in metadata}
    measurements, _ = read_tsv(long_path, {"unit_id", "feature_id", "value", "status"})
    wanted = set(spec["feature_ids"])
    pairs = set()
    joined = []
    for row in measurements:
        key = (row["unit_id"], row["feature_id"])
        if key in pairs:
            raise ValueError("Duplicate unit/feature measurement: " + str(key))
        pairs.add(key)
        if row["unit_id"] not in by_id:
            raise ValueError("Measurement has no metadata: " + row["unit_id"])
        if row["feature_id"] not in wanted:
            raise ValueError("Measurement feature absent from frozen hypothesis family: " + row["feature_id"])
        if row["status"] not in STATUSES:
            raise ValueError("Unrecognized measurement status: " + row["status"])
        if row["status"] == "assessed":
            try:
                value = float(row["value"])
            except ValueError as error:
                raise ValueError("Assessed measurement requires a number") from error
            if not math.isfinite(value) or value < 0:
                raise ValueError("Assessed measurement must be finite and nonnegative")
            if spec["outcome"] in {"carriage", "prevalence"} and value not in {0, 1}:
                raise ValueError("Binary outcome must be exactly 0 or 1")
            if spec["outcome"] == "count" and value != int(value):
                raise ValueError("Count outcome must be integer")
        elif row["value"].lower() not in UNKNOWN:
            raise ValueError("Unassessed/failed measurement cannot carry a numeric outcome")
        joined.append({**by_id[row["unit_id"]], **{k: row[k] for k in ("feature_id", "value", "status")}, "measurement_reason": row.get("reason", "")})
    # Preserve the declared denominator without inferring absence from a missing record.
    for fid in spec["feature_ids"]:
        for uid in ids:
            if (uid, fid) not in pairs:
                joined.append({**by_id[uid], "feature_id": fid, "value": "", "status": "not_assessed", "measurement_reason": "measurement_record_missing"})
    output.mkdir(parents=True, exist_ok=True)
    fields = meta_fields + ["feature_id", "value", "status", "measurement_reason"]
    write_tsv(output / "normalized_long.tsv", fields, sorted(joined, key=lambda row: (row["feature_id"], row["unit_id"])))
    (output / "analysis_spec.json").write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema_version": "1.0", "analysis_id": spec["analysis_id"],
        "scientific_status": spec["scientific_status"],
        "real_data_validated": False, "scientifically_calibrated": False,
        "input_hashes": {"long": sha(long_path), "metadata": sha(metadata_path), "spec": sha(spec_path)},
        "declared_feature_count": len(wanted), "metadata_units": len(ids), "provided_measurements": len(measurements),
        "inserted_not_assessed_records": len(joined) - len(measurements),
        "catalog_version": spec["catalog_version"], "measurement_profile": spec["measurement_profile"],
        "detection_profile": spec["detection_profile"],
        "lineage_claim": "provided_covariate" if spec.get("lineage_column") in [item["name"] for item in spec["covariates"]] else "not_adjusted",
    }
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return spec


def report(output):
    output = Path(output)
    estimates, fields = read_tsv(output / "effects.tsv")
    links = ["effects.tsv", "analysis_sample_audit.tsv", "design_report.tsv", "host_adjustment_comparison.tsv",
             "analysis_spec.json", "input_manifest.json", "sessionInfo.txt"]
    table = "<table><tr>" + "".join("<th>" + html.escape(f) + "</th>" for f in fields) + "</tr>"
    for row in estimates:
        table += "<tr>" + "".join("<td>" + html.escape(row[f]) + "</td>" for f in fields) + "</tr>"
    table += "</table>"
    text = ('<!doctype html><html><meta charset="utf-8"><title>B2 development analysis</title>'
            '<style>body{font:14px sans-serif;margin:24px}td,th{padding:4px;border:1px solid #bbb}table{border-collapse:collapse}</style>'
            '<h1>B2 development analysis</h1><p>Exploratory templates. Synthetic checks are not biological validation. '
            'Unknown/failed/unassessed outcomes are never filled as absence. See exact denominators and exclusions.</p>')
    text += "<ul>" + "".join('<li><a href="' + p + '">' + p + "</a></li>" for p in links) + "</ul>"
    (output / "association_report.html").write_text(text + table + "</html>")
    inventory = {p.name: sha(p) for p in output.iterdir() if p.is_file() and p.name != "output_manifest.json"}
    (output / "output_manifest.json").write_text(json.dumps({"files": inventory, "scientific_status": "development_exploratory"}, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--long", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--spec", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--schema", type=Path)
    parser.add_argument("--r-script", type=Path, default=ROOT / "bin/association_models.R")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    out = args.outdir.resolve()
    for source in (args.long, args.metadata, args.spec):
        if source.resolve().is_relative_to(out):
            parser.error("Output must not contain any input")
    if out.exists() and any(out.iterdir()):
        parser.error("Refusing to overwrite a nonempty analysis output")
    prepare(args.long, args.metadata, args.spec, out, args.schema)
    if not args.prepare_only:
        subprocess.run(["Rscript", str(args.r_script), str(out)], check=True)
        report(out)


if __name__ == "__main__":
    main()
