"""Aggregate B1 tables without silently filling failed/unknown units with zero."""
from __future__ import annotations
import html
import math
from pathlib import Path
from cohort_common import read_json, require, table, write_json, write_table
from cohort_measure import DETECTION_FIELDS, HOST_FIELDS, MEASUREMENT_FIELDS, apply_detection, _detection_profiles
from cohort_tools import METAPHLAN_FIELDS

def _cell(value):
    return "NA" if value is None else str(value)


def _unique(rows, keys, label):
    indexed = {}
    for row in rows:
        key = tuple(row.get(k) for k in keys)
        require(all(v not in (None, "", "NA") for v in key), f"Missing {label} key: {key}")
        require(key not in indexed, f"Duplicate {label} key: {key}")
        indexed[key] = row
    return indexed


def _number(row, key, integer=False, optional=False):
    value = row.get(key)
    if optional and value in (None, "", "NA"):
        return None
    try:
        number = float(value)
    except (ValueError, TypeError):
        require(False, f"Invalid numeric B1 field: {key}={value}")
    require(math.isfinite(number) and number >= 0, f"Invalid numeric B1 field: {key}={value}")
    require(not integer or number.is_integer(), f"Non-integer B1 count: {key}={value}")
    return int(number) if integer else number


def _validate_results(root, units, measurements, detections, hosts, manifests):
    """Reject partial/mixed result sets before emitting any successful summary."""
    reference = read_json(root / "reference/reference_manifest.json")
    profile = read_json(root / "measurement_profile.json")
    detection_config = read_json(root / "detection_profiles.json")
    profiles = _detection_profiles(detection_config)
    target_key = 'reference_id' if reference.get('measurement_contract')=='research_quantification_group_v1' else 'votu_id'
    targets = _unique([r for r in reference["references"] if r["reference_kind"] == "target"],
                      [target_key], "reference target")
    reference_hosts = {r["host_id"] for r in reference["references"] if r["reference_kind"] == "host"}
    for manifest in manifests:
        unit = units[manifest["unit"]["unit_id"]]
        require(manifest["unit"] == unit, "B1 unit manifest differs from prepared input")
        require(manifest.get("reference_version") == reference["mapping_reference_version"], "B1 reference version mismatch")
        require(manifest.get("profile") == profile, "B1 measurement profile mismatch")
    measured = _unique(measurements, ["unit_id", target_key], "measurement")
    detected = _unique(detections, ["unit_id", target_key, "detection_profile_id"], "detection")
    host_index = _unique(hosts, ["unit_id", "host_id"], "host")
    # Cardinality + per-key membership avoids materializing another unit x target set.
    require(len(measured) == len(units) * len(targets), "Missing/extra B1 measurement rows")
    require(len(detected) == len(measured) * len(profiles), "Missing/extra B1 detection rows")
    require(len(host_index) == len(units) * len(reference_hosts), "Missing/extra B1 host rows")
    for (unit_id, votu), row in measured.items():
        require(unit_id in units and (votu,) in targets, "Unknown B1 unit/vOTU row")
        unit, target = units[unit_id], targets[(votu,)]
        identity = {k: unit.get(k) for k in ["unit_id", "biological_sample_id", "library_id", "subject_id", "material_type", "layout"]}
        identity.update({k: reference[k] for k in ["catalog_version", "catalog_scope", "mapping_reference_version"]})
        identity.update(profile_id=profile["profile_id"], reference_id=target["reference_id"], candidate_id=target["candidate_id"])
        require(all(row.get(k) == _cell(v) for k, v in identity.items()), "B1 row identity/provenance mismatch")
        numeric = dict(row)
        for key in ["target_length", "assessable_length", "mapped_read_count", "assigned_fragment_count",
                    "ambiguous_fragment_count", "covered_bases_full", "covered_bases_assessable", "clean_fragments_denominator"]:
            numeric[key] = _number(row, key, integer=True)
        for key in ["mean_depth_full", "breadth_full"]:
            numeric[key] = _number(row, key)
        for key in ["mean_depth_assessable", "breadth_assessable", "fragment_rpkm_all_clean", "unique_mappability"]:
            numeric[key] = _number(row, key, optional=True)
        if row.get('informative_status'):
            for key in ('informative_length','covered_bases_informative'):numeric[key]=_number(row,key,integer=True)
            for key in ('breadth_informative','mean_depth_informative'):numeric[key]=_number(row,key,optional=True)
            require(numeric['covered_bases_informative']<=numeric['informative_length']<=numeric['target_length'],'Informative coverage exceeds length')
            if numeric['informative_length']:
                require(numeric['breadth_informative'] is not None and math.isclose(numeric['breadth_informative'],numeric['covered_bases_informative']/numeric['informative_length'],rel_tol=1e-9),'Informative breadth/count mismatch')
        require(numeric["target_length"] == target["length"] and numeric["target_length"] > 0, "B1 target length mismatch")
        require(numeric["clean_fragments_denominator"] == unit["counts"]["fragments"], "B1 clean denominator mismatch")
        require(numeric["assigned_fragment_count"] <= numeric["clean_fragments_denominator"], "B1 assigned count exceeds denominator")
        require(numeric["mapped_read_count"] == numeric["assigned_fragment_count"] * (2 if unit["layout"] == "PE" else 1), "B1 read/fragment count mismatch")
        require(numeric["covered_bases_assessable"] <= numeric["assessable_length"] <= numeric["target_length"]
                and numeric["covered_bases_full"] <= numeric["target_length"], "B1 coverage exceeds length")
        for key in ["breadth_full", "breadth_assessable", "unique_mappability"]:
            require(numeric[key] is None or numeric[key] <= 1, f"B1 fraction outside 0..1: {key}")
        require(math.isclose(numeric["breadth_full"], numeric["covered_bases_full"] / numeric["target_length"],
                             rel_tol=1e-9, abs_tol=1e-12), "B1 breadth/count mismatch")
        denominator = numeric["clean_fragments_denominator"]
        amount = numeric["fragment_rpkm_all_clean"]
        require((amount is None and denominator == 0) or
                (amount is not None and denominator > 0 and math.isclose(
                    amount, numeric["assigned_fragment_count"] * 1e9 / (numeric["target_length"] * denominator),
                    rel_tol=1e-9, abs_tol=1e-12)), "B1 RPKM/denominator mismatch")
        expected_status = "assessed" if unit["counts"]["fragments"] and unit["eligibility"]["quantification"]["status"] in {"eligible", "provisional"} else "not_assessed"
        require(row["measurement_status"] == expected_status, "B1 assessment status mismatch")
        for expected in apply_detection(numeric, detection_config):
            key = (unit_id, votu, expected["detection_profile_id"])
            require(key in detected, "Missing B1 target/profile detection")
            actual = detected[key]
            require(all(actual.get(k) == _cell(expected.get(k)) for k in DETECTION_FIELDS), "B1 detection differs from measurement/profile")
    for (unit_id, host_id), row in host_index.items():
        require(unit_id in units and host_id in reference_hosts, "Unknown B1 host row")
        require(all(row.get(k) == _cell(units[unit_id].get(k)) for k in
                    ["biological_sample_id", "library_id", "material_type"]), "B1 host identity mismatch")


def aggregate(prepared, measurement_dirs, metaphlan_dirs, outdir, metaphlan_disabled_reason="not_enabled"):
    root, out = Path(prepared), Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    clean = read_json(root / "clean_reads_manifest.json")
    units = {u["unit_id"]: u for u in clean["units"]}
    require(len(units) == len(clean["units"]), "Duplicate prepared B1 unit")
    measurements, detections, hosts, manifests = [], [], [], []
    seen = set()
    for folder in measurement_dirs:
        folder = Path(folder)
        manifest = read_json(folder / "measurement_manifest.json")
        identifier = manifest["unit"]["unit_id"]
        require(identifier in units and identifier not in seen, "Unknown/duplicate B1 result unit")
        seen.add(identifier)
        manifests.append(manifest)
        for filename, destination in [("measurement_long.tsv", measurements), ("detection_long.tsv", detections),
                                      ("host_backbone.tsv", hosts)]:
            rows = table(folder / filename)
            require(all(row.get("unit_id") == identifier for row in rows), "B1 row does not belong to result unit")
            destination.extend(rows)
    require(seen == set(units), "Missing B1 unit results; refusing incomplete success")
    _validate_results(root, units, measurements, detections, hosts, manifests)
    if read_json(root/'reference/reference_manifest.json').get('measurement_contract')=='research_quantification_group_v1':
        require(not metaphlan_dirs,'Research host-core route does not merge unrelated MetaPhlAn abundance')
        from cohort_research_report import aggregate_groups
        return aggregate_groups(root,out,units,measurements,detections,hosts,manifests)
    species, mp_status = [], {}
    for folder in metaphlan_dirs:
        folder = Path(folder)
        receipt = read_json(folder / "metaphlan_manifest.json")
        identifier = receipt["unit_id"]
        require(identifier in units and identifier not in mp_status, "Unknown/duplicate MetaPhlAn result")
        rows = table(folder / "metaphlan_species.tsv")
        require(all(row.get("unit_id") == identifier for row in rows), "MetaPhlAn row does not belong to result unit")
        _unique(rows, ["unit_id", "clade_name"], "MetaPhlAn")
        species.extend(rows)
        mp_status[identifier] = receipt["status"]
    measurements.sort(key=lambda row: (row["unit_id"], row["votu_id"]))
    detections.sort(key=lambda row: (row["unit_id"], row["votu_id"], row["detection_profile_id"]))
    hosts.sort(key=lambda row: (row["unit_id"], row["host_id"]))
    species.sort(key=lambda row: (row["unit_id"], row["clade_name"]))
    manifests.sort(key=lambda row: row["unit"]["unit_id"])
    write_table(out / "measurement_long.tsv", MEASUREMENT_FIELDS, measurements)
    write_table(out / "detection_long.tsv", DETECTION_FIELDS, detections)
    write_table(out / "host_backbone.tsv", HOST_FIELDS, hosts)
    write_table(out / "metaphlan_species.tsv", METAPHLAN_FIELDS, species)
    status = []
    metadata = []
    for identifier, unit in sorted(units.items()):
        metadata_row = dict(unit.get("metadata", {}))
        metadata_row.update({key: unit.get(key) for key in ["unit_id", "biological_sample_id", "library_id",
                             "subject_id", "material_type", "layout", "platform", "molecule"]})
        metadata_row.update({"preprocessing_id": clean["preprocessing_id"],
                             "quantification_eligibility": unit["eligibility"]["quantification"]["status"],
                             "host_adjustment_eligibility": unit["eligibility"]["host_adjustment"]["status"]})
        metadata.append(metadata_row)
        status.append({"unit_id": identifier, "measurement": "completed", "metaphlan": mp_status.get(identifier,
                       "not_applicable" if unit["material_type"] != "bulk_metagenome" else metaphlan_disabled_reason),
                       "real_data_validated": False, "scientifically_calibrated": False})
    fields = sorted(set().union(*(row.keys() for row in metadata))) if metadata else ["unit_id"]
    write_table(out / "unit_metadata.tsv", fields, metadata)
    write_table(out / "unit_status.tsv", ["unit_id", "measurement", "metaphlan", "real_data_validated", "scientifically_calibrated"], status)
    by_key = {(r["unit_id"], r["votu_id"]): r for r in measurements}
    profiles = read_json(root / "detection_profiles.json")["profiles"]
    for profile in profiles:
        profile_id = profile["profile_id"]
        prevalence, positive = [], []
        for row in detections:
            if row["detection_profile_id"] != profile_id:
                continue
            assessed = row["inference_status"] == "assessed"
            prevalence.append({"unit_id": row["unit_id"], "feature_id": row["votu_id"],
                               "value": row["value"] if assessed else "", "status": "assessed" if assessed else "unknown", "reason": row["inference_reason"]})
            amount = by_key[(row["unit_id"], row["votu_id"])]["fragment_rpkm_all_clean"]
            positive.append({"unit_id": row["unit_id"], "feature_id": row["votu_id"],
                             "value": amount if assessed and row["detection_status"] == "detected" else "",
                             "status": "assessed" if assessed and row["detection_status"] == "detected" else "unknown", "reason": row["inference_reason"]})
        write_table(out / f"b2_prevalence_{profile_id}.tsv", ["unit_id", "feature_id", "value", "status", "reason"], prevalence)
        write_table(out / f"b2_positive_abundance_{profile_id}.tsv", ["unit_id", "feature_id", "value", "status", "reason"], positive)
        features = sorted({r["feature_id"] for r in prevalence})
        unit_values = {identifier: {} for identifier in units}
        for row in prevalence:
            unit_values[row["unit_id"]][row["feature_id"]] = row["value"]
        matrix = [{"unit_id": identifier, **unit_values[identifier]} for identifier in sorted(units)]
        write_table(out / f"detection_matrix_{profile_id}.tsv", ["unit_id", *features], matrix)
    for metric in ["mean_depth_full", "breadth_full", "fragment_rpkm_all_clean"]:
        features = sorted({r["votu_id"] for r in measurements})
        matrix = []
        for identifier in sorted(units):
            row = {"unit_id": identifier}
            for feature in features:
                measurement = by_key.get((identifier, feature))
                row[feature] = measurement[metric] if measurement and measurement["measurement_status"] == "assessed" else None
            matrix.append(row)
        write_table(out / f"{metric}_matrix.tsv", ["unit_id", *features], matrix)
    write_json(out / "cohort_manifest.json", {"schema_version": 1, "stage": "B1", "build_state": "development",
               "prepared_manifest": read_json(root / "prepared_manifest.json"), "units": len(units),
               "measurement_rows": len(measurements), "detection_rows": len(detections),
               "unit_manifests": manifests, "metaphlan_status": status,
               "aggregation_integrity": {"status": "passed", "contract": "b1_complete_identity_numeric_detection_v1",
                   "checks": ["unit_ownership", "unique_keys", "complete_target_profile_sets", "reference_and_profile_identity",
                              "numeric_counts_and_denominators", "detection_recomputed_from_measurements"]},
               "b2_positive_abundance_unit": "fragment_rpkm_all_clean_on_detected_units_only",
               "real_data_validated": False, "scientifically_calibrated": False})
    (out / "OUTPUTS.md").write_text(
        "# B1 output guide\n\n"
        "Start with report_cohort.html, unit_status.tsv and cohort_manifest.json. "
        "A unit is a readset or explicitly merged library, not necessarily an independent subject.\n\n"
        "| File | Grain / use |\n|---|---|\n"
        "| measurement_long.tsv | One unit_id x votu_id; continuous measurements before detection filtering |\n"
        "| detection_long.tsv | One unit_id x votu_id x detection_profile_id; provisional detection and inference eligibility |\n"
        "| *_matrix.tsv | Rows are unit_id; feature columns are votu_id, sorted by ID |\n"
        "| b2_prevalence_*.tsv | unit_id x feature_id; assessed 0/1 or unknown with blank value |\n"
        "| b2_positive_abundance_*.tsv | RPKM only for detected, inference-assessable units; other values remain blank |\n"
        "| host_backbone.tsv | unit_id x host_id; masked host depth and breadth, not predicted viral-host pairs |\n"
        "| metaphlan_species.tsv | unit_id x clade_name; species/unclassified relative abundance, percent and fraction |\n"
        "| unit_metadata.tsv / unit_status.tsv | Unit identity, eligibility and execution state |\n\n"
        "## Measurement definitions\n\n"
        "assigned_fragment_count counts each retained pair once (or one SE read). mapped_read_count counts its mates. "
        "Ambiguous fragments do not enter assigned counts or coverage; their support is retained separately. "
        "mean_depth_full is observed-fragment base depth: overlapping mates count once; deletions and unsequenced gaps do not count. "
        "breadth_full = covered_bases_full / target_length, a 0..1 fraction, not sequencing depth or relative abundance. "
        "Assessable metrics instead use assessable_length after masks/end rules. "
        "fragment_rpkm_all_clean = assigned_fragment_count * 1e9 / (target_length * clean_fragments_denominator); "
        "unassigned/ambiguous fragments remain in the clean denominator. This is not absolute abundance, copy number, or TMM.\n\n"
        "## IDs, missingness and limits\n\n"
        "Join by IDs, never row order. reference_id, candidate_id and votu_id are distinct; follow the frozen reference map. "
        "Multiple source occurrences or host candidates do not justify duplicating representative abundance. "
        "NA/blank is unavailable; unknown/not_assessed/failed is not biological zero. An assessed zero means not detected by this assay. "
        "Shared-only signal remains not_resolvable; B2 receives unknown. A zero clean denominator remains unassessable. "
        "Missing units/targets/profiles, duplicate keys or inconsistent provenance now fail aggregation.\n\n"
        "cohort_manifest.json retains prepared and unit manifests plus aggregation_integrity checks. "
        "Reference/catalog versions and measurement/detection profile IDs are recorded in the long tables. "
        "Engineering completion does not establish scientific calibration. Detection settings remain provisional; "
        "VLP and bulk have different biases. Detection does not establish integration, a source locus, lifestyle, AMG activity or causality.\n",
        encoding="utf-8")
    rows = "".join("<tr>" + "".join(f"<td>{html.escape(str(row.get(k, '')))}</td>" for k in
                   ["unit_id", "measurement", "metaphlan"]) + "</tr>" for row in status)
    (out / "report_cohort.html").write_text(
        "<!doctype html><html lang='en'><meta charset='utf-8'><title>B1 cohort measurements</title>"
        "<style>body{font:16px system-ui;max-width:1050px;margin:40px auto}td,th{padding:8px;border:1px solid #ddd}table{border-collapse:collapse}</style>"
        "<h1>B1: development measurements</h1><p>Units are readsets or explicitly merged libraries, not automatically independent patients. "
        "Detection thresholds are provisional. Shared-only signal is not resolvable and exports as unknown, not absence. Relative signals are not absolute viral counts, integration calls or cargo linkage.</p>"
        f"<p>{len(units)} measurement units; {len(measurements)} unit/target rows.</p>"
        "<table><tr><th>Unit</th><th>Measurement</th><th>MetaPhlAn</th></tr>" + rows + "</table>"
        "<p><a href='OUTPUTS.md'>Output definitions and usage</a> · <a href='measurement_long.tsv'>Continuous measurements</a> · <a href='detection_long.tsv'>Detection</a> · "
        "<a href='host_backbone.tsv'>Host backbone</a> · <a href='metaphlan_species.tsv'>Species relative abundance</a> · "
        "<a href='cohort_manifest.json'>Manifest</a></p></html>")
    return status
