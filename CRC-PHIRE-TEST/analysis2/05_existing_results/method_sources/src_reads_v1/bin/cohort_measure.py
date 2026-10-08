"""Queryname-grouped SAM/BAM counting with explicit fragment and CIGAR semantics."""
from __future__ import annotations
from array import array
from collections import Counter
from contextlib import contextmanager, ExitStack
from itertools import groupby, zip_longest
from pathlib import Path
import re
import subprocess
from cohort_common import ContractError, digest, fasta, read_json, require, validate_profile, write_json, write_table, sha256

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
BASES = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
MEASUREMENT_FIELDS = [
    "unit_id", "biological_sample_id", "library_id", "subject_id", "material_type", "layout",
    "catalog_version", "catalog_scope", "mapping_reference_version", "profile_id",
    "reference_id", "votu_id", "candidate_id", "virus_scope", "target_length", "assessable_length",
    "mapped_read_count", "assigned_fragment_count", "ambiguous_fragment_count", "covered_bases_full", "covered_bases_assessable",
    "mean_depth_full", "mean_depth_assessable", "breadth_full", "breadth_assessable",
    "fragment_rpkm_all_clean", "clean_fragments_denominator", "measurement_status", "reason",
    "measurement_target_level", "resolution_status", "unique_mappability", "mappability_status",
    "competition_evidence", "sample_in_discovery_set"
]
HOST_FIELDS = ["unit_id", "biological_sample_id", "library_id", "material_type", "host_id",
               "assessable_length", "covered_bases", "mean_backbone_depth", "backbone_breadth",
               "assigned_fragment_count", "host_adjustment_eligibility", "status", "reason"]
DETECTION_FIELDS = MEASUREMENT_FIELDS[:11] + [
    "votu_id", "detection_profile_id", "detection_profile_status", "detection_status", "value", "reason",
    "resolution_status", "ambiguous_fragment_count", "inference_status", "inference_reason"
]

RESEARCH_FIELDS = ['quantification_group','biological_votu_ids','reference_host_id','host_assignment_status']
MEASUREMENT_FIELDS += RESEARCH_FIELDS
MEASUREMENT_FIELDS += ['informative_length','covered_bases_informative','breadth_informative','mean_depth_informative','informative_status']
DETECTION_FIELDS += RESEARCH_FIELDS

@contextmanager
def sam_stream(path):
    path = Path(path)
    if path.suffix.lower() in {".bam", ".cram"}:
        require(path.suffix.lower() != ".cram", "CRAM needs explicit reference support; use BAM/SAM")
        process = subprocess.Popen(["samtools", "view", "-h", str(path)], stdout=subprocess.PIPE, text=True)
        try:
            yield process.stdout
        finally:
            process.stdout.close()
            result = process.wait()
            require(result == 0, "samtools view failed")
    else:
        with path.open() as handle:
            yield handle

def parse_sam(line):
    cells = line.rstrip("\n").split("\t")
    require(len(cells) >= 11, "Malformed SAM record")
    record = {"qname": cells[0], "flag": int(cells[1]), "rname": cells[2], "pos": int(cells[3])-1,
              "mapq": int(cells[4]), "cigar": cells[5], "rnext": cells[6], "pnext": int(cells[7])-1,
              "tlen": int(cells[8]), "seq": cells[9].upper(), "qual": cells[10]}
    record["mate"] = 1 if record["flag"] & 64 else (2 if record["flag"] & 128 else 0)
    record["tags"] = {x[:2]: x[5:] for x in cells[11:] if len(x) > 5}
    return record

def sam_groups(handle):
    def records():
        for line in handle:
            if not line.strip() or line.startswith("@"):
                continue
            yield parse_sam(line)
    for name, rows in groupby(records(), key=lambda r: r["qname"]):
        yield name, list(rows)

def fill_secondary_sequence(record, primary_by_mate):
    if record["seq"] != "*":
        return record
    require(record["mate"] in primary_by_mate, "Secondary alignment lacks reconstructable read")
    primary = primary_by_mate[record["mate"]]
    require(primary["seq"] != "*", "Primary alignment lacks sequence")
    record = dict(record)
    record["seq"], record["qual"] = primary["seq"], primary["qual"]
    if bool(record["flag"] & 16) != bool(primary["flag"] & 16):
        record["seq"] = record["seq"].translate(BASES)[::-1]
        record["qual"] = record["qual"][::-1]
    return record

def alignment_metrics(record, references, profile):
    if record["flag"] & 4:
        return None
    require(record["rname"] in references, f"SAM has unknown reference: {record['rname']}")
    operations = CIGAR.findall(record["cigar"])
    require(operations and "".join(n+op for n, op in operations) == record["cigar"], "Malformed CIGAR")
    sequence = record["seq"]
    require(sequence != "*" and record["qual"] != "*" and len(record["qual"]) == len(sequence), "Sequence/quality missing")
    reference = references[record["rname"]]
    qpos, rpos = 0, record["pos"]
    require(rpos >= 0, "Negative alignment position")
    matched, columns, aligned_query, hard_clipped = 0, 0, 0, 0
    covered = set()
    for count_text, op in operations:
        count = int(count_text)
        require(count > 0, "Zero-length CIGAR operation")
        if op in "M=X":
            require(qpos + count <= len(sequence) and rpos + count <= len(reference), "CIGAR outside sequence/reference")
            for offset in range(count):
                qbase, rbase = sequence[qpos+offset], reference[rpos+offset]
                matched += int(qbase in "ACGT" and qbase == rbase)
                if (ord(record["qual"][qpos+offset])-33 >= profile["base_quality_min"]
                        and qbase in "ACGT" and rbase in "ACGT"):
                    covered.add(rpos+offset)
            columns += count
            aligned_query += count
            qpos += count
            rpos += count
        elif op == "I":
            qpos += count
            aligned_query += count
            columns += count
        elif op == "D":
            rpos += count
            columns += count
        elif op == "N":
            rpos += count  # skipped reference never contributes observed coverage
        elif op == "S":
            qpos += count
        elif op == "H":
            hard_clipped += count
        elif op == "P":
            pass
        require(qpos <= len(sequence) and rpos <= len(reference), "CIGAR outside bounds")
    require(qpos == len(sequence) and columns > 0, "CIGAR/read length mismatch")
    return {"identity": matched/columns, "aligned_query_fraction": aligned_query/(len(sequence)+hard_clipped),
            "covered": covered, "location": (record["rname"], record["pos"], record["cigar"], bool(record["flag"] & 16)),
            "aligned_query": aligned_query}

def passes_sequence(metrics, profile):
    return metrics and metrics["identity"] >= profile["identity_min"] and metrics["aligned_query_fraction"] >= profile["aligned_query_fraction_min"]

def pair_valid(first, second, profile):
    if not (first["flag"] & 2 and second["flag"] & 2):
        return False
    if first["rname"] != second["rname"]:
        return False
    if first["rnext"] not in {"=", second["rname"]} or second["rnext"] not in {"=", first["rname"]}:
        return False
    if first["pnext"] != second["pos"] or second["pnext"] != first["pos"]:
        return False
    if first["tlen"] != -second["tlen"] or not profile["insert_min"] <= abs(first["tlen"]) <= profile["insert_max"]:
        return False
    left, right = sorted([first, second], key=lambda r: (r["pos"], r["mate"]))
    directions = (bool(left["flag"] & 16), bool(right["flag"] & 16))
    return directions == {"fr": (False, True), "rf": (True, False), "ff": (False, False)}[profile["pair_orientation"]]

def classify_fragment(records, references, profile, layout, primary_records=None):
    primary = [r for r in records if not r["flag"] & (256 | 2048)]
    expected = {1, 2} if layout == "PE" else {0}
    require(len(primary) == len(expected) and {r["mate"] for r in primary} == expected,
            "Duplicate QNAME, missing mate, or SAM layout mismatch")
    primary_by_mate = {r["mate"]: r for r in primary}
    ordinary = None
    if primary_records is not None:
        ordinary_rows = [r for r in primary_records if not r["flag"] & (256 | 2048)]
        require(len(ordinary_rows) == len(expected)
                and {r["mate"] for r in ordinary_rows} == expected,
                "Duplicate default-primary record, missing mate, or primary/audit layout mismatch")
        require({r["qname"] for r in ordinary_rows} == {r["qname"] for r in primary}
                and len({r["qname"] for r in primary}) == 1,
                "Primary/audit QNAME identity mismatch")
        ordinary = {r["mate"]: r for r in ordinary_rows}
    if any(r["flag"] & 4 for r in primary):
        return "unmapped_or_incomplete_pair", None, set(), []
    if any(r["flag"] & 512 for r in primary):
        return "sequencing_qc_fail", None, set(), []
    primary_metrics = {r["mate"]: alignment_metrics(r, references, profile) for r in primary}
    if not all(passes_sequence(m, profile) for m in primary_metrics.values()):
        return "sequence_filter", None, set(), []
    passing = []
    reported = Counter()
    for raw in records:
        if raw["flag"] & (4 | 2048):
            continue
        reported[raw["mate"]] += 1
        row = fill_secondary_sequence(raw, primary_by_mate)
        metrics = alignment_metrics(row, references, profile)
        if passes_sequence(metrics, profile):
            passing.append((row, metrics))
    alternatives = [(r, m) for r, m in passing if m["location"] != primary_metrics[r["mate"]]["location"]]
    if alternatives or any(n >= profile["report_k"] for n in reported.values()):
        ids = sorted({r["rname"] for r, _ in passing})
        return "ambiguous_or_reporting_saturated", None, set(), ids
    if ordinary is not None:
        joined = []
        for row in primary:
            other = ordinary[row["mate"]]
            if any(row[key] != other[key] for key in ["rname", "pos", "cigar"]) or bool(row["flag"] & 16) != bool(other["flag"] & 16):
                return "primary_alignment_not_reconciled", None, set(), []
            joined.append({**row, "mapq": other["mapq"]})
        primary = joined
    if any(r["mapq"] == 255 for r in primary):
        return "mapq_unavailable", None, set(), []
    if any(r["mapq"] < profile["mapq_min"] for r in primary):
        return "mapq_filter", None, set(), []
    if layout == "PE" and not pair_valid(primary_by_mate[1], primary_by_mate[2], profile):
        return "discordant_pair", None, set(), []
    reference = primary[0]["rname"]
    positions = set().union(*(m["covered"] for m in primary_metrics.values()))
    end = profile["end_exclusion"]
    positions = {p for p in positions if end <= p < len(references[reference])-end}
    return "assigned", reference, positions, []

def _detection_profiles(data):
    require(data.get("schema_version") == 1 and isinstance(data.get("profiles"), list), "Invalid detection profiles")
    seen = set()
    for row in data["profiles"]:
        required = {"profile_id", "status", "breadth_field", "breadth_min", "assigned_fragments_min", "mean_depth_min",
                    "require_mappability", "unique_mappability_min"}
        require(required <= row.keys(), "Incomplete detection profile")
        require(row["profile_id"] not in seen, "Duplicate detection profile")
        seen.add(row["profile_id"])
        require(row["status"] == "provisional", "Uncalibrated implementation only supports provisional detection profiles")
        require(row["breadth_field"] in {"breadth_full", "breadth_assessable", "breadth_informative"}, "Unsupported breadth denominator")
        if row['breadth_field']=='breadth_informative':
            require(isinstance(row.get('informative_length_min'),int) and row['informative_length_min']>0,'Explicit minimum informative length required')
        require(0 <= row["breadth_min"] <= 1 and row["assigned_fragments_min"] >= 0 and row["mean_depth_min"] >= 0,
                "Invalid detection threshold")
        require(0 <= row["unique_mappability_min"] <= 1, "Invalid mappability threshold")
    return data["profiles"]

def apply_detection(row, profiles):
    result = []
    for profile in _detection_profiles(profiles):
        out = {key: row.get(key) for key in DETECTION_FIELDS}
        out.update(detection_profile_id=profile["profile_id"], detection_profile_status="provisional")
        if row["measurement_status"] != "assessed" or not row["assessable_length"]:
            status, reason, value = "unassessable", row["reason"] or "zero_assessable_length", None
        elif profile["require_mappability"] and (row["mappability_status"] != "applicable"
              or row["unique_mappability"] is None or row["unique_mappability"] < profile["unique_mappability_min"]):
            status, reason, value = "unassessable", "mappability_requirement", None
        elif profile['breadth_field']=='breadth_informative' and (row.get('informative_status')!='applicable' or row.get('informative_length',0)<profile['informative_length_min']):
            status,reason,value='unassessable','informative_domain_or_length_requirement',None
        else:
            detected = (row[profile["breadth_field"]] >= profile["breadth_min"]
                        and row["assigned_fragment_count"] >= profile["assigned_fragments_min"]
                        and row['mean_depth_informative' if profile['breadth_field']=='breadth_informative' else 'mean_depth_full'] >= profile["mean_depth_min"])
            status, reason, value = ("detected", "provisional_threshold_met", 1) if detected else ("not_detected", "provisional_threshold_not_met", 0)
        # A resolved-signal zero is not evidence of absence when shared support
        # cannot be assigned to this reference. Keep raw measurements unchanged.
        if status == "not_detected" and row.get("ambiguous_fragment_count", 0):
            if row["assigned_fragment_count"] == 0:
                status, value, reason = "not_resolvable", None, "ambiguous_only_not_detected_at_resolved_signal"
            else:
                status, reason = "not_detected_at_resolved_signal", "shared_signal_excluded_below_resolved_detection_threshold"
        inference = "assessed" if status in {"detected", "not_detected"} else "unknown"
        inference_reason = "provisional_detection_profile" if inference == "assessed" else reason
        # This is an inference eligibility flag, not a new detection threshold.
        # Inapplicable PE tile estimates are reported without pretending calibration.
        if (row.get("mappability_status") == "applicable" and row.get("unique_mappability") is not None
                and row["unique_mappability"] < profile["unique_mappability_min"]):
            inference, inference_reason = "unknown", "insufficient_unique_mappability_for_b2"
        out.update(detection_status=status, value=value, reason=reason,
                   inference_status=inference, inference_reason=inference_reason)
        result.append(out)
    return result

def measure(sam_path, reference_dir, unit_file, profile_file, detections_file, outdir, mappability_file=None, primary_sam_path=None):
    manifest = read_json(Path(reference_dir) / "reference_manifest.json")
    references = fasta(Path(reference_dir) / "reference.fna")
    units = read_json(unit_file)
    profile = validate_profile(read_json(profile_file))
    detections = read_json(detections_file)
    _detection_profiles(detections)
    mappability = {};evidence={}
    if mappability_file:
        evidence = read_json(mappability_file)
        if evidence.get('schema_version')==2:
            require(evidence.get('reference_fasta_sha256')==sha256(Path(reference_dir)/'reference.fna'),'Informative reference content differs')
        require(evidence.get("mapping_reference_version") == manifest["mapping_reference_version"],
                "Mappability reference version mismatch")
        require(evidence.get("measurement_profile") == profile, "Mappability measurement profile mismatch")
        mappability = {r["reference_id"]: r for r in evidence["references"]}
        targets = {r["reference_id"] for r in manifest["references"] if r["reference_kind"] == "target"}
        require(len(mappability) == len(evidence["references"]) and set(mappability) == targets,
                "Mappability target set incomplete/duplicate/unknown")
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    lengths = {r["reference_id"]: r["length"] for r in manifest["references"]}
    require(set(lengths) == set(references), "Reference manifest/FASTA mismatch")
    coverage = {key: array("I", [0]) * length for key, length in lengths.items()}
    counts = Counter()
    assigned = Counter()
    ambiguous = Counter()
    groups = 0
    read_lengths = Counter()
    pair_spans=Counter()
    with (out / "ambiguity.tsv").open("w") as ambiguity_handle:
        ambiguity_handle.write("unit_id\tfragment_id\tambiguity_group_id\treference_ids\n")
        with ExitStack() as stack:
            handle = stack.enter_context(sam_stream(sam_path))
            if primary_sam_path:
                primary_handle = stack.enter_context(sam_stream(primary_sam_path))
                groups_iter = zip_longest(sam_groups(handle), sam_groups(primary_handle))
            else:
                groups_iter = ((item, None) for item in sam_groups(handle))
            for audit_group, primary_group in groups_iter:
                require(audit_group is not None and (not primary_sam_path or primary_group is not None), "Primary/audit group count differs")
                name, records = audit_group
                require(primary_group is None or primary_group[0] == name, "Primary/audit QNAME sort mismatch")
                ordinary = primary_group[1] if primary_group else None
                groups += 1
                for record in records:
                    if not record["flag"] & (256 | 2048) and record["seq"] != "*":
                        read_lengths[len(record["seq"])] += 1
                        if units['layout']=='PE' and record['mate']==1 and record['flag'] & 2 and not record['flag'] & 4:pair_spans[abs(record['tlen'])]+=1
                reason, target, positions, alternatives = classify_fragment(records, references, profile, units["layout"], ordinary)
                counts[reason] += 1
                if reason == "assigned":
                    assigned[target] += 1
                    for position in positions:
                        coverage[target][position] += 1
                if alternatives:
                    group_id = "amb_" + digest(alternatives)[:20]
                    ambiguity_handle.write(f"{units['unit_id']}\t{name}\t{group_id}\t{','.join(alternatives)}\n")
                    for target in alternatives:
                        ambiguous[target] += 1
    require(groups == units["counts"]["fragments"], "SAM queryname groups differ from clean fragment denominator")
    identity = {k: units.get(k) for k in ["unit_id", "biological_sample_id", "library_id", "subject_id", "material_type", "layout"]}
    total = units["counts"]["fragments"]
    qstatus = units["eligibility"]["quantification"]["status"]
    assessed = total > 0 and qstatus in {"eligible", "provisional"}
    rows, host_rows, detection_rows = [], [], []
    for reference in manifest["references"]:
        identifier = reference["reference_id"]
        cov = coverage[identifier]
        segments = [(max(start, profile["end_exclusion"]), min(end, len(cov)-profile["end_exclusion"]))
                    for start, end in reference["assessable_intervals"]]
        assessable_length = sum(max(0, end-start) for start, end in segments)
        covered_full = sum(depth >= profile["breadth_min_depth"] for depth in cov)
        covered_assessable = sum(cov[pos] >= profile["breadth_min_depth"] for start, end in segments for pos in range(start, end))
        sum_assessable = sum(cov[pos] for start, end in segments for pos in range(start, end))
        if reference["reference_kind"] == "host":
            host_rows.append({**identity, "host_id": reference["host_id"], "assessable_length": assessable_length,
                              "covered_bases": covered_assessable, "depth_sum": sum_assessable,
                              "assigned_fragment_count": assigned[identifier]})
            continue
        if reference["reference_kind"] != "target":
            continue
        maprow = mappability.get(identifier, {})
        mapfraction = maprow.get("unique_tile_fraction")
        applicable = set(read_lengths) == {profile["mappability_read_length"]} and units["layout"] == "SE"
        mapstatus = "applicable" if maprow.get("status") == "estimated" and applicable else ("reference_estimate_only" if maprow else "not_assessed")
        info=[];infostatus='not_assessed'
        if evidence.get('schema_version')==2:
            domain=maprow.get('layouts',{}).get(units['layout'],{})
            applicable=(set(read_lengths)=={domain.get('read_length')} and (units['layout']=='SE' or set(pair_spans)=={domain.get('fragment_length')}))
            mapfraction=domain.get('unique_tile_fraction')
            mapstatus='applicable' if domain.get('status')=='complete' and applicable else 'reference_estimate_only'
            infostatus=mapstatus
            if applicable and domain.get('status')=='complete':
                for start,end in domain['informative_intervals']:
                    require(0<=start<end<=len(cov),'Informative coordinate outside reference')
                    info.extend(pos for pos in range(start,end) if profile['end_exclusion']<=pos<len(cov)-profile['end_exclusion'])
                require(len(info)==len(set(info)),'Overlapping informative intervals')
        info_covered=sum(cov[pos]>=profile['breadth_min_depth'] for pos in info)
        discovery = manifest.get("discovery_samples", [])
        discovery_ids = [r.get("biological_sample_id") if isinstance(r, dict) else r for r in discovery]
        row = {**identity, "catalog_version": manifest["catalog_version"], "catalog_scope": manifest["catalog_scope"],
               "mapping_reference_version": manifest["mapping_reference_version"], "profile_id": profile["profile_id"],
               "reference_id": identifier, "votu_id": reference["votu_id"], "candidate_id": reference["candidate_id"],
               "virus_scope": reference.get("virus_scope", "unknown"), "target_length": len(cov), "assessable_length": assessable_length,
               "mapped_read_count": assigned[identifier] * (2 if units["layout"] == "PE" else 1),
               "assigned_fragment_count": assigned[identifier], "ambiguous_fragment_count": ambiguous[identifier],
               "covered_bases_full": covered_full,
               "covered_bases_assessable": covered_assessable, "mean_depth_full": sum(cov)/len(cov),
               "mean_depth_assessable": sum_assessable/assessable_length if assessable_length else None,
               "breadth_full": covered_full/len(cov), "breadth_assessable": covered_assessable/assessable_length if assessable_length else None,
               "fragment_rpkm_all_clean": assigned[identifier]*1e9/(len(cov)*total) if total else None,
               "clean_fragments_denominator": total, "measurement_status": "assessed" if assessed else "not_assessed",
               "reason": "" if assessed else "zero_clean_fragments_or_quantification_ineligible",
               "measurement_target_level": "votu_representative", "resolution_status": ("not_resolvable_ambiguous_only" if not assigned[identifier] else "resolved_and_ambiguous_signal") if ambiguous[identifier] else "reference_limited",
               "unique_mappability": mapfraction, "mappability_status": mapstatus,
               "informative_length":len(info),"covered_bases_informative":info_covered,"breadth_informative":info_covered/len(info) if info else None,
               "mean_depth_informative":sum(cov[pos] for pos in info)/len(info) if info else None,"informative_status":infostatus,
               "competition_evidence": manifest["competition_evidence"],
               "sample_in_discovery_set": (units["biological_sample_id"] in discovery_ids) if manifest.get("discovery_provenance_status") == "provided" else "unknown"}
        if manifest.get('measurement_contract')=='research_quantification_group_v1':
            row.update({k:reference.get(k,'') for k in RESEARCH_FIELDS})
            row['measurement_target_level']='quantification_group'
        rows.append(row)
        detection_rows.extend(apply_detection(row, detections))
    aggregated = []
    host_ids = sorted({r["host_id"] for r in host_rows})
    for host in host_ids:
        group = [r for r in host_rows if r["host_id"] == host]
        length = sum(r["assessable_length"] for r in group)
        covered = sum(r["covered_bases"] for r in group)
        eligible = units["eligibility"]["host_adjustment"]
        status = "assessed" if length and total and units["material_type"] == "bulk_metagenome" else "not_assessed"
        aggregated.append({**identity, "host_id": host, "assessable_length": length, "covered_bases": covered,
                           "mean_backbone_depth": sum(r["depth_sum"] for r in group)/length if length else None,
                           "backbone_breadth": covered/length if length else None,
                           "assigned_fragment_count": sum(r["assigned_fragment_count"] for r in group),
                           "host_adjustment_eligibility": eligible["status"], "status": status,
                           "reason": eligible.get("reason", "") or ("not_bulk_or_no_assessable_backbone" if status != "assessed" else "")})
    write_table(out / "measurement_long.tsv", MEASUREMENT_FIELDS, rows)
    write_table(out / "detection_long.tsv", DETECTION_FIELDS, detection_rows)
    write_table(out / "host_backbone.tsv", HOST_FIELDS, aggregated)
    write_json(out / "measurement_manifest.json", {"schema_version": 1, "unit": units, "profile": profile,
               "reference_version": manifest["mapping_reference_version"], "fragment_filter_counts": dict(counts),
               "input_fragments": total, "observed_sam_fragments": groups, "read_lengths": dict(read_lengths),
               "assigned_fragments": sum(assigned.values()), "raw_reference_counts": dict(assigned),
               "mappability_sha256":sha256(mappability_file) if mappability_file else None,
               "pair_span_counts":dict(pair_spans),"real_data_validated": False, "scientifically_calibrated": False})
    return rows
