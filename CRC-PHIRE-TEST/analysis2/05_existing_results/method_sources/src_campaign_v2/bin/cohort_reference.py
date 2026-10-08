"""Catalog/host adapter. Preserve candidate provenance; mapping has its own IDs."""
from __future__ import annotations
import hashlib
from pathlib import Path
from cohort_common import (child_path, digest, fasta, guard_output, read_json, require,
                           sha256, table, validate_clean_manifest, validate_profile,
                           write_fasta, write_json, write_table)

def _ranges(sequence):
    result = []
    start = None
    for i, base in enumerate(sequence + "N"):
        if base in "ACGT" and start is None:
            start = i
        elif base not in "ACGT" and start is not None:
            result.append([start, i])
            start = None
    return result

def _scope_ok(row, scope):
    if scope == "combined_exploratory":
        return True
    if scope == "virome_discovery":
        if row.get("source_kind") in {"legacy_bacterial_genome", "bacterial_genome"}:
            return False
        return (row.get("catalog_scope") == scope
                or row.get("source_kind") in {"metagenome_assembly", "virome_assembly", "viral_assembly", "independent_fasta"}
                or row.get("material_type") in {"bulk_metagenome", "vlp"}
                or row.get("source_type") in {"metagenome", "viral_contig", "virome"})
    # Source identity, host QC and viral eligibility are separate requirements.
    # A species/reference-panel analysis is a later B2 eligibility decision.
    host_qc = str(row.get("source_genome_eligibility", row.get("main_catalog_eligible", "unknown"))).lower()
    source_valid = (row.get("source_kind") in {"bacterial_genome", "legacy_bacterial_genome"}
                    and row.get("candidate_type") == "provirus_locus"
                    and bool(row.get("source_host_genome"))
                    and row.get("source_host_genome") not in {"unknown", "NA"}
                    and row.get("integration_evidence") not in {None, "", "unknown", "uncertain"}
                    and row.get("slice_verification") in {"verified_this_run", "verified"}
                    and host_qc in {"eligible", "true"})
    if not source_valid:
        return False
    viral_status = row.get("viral_catalog_eligibility")
    if viral_status is None:
        # Legacy v0.2 representative membership alone is not V1 viral quality.
        # Reuse the shared policy instead of creating a divergent cutoff copy.
        from viral_catalog import eligibility
        viral_status, _ = eligibility({**row, "virus_scope": row.get("virus_scope") or "unknown",
                                       "length": row.get("length", 0), "completeness": row.get("completeness"),
                                       "contamination": row.get("contamination")})
    return viral_status == "eligible"

def prepare_reference(catalog_dir, outdir, scope, host_references=None, decoy_fasta=None):
    require(scope in {"host_resolved_prophage", "virome_discovery", "combined_exploratory"}, "Explicit catalog_scope required")
    catalog = Path(catalog_dir).resolve()
    inputs = [catalog]
    if host_references:
        inputs.append(Path(host_references))
    if decoy_fasta:
        inputs.append(Path(decoy_fasta))
    out = guard_output(outdir, inputs)
    source_manifest = read_json(catalog / "catalog_manifest.json")
    version = source_manifest.get("catalog_version")
    require(version, "Missing catalog_version")
    reps_file = catalog / "vOTU_representatives.tsv"
    seq_file = catalog / "vOTU_representatives.fna"
    sequences = fasta(seq_file)
    representatives = table(reps_file)
    master_file = catalog / "viral_sequence_master.tsv"
    if not master_file.exists():
        master_file = catalog / "prophage_master.tsv"
    master_rows = table(master_file)
    occurrence_file = catalog / "source_occurrences.tsv"
    occurrences = table(occurrence_file) if occurrence_file.exists() else master_rows
    masters = {}
    for row in master_rows:
        key = row.get("candidate_id") or row.get("viral_sequence_id")
        require(key and key not in masters, "Duplicate/missing candidate master ID")
        masters[key] = row
    records = {}
    rows = []
    excluded = []
    aliases = set()
    for rep in sorted(representatives, key=lambda r: r["votu_id"]):
        candidate = rep.get("representative_id") or rep.get("candidate_id") or rep.get("viral_sequence_id")
        require(candidate in sequences and candidate in masters, f"Representative header is not mapped: {candidate}")
        require(rep.get("catalog_version", version) == version, "Mixed catalog versions")
        require(rep["votu_id"] not in aliases, "Multiple representatives for one vOTU")
        aliases.add(rep["votu_id"])
        master = masters[candidate]
        source_rows = [row for row in occurrences if row.get("viral_sequence_id") == master.get("viral_sequence_id")] if master.get("viral_sequence_id") else [master]
        in_scope = _scope_ok(master, scope)
        if scope == "host_resolved_prophage" and source_rows:
            in_scope = any(_scope_ok({**master, **row}, scope)
                           and str(row.get("host_resolved_prophage", "true")).lower() == "true" for row in source_rows)
        if scope == "virome_discovery" and occurrence_file.exists():
            in_scope = any(_scope_ok(row, scope) for row in source_rows)
        if not in_scope:
            excluded.append({"candidate_id": candidate, "votu_id": rep["votu_id"], "reason": "host_source_QC_viral_or_slice_eligibility_not_confirmed" if scope == "host_resolved_prophage" else "outside_requested_scope"})
            continue
        sequence = sequences[candidate]
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        if master.get("sequence_sha256"):
            require(master["sequence_sha256"] == sequence_hash, f"Representative hash mismatch: {candidate}")
        identifier = "target_" + digest([version, rep["votu_id"], sequence_hash])[:20]
        records[identifier] = sequence
        rows.append({"reference_id": identifier, "reference_kind": "target", "catalog_version": version,
                     "votu_id": rep["votu_id"], "candidate_id": candidate,
                     "viral_sequence_id": master.get("viral_sequence_id", "seq_" + sequence_hash),
                     "source_host_genome": master.get("source_host_genome", "unknown"),
                     "candidate_type": master.get("candidate_type", "unknown"),
                     "virus_scope": master.get("virus_scope") or "unknown",
                     "length": len(sequence), "sequence_sha256": sequence_hash,
                     "assessable_intervals": _ranges(sequence), "masked_bases": 0})
    mask_rows = []
    host_evidence = []
    if host_references:
        for host in table(host_references):
            require(host.get("host_id") and host.get("fasta"), "Host TSV needs host_id,fasta,mask_bed")
            path = child_path(Path(host_references).parent, host["fasta"])
            source = fasta(path)
            masks = {key: [] for key in source}
            # Every known locus for this host is masked, including candidates outside the chosen scope.
            for row in occurrences:
                if row.get("source_host_genome", row.get("genome_id")) != host["host_id"]:
                    continue
                key = row.get("original_contig_id") if row.get("original_contig_id") in source else row.get("contig_id")
                if row.get("start0") not in {None, "", "NA"}:
                    require(key in source, f"Host locus contig missing from supplied host reference: {host['host_id']}")
                    start, end = int(row["start0"]), int(row["end0"])
                    require(0 <= start < end <= len(source[key]), "Invalid host locus mask coordinates")
                    locus = source[key][start:end]
                    if row.get("strand") == "-":
                        locus = locus.translate(str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN"))[::-1]
                    expected = row.get("sequence_sha256")
                    require(expected, "Host catalog locus lacks sequence hash; provide a traceable catalog")
                    require(hashlib.sha256(locus.encode()).hexdigest() == expected,
                            "Host locus sequence hash mismatch; supplied host reference differs from source")
                    masks[key].append((start, end, "catalog_locus_sequence_verified"))
            if host.get("mask_bed") and host["mask_bed"] not in {"NA", "unknown"}:
                bed = child_path(Path(host_references).parent, host["mask_bed"])
                for line in bed.read_text().splitlines():
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    require(len(parts) >= 3 and parts[0] in source, "Mask BED contig mismatch")
                    masks[parts[0]].append((int(parts[1]), int(parts[2]), "explicit_mask"))
            for name, sequence in source.items():
                identifier = "host_" + digest([host["host_id"], name, sha256(path)])[:20]
                require(identifier not in records, "Duplicate host reference")
                bases = list(sequence)
                for start, end, reason in masks[name]:
                    require(0 <= start < end <= len(sequence), "Invalid host mask interval")
                    bases[start:end] = "N" * (end-start)
                    mask_rows.append({"reference_id": identifier, "start0": start, "end0": end, "reason": reason,
                                      "host_id": host["host_id"], "original_contig_id": name})
                sequence = "".join(bases)
                records[identifier] = sequence
                rows.append({"reference_id": identifier, "reference_kind": "host", "host_id": host["host_id"],
                             "original_contig_id": name, "length": len(sequence),
                             "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                             "assessable_intervals": _ranges(sequence), "masked_bases": sequence.count("N")})
            host_evidence.append({"host_id": host["host_id"], "input_fasta": str(path), "sha256": sha256(path),
                                  "mask_basis": "catalog_loci_and_explicit_bed", "mask_completeness": "not_exhaustive"})
    if decoy_fasta:
        for name, sequence in fasta(decoy_fasta).items():
            identifier = "decoy_" + digest([name, sequence])[:20]
            records[identifier] = sequence
            rows.append({"reference_id": identifier, "reference_kind": "decoy", "original_contig_id": name,
                         "length": len(sequence), "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                         "assessable_intervals": _ranges(sequence), "masked_bases": 0})
    manifest = {"schema_version": 1, "catalog_version": version, "catalog_scope": scope,
                "catalog_manifest_sha256": sha256(catalog / "catalog_manifest.json"),
                "representatives_sha256": sha256(seq_file), "representative_table_sha256": sha256(reps_file),
                "source_occurrences_sha256": sha256(occurrence_file) if occurrence_file.exists() else None,
                "source_master_sha256": sha256(master_file), "representative_strategy": "one_explicit_representative_per_votu",
                "discovery_samples": source_manifest.get("discovery_samples", []),
                "discovery_provenance_status": "provided" if "discovery_samples" in source_manifest else "unknown_legacy_catalog",
                "specificity_scope": "supplied_references_only",
                "host_scope_contract": "bacterial_source_provirus_integration_verified_slice_host_QC_eligible_viral_eligible_no_species_panel",
                "competition_evidence": "supplied_host_or_decoy" if host_evidence or decoy_fasta else "insufficient_no_host_or_decoy",
                "host_inputs": host_evidence, "references": rows,
                "decoy_sha256": sha256(decoy_fasta) if decoy_fasta else None}
    manifest["mapping_reference_version"] = "mapref_" + digest(manifest)[:20]
    write_fasta(out / "reference.fna", records)
    manifest["reference_fasta_sha256"] = sha256(out / "reference.fna")
    write_json(out / "reference_manifest.json", manifest)
    write_table(out / "reference_map.tsv", ["reference_id", "reference_kind", "catalog_version", "votu_id", "candidate_id",
                "viral_sequence_id", "source_host_genome", "candidate_type", "virus_scope", "host_id", "original_contig_id",
                "length", "sequence_sha256", "masked_bases"], rows)
    write_table(out / "host_mask.bed", ["reference_id", "start0", "end0", "reason", "host_id", "original_contig_id"], mask_rows)
    write_table(out / "scope_exclusions.tsv", ["candidate_id", "votu_id", "reason"], excluded)
    return manifest

def prepare_inputs(clean_reads_manifest, catalog_dir, outdir, scope, measurement_profile, detection_profiles,
                   host_references=None, decoy_fasta=None, research_reference_config=None):
    manifest = validate_clean_manifest(read_json(clean_reads_manifest))
    profile = validate_profile(read_json(measurement_profile))
    out = guard_output(outdir, [clean_reads_manifest, catalog_dir])
    if research_reference_config:
        require(not host_references and not decoy_fasta,'Research config owns host core and decoys')
        from cohort_research_reference import prepare
        reference=prepare(catalog_dir,out/'reference',research_reference_config)
    else:
        reference = prepare_reference(catalog_dir, out / "reference", scope, host_references, decoy_fasta)
    (out / "units").mkdir(exist_ok=True)
    for unit in manifest["units"]:
        write_json(out / "units" / (unit["unit_id"] + ".json"), unit)
    write_json(out / "clean_reads_manifest.json", manifest)
    write_json(out / "measurement_profile.json", profile)
    write_json(out / "detection_profiles.json", read_json(detection_profiles))
    receipt = {"schema_version": 1, "preprocessing_id": manifest["preprocessing_id"],
               "mapping_reference_version": reference["mapping_reference_version"],
               "clean_reads_manifest_sha256": sha256(clean_reads_manifest), "measurement_profile_sha256": sha256(measurement_profile),
               "detection_profiles_sha256": sha256(detection_profiles),
               "units": len(manifest["units"]), "real_data_validated": False, "scientifically_calibrated": False}
    write_json(out / "prepared_manifest.json", receipt)
    return receipt
