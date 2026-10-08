#!/usr/bin/env python3
"""B1 production CLI. Tools run in pinned Nextflow containers; prepare is dependency-free."""
from __future__ import annotations
import argparse
import json
import sys
from cohort_common import ContractError
from cohort_measure import measure
from cohort_reference import prepare_inputs, prepare_reference
from cohort_report import aggregate
from cohort_tools import build_index, map_reads, mappability, metaphlan

def parser():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    for option in ["clean-reads-manifest", "catalog-dir", "outdir", "catalog-scope", "measurement-profile", "detection-profiles"]:
        prep.add_argument("--" + option, required=True)
    prep.add_argument("--research-reference-config")
    prep.add_argument("--host-references")
    prep.add_argument("--decoy-fasta")
    ref = commands.add_parser("reference")
    for option in ["catalog-dir", "outdir", "catalog-scope"]:
        ref.add_argument("--" + option, required=True)
    ref.add_argument("--host-references")
    ref.add_argument("--decoy-fasta")
    index = commands.add_parser("index")
    index.add_argument("--reference-dir", required=True)
    index.add_argument("--outdir", required=True)
    index.add_argument("--threads", type=int, default=3)
    mapping = commands.add_parser("map")
    for option in ["index-dir", "unit", "profile", "reads1", "outdir"]:
        mapping.add_argument("--" + option, required=True)
    mapping.add_argument("--reads2")
    mapping.add_argument("--threads", type=int, default=3)
    counts = commands.add_parser("measure")
    for option in ["sam", "reference-dir", "unit", "profile", "detections", "outdir"]:
        counts.add_argument("--" + option, required=True)
    counts.add_argument("--mappability")
    counts.add_argument("--primary-sam")
    tile = commands.add_parser("mappability")
    for option in ["index-dir", "profile", "outdir"]:
        tile.add_argument("--" + option, required=True)
    tile.add_argument("--threads", type=int, default=3)
    host = commands.add_parser("metaphlan")
    for option in ["unit", "reads1", "database", "index", "database-manifest", "outdir"]:
        host.add_argument("--" + option, required=True)
    host.add_argument("--reads2")
    host.add_argument("--threads", type=int, default=3)
    report = commands.add_parser("report")
    report.add_argument("--prepared", required=True)
    report.add_argument("--measurements", nargs="*", default=[])
    report.add_argument("--metaphlan-results", nargs="*", default=[])
    report.add_argument("--metaphlan-disabled-reason", default="not_enabled")
    report.add_argument("--outdir", required=True)
    return p

def main():
    a = parser().parse_args()
    if a.command == "prepare":
        result = prepare_inputs(a.clean_reads_manifest, a.catalog_dir, a.outdir, a.catalog_scope,
                                a.measurement_profile, a.detection_profiles, a.host_references, a.decoy_fasta, a.research_reference_config)
    elif a.command == "reference":
        result = prepare_reference(a.catalog_dir, a.outdir, a.catalog_scope, a.host_references, a.decoy_fasta)
    elif a.command == "index":
        result = build_index(a.reference_dir, a.outdir, a.threads)
    elif a.command == "map":
        result = str(map_reads(a.index_dir, a.unit, a.profile, a.reads1, a.reads2, a.outdir, a.threads))
    elif a.command == "measure":
        result = measure(a.sam, a.reference_dir, a.unit, a.profile, a.detections, a.outdir, a.mappability, a.primary_sam)
    elif a.command == "mappability":
        result = mappability(a.index_dir, a.profile, a.outdir, a.threads)
    elif a.command == "metaphlan":
        result = metaphlan(a.unit, a.reads1, a.reads2, a.database, a.index, a.database_manifest, a.outdir, a.threads)
    else:
        result = aggregate(a.prepared, a.measurements, a.metaphlan_results, a.outdir, a.metaphlan_disabled_reason)
    print(json.dumps({"command": a.command, "status": "completed", "result_rows": len(result) if isinstance(result, list) else None}))

if __name__ == "__main__":
    try:
        main()
    except (ContractError, OSError, ValueError) as error:
        print(f"B1 ERROR: {error}", file=sys.stderr)
        sys.exit(1)
