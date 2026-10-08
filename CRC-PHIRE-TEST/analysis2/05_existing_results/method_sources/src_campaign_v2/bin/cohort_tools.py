"""External B1 tool adapters: explicit argv, checked exits, no automatic DB download."""
from __future__ import annotations
from collections import defaultdict
from itertools import zip_longest
from pathlib import Path
import shutil
import subprocess
from cohort_common import (digest, fasta, read_json, require, sha256, text_open,
                           validate_profile, write_fasta, write_json, write_table)
from cohort_measure import alignment_metrics, fill_secondary_sequence, passes_sequence, sam_groups

def run_checked(argv, logfile, cwd=None):
    with Path(logfile).open("w") as log:
        subprocess.run([str(x) for x in argv], stdout=log, stderr=subprocess.STDOUT, cwd=cwd, check=True)

def version(argv):
    result = subprocess.run(argv, capture_output=True, text=True, check=True)
    return (result.stdout + result.stderr).strip()

def build_index(reference_dir, outdir, threads=3):
    source = Path(reference_dir)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    manifest = read_json(source / "reference_manifest.json")
    require(sha256(source / "reference.fna") == manifest["reference_fasta_sha256"], "Mapping FASTA checksum mismatch")
    shutil.copy2(source / "reference.fna", out / "reference.fna")
    shutil.copy2(source / "reference_manifest.json", out / "reference_manifest.json")
    for name in ["reference_map.tsv", "host_mask.bed", "scope_exclusions.tsv"]:
        shutil.copy2(source / name, out / name)
    if (source/'quantification_group_members.tsv').is_file():
        shutil.copy2(source/'quantification_group_members.tsv',out/'quantification_group_members.tsv')
    command = ["bowtie2-build", "--threads", str(threads), str(out / "reference.fna"), str(out / "index")]
    if manifest["references"]:
        run_checked(command, out / "bowtie2-build.log")
        require(len(list(out.glob("index*.bt2*"))) == 6, "Incomplete Bowtie2 index")
    write_json(out / "index_manifest.json", {"command": command, "status": "built" if manifest["references"] else "empty_reference",
               "reference_version": manifest["mapping_reference_version"], "bowtie2": version(["bowtie2", "--version"]),
               "index_files": {p.name: sha256(p) for p in sorted(out.glob("index*.bt2*"))}})

def fastq_records(path):
    with text_open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq, plus, qual = handle.readline().strip(), handle.readline(), handle.readline().strip()
            require(header.startswith("@") and plus.startswith("+") and seq and len(seq) == len(qual), "Invalid FASTQ")
            name = header[1:].split()[0]
            if name.endswith("/1") or name.endswith("/2"):
                name = name[:-2]
            yield name, seq, qual

def empty_alignment(unit, read1, read2, path):
    with Path(path).open("w") as handle:
        handle.write("@HD\tVN:1.6\tSO:queryname\n")
        first = fastq_records(read1)
        if unit["layout"] == "PE":
            for a, b in zip_longest(first, fastq_records(read2)):
                require(a and b and a[0] == b[0], "PE reads not synchronized")
                for record, flag in [(a, 77), (b, 141)]:
                    handle.write(f"{record[0]}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{record[1]}\t{record[2]}\n")
        else:
            for name, seq, qual in first:
                handle.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\n")

def map_reads(index_dir, unit_file, profile_file, read1, read2, outdir, threads=3):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    unit, profile = read_json(unit_file), validate_profile(read_json(profile_file))
    for key, actual in [("reads_1", read1), ("reads_2", read2)]:
        if actual:
            require(unit.get("checksums", {}).get(unit[key]) == sha256(actual), "Clean read checksum changed after preparation")
    reference = read_json(Path(index_dir) / "reference_manifest.json")
    sam = out / "alignments.sam"
    command = ["bowtie2", "--end-to-end", "--very-sensitive", "-k", str(profile["report_k"]),
               "--seed", "0", "-p", str(threads), "-x", str(Path(index_dir) / "index"), "-S", str(sam)]
    if unit["layout"] == "PE":
        require(read2, "PE second file missing")
        command += ["--no-mixed", "--no-discordant", "--"+profile["pair_orientation"],
                    "-I", str(profile["insert_min"]), "-X", str(profile["insert_max"]), "-1", str(read1), "-2", str(read2)]
    else:
        command += ["-U", str(read1)]
    primary_sam = out / "alignments.primary.sam"
    primary_command = list(command)
    report_at = primary_command.index("-k")
    del primary_command[report_at:report_at+2]
    primary_command[primary_command.index("-S")+1] = str(primary_sam)
    if reference["references"] and unit["counts"]["fragments"]:
        run_checked(primary_command, out / "bowtie2-primary.log")
        run_checked(command, out / "bowtie2-audit.log")
    else:
        empty_alignment(unit, read1, read2, sam)
        shutil.copy2(sam, primary_sam)
    bam = out / "alignments.qname.bam"
    sort = ["samtools", "sort", "-n", "-@", str(max(0, threads-1)), "-m", "256M", "-o", str(bam), str(sam)]
    run_checked(sort, out / "samtools-sort.log")
    sam.unlink()
    primary_bam = out / "alignments.primary.qname.bam"
    primary_sort = ["samtools", "sort", "-n", "-@", str(max(0, threads-1)), "-m", "256M", "-o", str(primary_bam), str(primary_sam)]
    run_checked(primary_sort, out / "samtools-primary-sort.log")
    primary_sam.unlink()
    write_json(out / "mapping_manifest.json", {"schema_version": 1, "unit_id": unit["unit_id"], "command": command,
               "primary_command": primary_command, "primary_bam_sha256": sha256(primary_bam),
               "mapq_source": profile["mapq_source"], "sort_command": sort, "alignment_report_policy": "up_to_k_per_read_or_pair_secondary_retained_for_audit",
               "report_limit_is_exhaustive": False, "bam_sha256": sha256(bam),
               "bowtie2": version(["bowtie2", "--version"]), "samtools": version(["samtools", "--version"]),
               "clean_checksums": unit["checksums"], "profile": profile})
    return bam

def mappability(index_dir, profile_file, outdir, threads=3):
    if read_json(profile_file).get('informative_policy'):
        from cohort_informative import build
        return build(index_dir,profile_file,outdir,threads)['references']
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    profile = validate_profile(read_json(profile_file))
    references = fasta(Path(index_dir) / "reference.fna")
    manifest = read_json(Path(index_dir) / "reference_manifest.json")
    targets = [r for r in manifest["references"] if r["reference_kind"] == "target"]
    reads, tile_meta = {}, {}
    potential = {}
    for target in targets:
        name, length = target["reference_id"], target["length"]
        read_length = profile["mappability_read_length"]
        starts = list(range(0, max(0, length-read_length+1), profile["mappability_stride"]))
        if length >= read_length:
            starts = sorted(set(starts + [length-read_length]))
        starts = [p for p in starts if set(references[name][p:p+read_length]) <= set("ACGT")]
        potential[name] = len(starts)
        maximum = profile["mappability_max_tiles_per_reference"]
        if len(starts) > maximum:
            indices = [round(i*(len(starts)-1)/(maximum-1)) for i in range(maximum)] if maximum > 1 else [len(starts)//2]
            starts = [starts[i] for i in indices]
        for start in starts:
            tile_id = "tile_" + digest([name, start, read_length])[:24]
            reads[tile_id] = references[name][start:start+read_length]
            tile_meta[tile_id] = {"reference_id": name, "start0": start, "end0": start+read_length}
    write_fasta(out / "tiles.fna", reads)
    command = ["bowtie2", "--end-to-end", "--very-sensitive", "-f", "-k", str(profile["report_k"]),
               "--seed", "0", "-p", str(threads), "-x", str(Path(index_dir) / "index"),
               "-U", str(out / "tiles.fna"), "-S", str(out / "tiles.sam")]
    evidence = []
    if reads:
        run_checked(command, out / "tiles-mapping.log")
        with (out / "tiles.sam").open() as handle:
            for name, records in sam_groups(handle):
                require(name in tile_meta, "Unexpected mappability query")
                primary = {r["mate"]: r for r in records if not r["flag"] & (256 | 2048)}
                locations = set()
                for raw in records:
                    if raw["flag"] & (4 | 2048):
                        continue
                    record = fill_secondary_sequence(raw, primary)
                    if record["qual"] == "*":  # Bowtie2 FASTA input: known synthetic quality for geometric audit
                        record["qual"] = "I"*len(record["seq"])
                    metrics = alignment_metrics(record, references, profile)
                    if passes_sequence(metrics, profile):
                        locations.add(metrics["location"])
                origin = tile_meta[name]
                self_hit = any(loc[0] == origin["reference_id"] and loc[1] == origin["start0"] for loc in locations)
                unique = len(locations) == 1 and self_hit and len(records) < profile["report_k"]
                evidence.append({"tile_id": name, **origin, "unique": unique, "self_hit": self_hit,
                                 "n_passing_locations": len(locations), "reporting_saturated": len(records) >= profile["report_k"]})
    by_ref = defaultdict(list)
    for row in evidence:
        by_ref[row["reference_id"]].append(row)
    rows = []
    for target in targets:
        name = target["reference_id"]
        results = by_ref[name]
        rows.append({"reference_id": name, "status": "estimated" if results else "not_assessed",
                     "unique_tile_fraction": sum(r["unique"] for r in results)/len(results) if results else None,
                     "n_tiles": len(results), "potential_stride_tiles": potential[name],
                     "uniform_subsampled": len(results) < potential[name],
                     "read_length": profile["mappability_read_length"]})
    write_table(out / "tile_evidence.tsv", ["tile_id", "reference_id", "start0", "end0", "unique", "self_hit",
                "n_passing_locations", "reporting_saturated"], evidence)
    write_json(out / "mappability.json", {"schema_version": 1, "algorithm": "synthetic_se_tiles_bowtie2_k_v1",
               "measurement_profile": profile, "mapping_reference_version": manifest["mapping_reference_version"],
               "scope": "supplied_references_only_no_unknown_decoy_guarantee", "command": command,
               "paired_end_mappability": "not_assessed", "references": rows})
    return rows

def metaphlan_preflight(database, index, manifest_file=None):
    require(index and index != "latest", "Fixed MetaPhlAn index required")
    db = Path(database)
    require(db.is_dir(), "MetaPhlAn database missing")
    require((db / (index + ".pkl")).is_file(), "MetaPhlAn metadata PKL missing")
    extensions = ["bt2l", "bt2"]
    complete = any(all((db / (index + "." + part + "." + ext)).is_file()
                      for part in ["1", "2", "3", "4", "rev.1", "rev.2"]) for ext in extensions)
    require(complete, "MetaPhlAn complete fixed Bowtie2 database required; automatic download disabled")
    require(manifest_file and Path(manifest_file).is_file(), "MetaPhlAn database manifest required")
    manifest = read_json(manifest_file)
    require(manifest.get("index") == index and manifest.get("database_files"), "Database manifest index/files missing")
    for name, expected in manifest["database_files"].items():
        require((db / name).is_file() and sha256(db / name) == expected, "MetaPhlAn database hash mismatch")
    return manifest

def parse_metaphlan(path, unit):
    fields = None
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#clade_name"):
            fields = line.lstrip("#").split("\t")
        elif line and not line.startswith("#"):
            require(fields and "relative_abundance" in fields, "Unsupported MetaPhlAn table")
            cells = dict(zip(fields, line.split("\t"), strict=True))
            clade = cells["clade_name"]
            if clade.split("|")[-1].startswith("s__") or clade.upper() == "UNCLASSIFIED":
                value = float(cells["relative_abundance"])
                require(0 <= value <= 100, "MetaPhlAn relative abundance must be percent")
                rows.append({"unit_id": unit["unit_id"], "biological_sample_id": unit["biological_sample_id"],
                             "library_id": unit["library_id"], "clade_name": clade,
                             "taxon_level": "unclassified" if clade.upper() == "UNCLASSIFIED" else "species",
                             "ncbi_tax_id": cells.get("NCBI_tax_id", ""), "relative_abundance_percent": value,
                             "relative_abundance_fraction": value/100, "status": "assessed",
                             "host_adjustment_eligibility": unit["eligibility"]["host_adjustment"]["status"]})
    require(fields is not None, "MetaPhlAn output header missing")
    return rows

METAPHLAN_FIELDS = ["unit_id", "biological_sample_id", "library_id", "clade_name", "taxon_level", "ncbi_tax_id",
                   "relative_abundance_percent", "relative_abundance_fraction", "status", "host_adjustment_eligibility"]

def metaphlan(unit_file, read1, read2, database, index, database_manifest, outdir, threads=3):
    unit = read_json(unit_file)
    require(unit["material_type"] == "bulk_metagenome", "MetaPhlAn host profiling here requires bulk material")
    dbmanifest = metaphlan_preflight(database, index, database_manifest)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    infile = str(read1) + ("," + str(read2) if unit["layout"] == "PE" else "")
    command = ["metaphlan", infile, "--input_type", "fastq", "--db_dir", str(database), "--index", index,
               "--offline", "--nproc", str(threads), "--bt2_ps", "very-sensitive", "--min_mapq_val", "5",
               "--mapout", str(out / "markers.mapout.bz2"), "-t", "rel_ab", "-o", str(out / "profile.tsv")]
    if unit["counts"]["fragments"]:
        run_checked(command, out / "metaphlan.log")
        rows = parse_metaphlan(out / "profile.tsv", unit)
        status = "assessed"
    else:
        rows, status = [], "not_assessed_zero_clean_fragments"
    write_table(out / "metaphlan_species.tsv", METAPHLAN_FIELDS, rows)
    write_json(out / "metaphlan_manifest.json", {"schema_version": 1, "unit_id": unit["unit_id"], "status": status,
               "command": command, "metaphlan": version(["metaphlan", "--version"]), "database_manifest": dbmanifest,
               "species_unit": "relative_abundance_fraction_and_percent_not_depth", "real_data_validated": False})
