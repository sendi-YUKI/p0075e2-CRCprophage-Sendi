#!/usr/bin/env python3
"""Fetch the pinned A0/A1 public assets with resumable, verified local manifests."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import sys
import tarfile
import time
import urllib.request

CHUNK = 4 * 1024 * 1024
RESERVE = 2 * 1024**3


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temp.replace(path)


def hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            sha.update(block)
            md5.update(block)
    return {"sha256": sha.hexdigest(), "md5": md5.hexdigest(), "bytes": path.stat().st_size}


def space(path, needed):
    available = shutil.disk_usage(path).free
    if available < needed + RESERVE:
        raise RuntimeError(f"Insufficient disk at {path}: {available} bytes free, need {needed} plus {RESERVE} reserve")


def validate_archive_hashes(path, asset, previous=None):
    result = hashes(path)
    if asset.get("download_bytes") and result["bytes"] != asset["download_bytes"]:
        raise RuntimeError(f"Size mismatch for {path}: {result['bytes']} != {asset['download_bytes']}")
    for name in ("md5", "sha256"):
        expected = asset.get("expected_" + name)
        if expected and result[name] != expected:
            raise RuntimeError(f"Publisher {name} mismatch for {path}; preserve file for inspection")
    if previous and previous.get("sha256") and result["sha256"] != previous["sha256"]:
        raise RuntimeError(f"Local locked SHA256 mismatch for {path}; preserve file for inspection")
    return result


def download(asset, destination, previous=None):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file():
        return validate_archive_hashes(destination, asset, previous)
    part = destination.with_name(destination.name + ".part")
    expected_size = asset.get("download_bytes")
    offset = part.stat().st_size if part.exists() else 0
    if expected_size and offset == expected_size:
        result = validate_archive_hashes(part, asset, previous)
        part.replace(destination)
        return result
    if expected_size and offset > expected_size:
        raise RuntimeError(f"Partial file larger than pinned asset: {part}; preserve for inspection")
    space(destination.parent, max(0, (expected_size or 8 * 1024**3) - offset))
    for attempt in range(1, 4):
        offset = part.stat().st_size if part.exists() else 0
        request = urllib.request.Request(asset["url"], headers={"User-Agent": "crc-phage-pipeline/0.1", "Accept-Encoding": "identity"})
        if offset:
            request.add_header("Range", f"bytes={offset}-")
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                if offset and response.status == 206:
                    if not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                        raise RuntimeError("Server returned an incorrect resumed byte range")
                    mode = "ab"
                elif response.status == 200:
                    mode, offset = "wb", 0
                else:
                    raise RuntimeError(f"Unexpected download HTTP status {response.status}")
                print(f"Downloading {asset['id']} from byte {offset}", flush=True)
                last_report = time.monotonic()
                with part.open(mode) as handle:
                    for block in iter(lambda: response.read(CHUNK), b""):
                        space(destination.parent, len(block))
                        handle.write(block)
                        offset += len(block)
                        if time.monotonic() - last_report >= 20:
                            print(f"{asset['id']}: {offset}/{expected_size or 'unknown'} bytes", flush=True)
                            last_report = time.monotonic()
            result = validate_archive_hashes(part, asset, previous)
            part.replace(destination)
            return result
        except (OSError, TimeoutError) as exc:
            if attempt == 3:
                raise
            print(f"Download attempt {attempt} failed: {exc}; resuming in {attempt * 2}s", file=sys.stderr, flush=True)
            time.sleep(attempt * 2)
    raise AssertionError("unreachable")


def path_within(root, relative):
    pure = PurePosixPath(relative)
    if pure.is_absolute() or any(p in ("..", "") for p in pure.parts) or "\\" in relative or ":" in relative:
        raise RuntimeError(f"Unsafe archive path: {relative}")
    target = root.joinpath(*pure.parts).resolve()
    if not target.is_relative_to(root.resolve()):
        raise RuntimeError(f"Archive path escapes target: {relative}")
    return target


def extract_database(archive, root, asset):
    destination = root / asset["relative_directory"]
    destination.mkdir(parents=True, exist_ok=True)
    database = destination / asset["archive_top_directory"]
    receipt = destination / "extraction_manifest.json"
    if receipt.exists():
        manifest = json.loads(receipt.read_text(encoding="utf-8"))
        if manifest["archive_sha256"] != hashes(archive)["sha256"]:
            raise RuntimeError(f"Existing extraction belongs to a different archive: {destination}")
        for entry in manifest["files"]:
            entry_path = path_within(destination, entry["path"])
            if not entry_path.is_file() or hashes(entry_path)["sha256"] != entry["sha256"]:
                raise RuntimeError(f"Extracted database file changed or missing: {entry_path}")
        return database, manifest
    # Never create/follow filesystem symlinks. Upstream geNomad has internal links;
    # resolve them within this archive and materialize verified regular-file copies.
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        member_by_name = {PurePosixPath(member.name).as_posix(): member for member in members}
        if len(member_by_name) != len(members):
            raise RuntimeError("Duplicate member paths in database archive")

        def regular_source(member, seen=None):
            seen = set() if seen is None else seen
            if member.name in seen:
                raise RuntimeError(f"Cyclic archive link: {member.name}")
            seen.add(member.name)
            if member.isfile():
                return member
            if not (member.issym() or member.islnk()):
                raise RuntimeError(f"Archive link does not lead to a regular file: {member.name}")
            link = PurePosixPath(member.linkname)
            if link.is_absolute() or ".." in link.parts or "\\" in member.linkname or ":" in member.linkname:
                raise RuntimeError(f"Unsafe archive link: {member.name} -> {member.linkname}")
            target_name = (PurePosixPath(member.name).parent / link if member.issym() else link).as_posix()
            path_within(destination, target_name)
            if PurePosixPath(target_name).parts[0] != asset["archive_top_directory"]:
                raise RuntimeError(f"Archive link leaves database root: {member.name}")
            if target_name not in member_by_name:
                raise RuntimeError(f"Archive link target is absent: {member.name} -> {target_name}")
            return regular_source(member_by_name[target_name], seen)

        link_sources = {}
        for member in members:
            path_within(destination, member.name)
            if not (member.isdir() or member.isfile() or member.issym() or member.islnk()):
                raise RuntimeError(f"Device/special file is not allowed in database archive: {member.name}")
            if PurePosixPath(member.name).parts[0] != asset["archive_top_directory"]:
                raise RuntimeError(f"Unexpected archive root: {member.name}")
            if member.issym() or member.islnk():
                link_sources[member.name] = regular_source(member)
        expanded = sum(member.size for member in members if member.isfile()) + sum(member.size for member in link_sources.values())
        space(destination, expanded)
        for member in members:
            target = path_within(destination, member.name)
            cursor = destination
            for component in PurePosixPath(member.name).parts:
                cursor = cursor / component
                if cursor.is_symlink():
                    raise RuntimeError(f"Refusing to overwrite symlink {cursor}")
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                with handle.extractfile(member) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output, CHUNK)
                os.chmod(target, 0o644)
        for name, source_member in link_sources.items():
            target = path_within(destination, name)
            source = path_within(destination, source_member.name)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            os.chmod(target, 0o644)
    for required in asset["required_paths"]:
        if not (database / required).exists():
            raise RuntimeError(f"Database required path missing: {database / required}")
    if asset["id"] == "genomad":
        version = (database / "version.txt").read_text(encoding="utf-8").strip()
        if version != asset["version"]:
            raise RuntimeError(f"Unexpected geNomad database version: {version}")
    files = [{"path": path.relative_to(destination).as_posix(), **hashes(path)} for path in sorted(database.rglob("*")) if path.is_file()]
    manifest = {"archive_sha256": hashes(archive)["sha256"], "extracted_utc": now(), "expanded_bytes": expanded, "archive_links_materialized": {name: member.name for name, member in link_sources.items()}, "files": files}
    save_json(receipt, manifest)
    return database, manifest


def extract_genome(archive, root, asset):
    destination = root / asset["fasta_filename"]
    temporary = destination.with_name(destination.name + ".part")
    space(root, asset["expanded_estimate_bytes"])
    with gzip.open(archive, "rb") as source, temporary.open("wb") as output:
        shutil.copyfileobj(source, output, CHUNK)
    identifiers, bases = [], 0
    with temporary.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split()[0])
            else:
                bases += len(line.strip())
    if identifiers != asset["sequence_accessions"] or bases != asset["expected_bases"]:
        raise RuntimeError(f"Reference identity mismatch for {asset['id']}: {identifiers}, {bases} bases")
    temporary.replace(destination)
    return destination, {"sequence_accessions": identifiers, "bases": bases, **hashes(destination)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--kind", choices=("genomes", "databases"), required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--assets", type=Path, default=Path(__file__).resolve().parents[1] / "assets" / "required_assets.json")
    args = parser.parse_args()
    catalog = json.loads(args.assets.read_text(encoding="utf-8"))
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    previous = json.loads(args.manifest.read_text(encoding="utf-8")) if args.manifest.exists() else {}
    if previous and (previous.get("kind") != args.kind):
        raise RuntimeError("Manifest kind differs; use a separate manifest per asset kind")
    manifest = {"schema_version": 1, "kind": args.kind, "updated_utc": now(), "root": str(root), "asset_catalog_sha256": hashes(args.assets)["sha256"], "assets": previous.get("assets", {})}
    save_json(args.manifest, manifest)
    sample_rows = []
    try:
        for asset in catalog[args.kind]:
            old = manifest["assets"].get(asset["id"], {})
            entry = {"id": asset["id"], "version": asset["version"], "url": asset["url"], "status": "downloading", "started_utc": now(), "archive": old.get("archive", {})}
            manifest["assets"][asset["id"]] = entry
            save_json(args.manifest, manifest)
            archive = root / "archives" / asset["filename"]
            # Refuse before a large download if the extraction allowance cannot fit.
            space(root, asset["expanded_estimate_bytes"] + max(0, asset["download_bytes"] - (archive.stat().st_size if archive.exists() else 0)))
            archive_hashes = download(asset, archive, old.get("archive"))
            entry["archive"] = {"path": str(archive), **archive_hashes}
            entry["status"] = "extracting"
            save_json(args.manifest, manifest)
            if args.kind == "databases":
                result, extraction = extract_database(archive, root, asset)
                entry.update({"path": str(result), "expanded_bytes": extraction["expanded_bytes"], "files": extraction["files"], "requires_generated_index": asset["requires_generated_index"], "status": "downloaded_index_pending" if asset["requires_generated_index"] else "ready"})
            else:
                result, sequence = extract_genome(archive, root, asset)
                entry.update({"path": str(result), "sequence": sequence, "status": "ready", "evidence_role": asset["evidence_role"], "evidence": asset["evidence"]})
                sample_rows.append({"genome_id": asset["id"], "fasta": result.relative_to(root).as_posix()})
            entry["completed_utc"] = now()
            save_json(args.manifest, manifest)
            print(f"{asset['id']}: {entry['status']} at {entry['path']}", flush=True)
        if args.kind == "genomes":
            with (root / "genomes.tsv").open("w", newline="", encoding="utf-8") as output:
                writer = csv.DictWriter(output, fieldnames=("genome_id", "fasta"), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(sample_rows)
        manifest["download_exit_code"] = 0
        manifest["updated_utc"] = now()
        save_json(args.manifest, manifest)
        return 0
    except Exception as exc:
        if "entry" in locals():
            entry.update({"status": "failed", "error": str(exc), "failed_utc": now()})
        manifest.update({"download_exit_code": 1, "updated_utc": now()})
        save_json(args.manifest, manifest)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
