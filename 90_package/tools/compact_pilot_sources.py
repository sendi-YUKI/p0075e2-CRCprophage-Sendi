#!/usr/bin/env python3
"""Compact historical pilot source trees without losing any source version.

Examples (all commands refuse to overwrite a nonempty output directory)::

    python compact_pilot_sources.py compact --source 05_existing_results/method_sources --output staging_sources
    python compact_pilot_sources.py verify --archive staging_sources --source 05_existing_results/method_sources
    python compact_pilot_sources.py restore --archive staging_sources --output restored_sources

The readable canonical tree is src_pilot6_measurement. Files identical at the
same relative path in other snapshots point to that tree; every differing or
missing file is stored in source_overlays/<snapshot>/<relative path>. The JSON
manifest records each original file, its storage path, byte count, and SHA256.
The restore command reconstructs all three snapshots, including empty folders,
then verifies their complete file inventories and byte hashes. No input files
are modified. Historical absolute paths inside source files are left intact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import tempfile


CANONICAL = "src_pilot6_measurement"
SNAPSHOTS = (CANONICAL, "src_reads_v1", "src_campaign_v2")
MANIFEST = "source_snapshot_manifest.json"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def inventory(root: Path) -> tuple[dict[str, dict], list[str]]:
    if not root.is_dir():
        raise ValueError(f"Missing snapshot directory: {root}")
    files = {}
    directories = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or (hasattr(path, "is_junction") and path.is_junction()):
            raise ValueError(f"Snapshot contains a link rather than a file/directory: {path}")
        relative = path.relative_to(root).as_posix()
        if path.is_dir():
            directories.append(relative)
        elif path.is_file():
            files[relative] = {"bytes": path.stat().st_size, "sha256": digest(path)}
        else:
            raise ValueError(f"Unsupported snapshot entry: {path}")
    return files, directories


def safe_child(root: Path, relative: str) -> Path:
    """Reject absolute/parent paths and links escaping a manifest's root."""
    posix = PurePosixPath(relative)
    if not relative or "\\" in relative or ":" in relative or posix.is_absolute():
        raise ValueError(f"Invalid manifest path: {relative!r}")
    if any(part in ("", ".", "..") for part in posix.parts):
        raise ValueError(f"Invalid manifest path: {relative!r}")
    path = root.joinpath(*posix.parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Manifest path escapes its root: {relative!r}")
    return path


def prepare_output(output: Path, source: Path) -> None:
    output = output.resolve()
    source = source.resolve()
    if output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError("Input and output directories must not contain one another")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output must be absent or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)


def load_manifest(archive: Path) -> dict:
    manifest = json.loads((archive / MANIFEST).read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("canonical_snapshot") != CANONICAL:
        raise ValueError("Unsupported source snapshot manifest")
    if set(manifest.get("snapshots", {})) != set(SNAPSHOTS):
        raise ValueError("Manifest must contain the three expected source snapshots")
    return manifest


def verify(archive: Path, source: Path | None = None) -> dict:
    """Verify stored bytes and, optionally, the original/reconstructed trees."""
    manifest = load_manifest(archive)
    checked_storage = {}
    referenced = 0
    for snapshot, record in manifest["snapshots"].items():
        expected = {}
        for relative, entry in record["files"].items():
            safe_child(archive, relative)
            stored = safe_child(archive, entry["stored_path"])
            stored_relative = stored.relative_to(archive).as_posix()
            allowed = stored_relative == f"{CANONICAL}/{relative}" or (
                snapshot != CANONICAL
                and stored_relative == f"source_overlays/{snapshot}/{relative}"
            )
            if not allowed:
                raise ValueError(f"Unexpected storage mapping: {snapshot}/{relative}")
            if stored_relative not in checked_storage:
                checked_storage[stored_relative] = {
                    "bytes": stored.stat().st_size,
                    "sha256": digest(stored),
                }
            expected[relative] = {"bytes": entry["bytes"], "sha256": entry["sha256"]}
            if checked_storage[stored_relative] != expected[relative]:
                raise ValueError(f"Stored file failed hash/size verification: {stored_relative}")
            referenced += 1
        for relative in record["directories"]:
            safe_child(archive, relative)
        if source is not None:
            actual, directories = inventory(source / snapshot)
            if actual != expected or directories != record["directories"]:
                missing = sorted(set(expected) - set(actual))
                extra = sorted(set(actual) - set(expected))
                changed = sorted(k for k in set(expected) & set(actual) if expected[k] != actual[k])
                raise ValueError(
                    f"Snapshot mismatch in {snapshot}: missing={missing}, extra={extra}, "
                    f"changed={changed}, directories_match={directories == record['directories']}"
                )
    return {
        "verified_original_files": referenced,
        "stored_files": len(checked_storage),
        "deduplicated_file_copies": referenced - len(checked_storage),
        "original_bytes": sum(
            f["bytes"] for s in manifest["snapshots"].values() for f in s["files"].values()
        ),
        "stored_bytes": sum(f["bytes"] for f in checked_storage.values()),
        "source_comparison": source is not None,
    }


def reconstruct(archive: Path, output: Path) -> dict:
    """Materialize and byte-verify all original snapshots without overwriting."""
    verify(archive)
    manifest = load_manifest(archive)
    prepare_output(output, archive)
    for snapshot, record in manifest["snapshots"].items():
        target = output / snapshot
        target.mkdir()
        for relative in record["directories"]:
            safe_child(target, relative).mkdir(parents=True, exist_ok=True)
        for relative, entry in record["files"].items():
            destination = safe_child(target, relative)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(safe_child(archive, entry["stored_path"]), destination)
    return verify(archive, output)


def compact(source: Path, output: Path) -> dict:
    inventories = {snapshot: inventory(source / snapshot) for snapshot in SNAPSHOTS}
    prepare_output(output, source)
    canonical_files = inventories[CANONICAL][0]
    manifest = {
        "schema_version": 1,
        "canonical_snapshot": CANONICAL,
        "content_policy": "Original file bytes retained; no source/path/license rewriting",
        "restoration": "Use compact_pilot_sources.py restore --archive ARCHIVE --output EMPTY_DIRECTORY",
        "snapshots": {},
    }
    for snapshot in SNAPSHOTS:
        files, directories = inventories[snapshot]
        entries = {}
        for relative, facts in files.items():
            if snapshot == CANONICAL or canonical_files.get(relative) == facts:
                stored_relative = f"{CANONICAL}/{relative}"
            else:
                stored_relative = f"source_overlays/{snapshot}/{relative}"
            destination = safe_child(output, stored_relative)
            if not destination.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source / snapshot / relative, destination)
            entries[relative] = {"stored_path": stored_relative, **facts}
        manifest["snapshots"][snapshot] = {"directories": directories, "files": entries}
    # Preserve empty directories in the readable canonical snapshot too.
    (output / CANONICAL).mkdir(exist_ok=True)
    for relative in inventories[CANONICAL][1]:
        safe_child(output / CANONICAL, relative).mkdir(parents=True, exist_ok=True)
    (output / MANIFEST).write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    result = verify(output, source)
    # Reconstruct into an isolated temporary sibling and verify all three trees.
    # TemporaryDirectory removes only its own newly created temporary directory.
    with tempfile.TemporaryDirectory(prefix="pilot-source-verify-", dir=output.parent) as temporary:
        reconstruct(output, Path(temporary) / "restored")
    result["byte_exact_reconstruction_verified"] = True
    result["snapshot_file_counts"] = {name: len(files) for name, (files, _) in inventories.items()}
    result["overlay_file_counts"] = {
        name: sum(f["stored_path"].startswith("source_overlays/") for f in s["files"].values())
        for name, s in manifest["snapshots"].items() if name != CANONICAL
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    actions = parser.add_subparsers(dest="action", required=True)
    command = actions.add_parser("compact", help="Create a canonical tree and exact version overlays")
    command.add_argument("--source", type=Path, required=True, help="Existing method_sources directory")
    command.add_argument("--output", type=Path, required=True, help="Absent or empty staging directory")
    command = actions.add_parser("verify", help="Verify hashes and optionally compare with complete snapshots")
    command.add_argument("--archive", type=Path, required=True)
    command.add_argument("--source", type=Path)
    command = actions.add_parser("restore", help="Reconstruct and verify all three historical snapshots")
    command.add_argument("--archive", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == "compact":
            result = compact(args.source.resolve(), args.output.resolve())
        elif args.action == "verify":
            result = verify(args.archive.resolve(), args.source.resolve() if args.source else None)
        else:
            result = reconstruct(args.archive.resolve(), args.output.resolve())
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"action": args.action, "status": "verified", **result}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
