#!/usr/bin/env python3
"""Bridge validated R2 assemblies to preserved A0-A4 and audited annotation crosswalks."""
import argparse
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import contextlib
import spark_runtime
import sys

from reads_io import dump, ident, sha, table, write_tsv

R=Path(__file__).resolve().parents[1]
HEAVY_LOCK=Path.home()/"scratch/crc-phage/core-heavy.lock"
GENE_FIELDS={"gene_id","genome_id","contig_id","start0","end0","strand"}


def read(path):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([path])
    return json.loads(Path(path).read_text())


def resolve(value,base):
    p=Path(value).expanduser()
    p=(p if p.is_absolute() else base/p).resolve()
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([p])
    return p


def guard(out,bridge,protected):
    for target in (out,bridge):
        target=Path(target).resolve()
        for source in protected:
            source=Path(source).resolve()
            if spark_runtime.enabled():
                from spark_fs_guard import reject_ipc_storage
                reject_ipc_storage([source])
            if target==source or target.is_relative_to(source) or source.is_relative_to(target):
                raise ValueError(f"Bridge/core output overlaps protected input: {target} / {source}")


def declared_sources(assembly):
    """Collect paths before writing a lock or receipt into potentially conflicting output."""
    paths=[]
    try:
        units=read(assembly/"assembly_manifest.json")["units"]
        for u in units:
            for key in ("assembly_fasta","genbank","contig_map"):
                if u.get(key):paths.append(resolve(u[key],assembly))
        for g in table(assembly/"genomes.tsv"):
            for key in ("fasta","genbank"):
                if g.get(key):paths.append(resolve(g[key],assembly))
    except (KeyError,ValueError,OSError,TypeError):
        # Actual validation below records a failed attempt; never invent missing units.
        pass
    return paths


def file_digest(path,label):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([path])
    if not path.is_file():raise ValueError(f"Missing {label}: {path}")
    return sha(path)


def validate_inputs(assembly,metadata,reference):
    if spark_runtime.enabled():
        from spark_fs_guard import reject_ipc_storage
        reject_ipc_storage([assembly/'assembly_manifest.json',assembly/'genomes.tsv',metadata])
    manifest=read(assembly/"assembly_manifest.json")
    if manifest.get("schema_version")!=1 or not isinstance(manifest.get("units"),list):
        raise ValueError("Unsupported R2 assembly manifest")
    units=manifest["units"];genomes=table(assembly/"genomes.tsv")
    if not genomes:raise ValueError("No successful bacterial-isolate assemblies; bulk/VLP use V1")
    if len({u["unit_id"] for u in units})!=len(units) or len({g["genome_id"] for g in genomes})!=len(genomes):
        raise ValueError("Duplicate assembly/genome IDs")
    byid={u["unit_id"]:u for u in units}
    expected={u["unit_id"] for u in units if u.get("material_type")=="bacterial_isolate" and u.get("status")=="success"}
    if {g["genome_id"] for g in genomes}!=expected:
        raise ValueError("genomes.tsv must contain every and only successful bacterial-isolate unit")
    banks=[];normalized=[];files={}
    for g in sorted(genomes,key=lambda row:row["genome_id"]):
        gid=ident(g["genome_id"],"genome_id");u=byid[gid]
        fasta=resolve(g["fasta"],assembly)
        if resolve(u["assembly_fasta"],assembly)!=fasta:raise ValueError("Manifest/genomes FASTA paths disagree")
        expected_sha=u.get("sha256","")
        if not re.fullmatch("[a-f0-9]{64}",expected_sha):raise ValueError("Missing/invalid assembly SHA256")
        if file_digest(fasta,"assembly")!=expected_sha:raise ValueError("Assembly checksum changed")
        files[str(fasta)]=expected_sha
        bank=resolve(g["genbank"],assembly) if g.get("genbank") else None
        unit_bank=resolve(u["genbank"],assembly) if u.get("genbank") else None
        if unit_bank is not None and unit_bank!=bank:raise ValueError("Manifest/genomes GBFF paths disagree")
        digest=""
        if bank is not None:
            expected_sha=u.get("genbank_sha256","")
            if not re.fullmatch("[a-f0-9]{64}",expected_sha):raise ValueError("Provided GBFF requires a recorded SHA256")
            digest=file_digest(bank,"GBFF")
            if digest!=expected_sha or (g.get("genbank_sha256") and g["genbank_sha256"]!=digest):
                raise ValueError("GBFF checksum changed")
            files[str(bank)]=digest
        banks.append(dict(genome_id=gid,genbank=str(bank) if bank else "",genbank_sha256=digest))
        normalized.append(dict(genome_id=gid,fasta=str(fasta)))
    meta=table(metadata)
    if not meta or any(not row.get("genome_id") for row in meta):
        raise ValueError("Metadata must declare genome_id rows")
    if len({r["genome_id"] for r in meta})!=len(meta):raise ValueError("Duplicate metadata genome IDs")
    if not expected.issubset({r["genome_id"] for r in meta}):raise ValueError("Metadata missing bridge genome IDs")
    files[str(metadata)]=file_digest(metadata,"metadata")
    for name in ("assembly_manifest.json","genomes.tsv"):
        files[str(assembly/name)]=file_digest(assembly/name,"R2 contract")
    if reference is not None:
        files[str(reference/"reference_manifest.json")]=file_digest(reference/"reference_manifest.json","reference manifest")
    binding=dict(assembly_result=str(assembly),metadata=str(metadata),reference_dir=str(reference) if reference else None,
                 input_files=files,genome_ids=sorted(expected))
    return normalized,banks,binding


def run_logged(command,log,cwd=None):
    with log.open("w") as handle:
        handle.write(json.dumps(list(map(str,command)))+"\n");handle.flush()
        subprocess.run(list(map(str,command)),check=True,cwd=cwd,stdout=handle,stderr=subprocess.STDOUT)


def core_receipt(out,previous_attempt):
    record=read(out/"pipeline_info/run_status.json")
    workflow=read(out/"pipeline_info/workflow_status.json")
    final=read(out/"core_data_manifest.json")
    if record.get("status")!="completed" or workflow.get("success") is not True or final.get("status")!="completed":
        raise ValueError("Core exited without completed run/workflow/aggregate evidence")
    if not record.get("attempt") or record["attempt"]==previous_attempt:
        raise ValueError("Core did not write a fresh attempt receipt")
    session=record.get("nextflow_session_id")
    if not session or session!=workflow.get("session_id"):
        raise ValueError("Core workflow session does not match its current run receipt")
    expected=final.get("outputs",{}).get("gene_table.tsv")
    if not expected or file_digest(out/"gene_table.tsv","canonical genes")!=expected:
        raise ValueError("Core canonical gene table missing or checksum mismatch")
    genes=table(out/"gene_table.tsv")
    if any(not GENE_FIELDS.issubset(row) for row in genes):raise ValueError("Canonical gene schema incomplete")
    if len({row["gene_id"] for row in genes})!=len(genes):raise ValueError("Duplicate canonical gene IDs")
    return genes,session


def crosswalks(out,banks,genes,image,attempt,resume):
    cross=out/"annotation_crosswalk";cross.mkdir(exist_ok=True)
    statuses=[];all_ids={b["genome_id"] for b in banks}
    if any(g["genome_id"] not in all_ids for g in genes):raise ValueError("Canonical gene genome foreign key mismatch")
    for bank in banks:
        gid=bank["genome_id"];row=dict(genome_id=gid,status="not_assessed",reason="no_GBFF",execution="not_run",crosswalk_sha256="")
        if bank["genbank"]:
            try:
                subset=[g for g in genes if g["genome_id"]==gid]
                if not subset:raise ValueError("No canonical genes for GBFF crosswalk")
                gene_file=cross/f"{gid}.canonical.tsv"
                write_tsv(gene_file,subset,list(subset[0]))
                result=cross/f"{gid}.crosswalk.tsv";receipt=cross/f"{gid}.receipt.json"
                inputs=dict(genes_sha256=sha(gene_file),genbank_sha256=bank["genbank_sha256"],container=image,
                            adapter_sha256=sha(R/"bin/reads_assembly.py"),reads_io_sha256=sha(R/"bin/reads_io.py"))
                if spark_runtime.enabled():
                    inputs.update(runtime="native-locked-prefix",native_lock_sha256=spark_runtime.environment("reads-bakta")[0]["lock_sha256"],guard_sha256=sha(R/"bin/spark_fs_guard.py"))
                prior=read(receipt) if receipt.exists() else {}
                cached=bool(resume and result.is_file() and prior.get("inputs")==inputs
                            and prior.get("status")=="success" and prior.get("output_sha256")==sha(result))
                if not cached:
                    if result.exists():
                        (attempt/f"{gid}.previous_crosswalk.tsv").write_bytes(result.read_bytes())
                        result.unlink()
                    mounts={str(R):"ro",str(out):"rw",str(Path(bank["genbank"]).parent):"ro"}
                    command=["docker","run","--rm","--platform","linux/amd64","--cpus","1","--memory","2g",
                             "--user",f"{os.getuid()}:{os.getgid()}","--workdir",str(cross)]
                    for path,mode in mounts.items():command+=["-v",f"{path}:{path}:{mode}"]
                    command += [image,"python3",str(R/"bin/reads_assembly.py"),"crosswalk","--genes",str(gene_file),
                                "--gbff",bank["genbank"],"--out",str(result)]
                    if spark_runtime.enabled():
                        spark_runtime.validate_output_inputs(out,[bank["genbank"]])
                        command=spark_runtime.command("reads-bakta",["python3",str(R/"bin/reads_assembly.py"),"crosswalk","--genes",str(gene_file),"--gbff",bank["genbank"],"--out",str(result)],writable=[cross])
                    run_logged(command,attempt/f"{gid}.crosswalk.log")
                    if not result.is_file() or not result.stat().st_size or "canonical_gene_id" not in result.open().readline().rstrip().split("\t"):
                        raise ValueError("Crosswalk command produced no valid output table")
                    dump(receipt,dict(status="success",inputs=inputs,output_sha256=sha(result)))
                (cross/f"{gid}.failure.json").unlink(missing_ok=True)
                row.update(status="assessed",reason="overlap_coordinates_no_forced_one_to_one",
                           execution="cached" if cached else "completed",crosswalk_sha256=sha(result))
            except Exception as error:
                row.update(status="failed",reason=str(error),execution="failed")
                dump(cross/f"{gid}.failure.json",row)
        statuses.append(row)
        write_tsv(cross/"status.tsv",statuses,["genome_id","status","reason","execution","crosswalk_sha256"])
    return statuses


def execute(args):
    assembly=Path(args.assembly_result).expanduser().resolve();out=Path(args.outdir).expanduser().resolve()
    metadata=Path(args.metadata).expanduser().resolve()
    reference=Path(args.reference_dir).expanduser().resolve() if args.reference_dir else None
    bridge=(out.parent/(out.name+"_bridge")).resolve()
    protected=[assembly,metadata,R/"bin",R/"assets"]+([reference] if reference else [])+declared_sources(assembly)
    try:guard(out,bridge,protected)
    except Exception as error:
        print(json.dumps(dict(status="failed",stage="unsafe_output_rejected_before_write",error=str(error))),file=sys.stderr)
        return 1
    bridge.mkdir(parents=True,exist_ok=True)
    with (bridge/"bridge.lock").open("a") as output_lock:
        try:fcntl.flock(output_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print("Bridge output is already active; existing receipt preserved",file=sys.stderr);return 1
        state=bridge/"bridge_manifest.json";previous=read(state) if state.exists() else {}
        stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        attempt=bridge/"attempts"/stamp;attempt.mkdir(parents=True)
        if previous:dump(attempt/"previous_bridge_manifest.json",previous)
        record=dict(status="running",stage="validation",attempt=stamp,outdir=str(out),actual_repo_path=str(R),
                    binding=previous.get("binding"),resource_limits=dict(cpus=6,memory_gib=12),
                    code_sha256=sha(Path(__file__)),real_data_validated=False,scientifically_calibrated=False)
        dump(state,record);core_verified=False
        try:
            if previous and not args.resume:raise ValueError("Bridge already exists; use --resume or a new output")
            if args.resume and not previous:raise ValueError("No bridge receipt exists for --resume")
            if out.exists() and any(out.iterdir()) and not args.resume:raise ValueError("Core output is not empty")
            normalized,banks,binding=validate_inputs(assembly,metadata,reference)
            if args.resume and previous.get("binding") and previous["binding"]!=binding:
                raise ValueError("Bridge resume inputs changed; use a new output")
            record["binding"]=binding
            genome_file=bridge/"genomes.tsv";bank_file=bridge/"genbanks.tsv"
            write_tsv(genome_file,normalized,["genome_id","fasta"])
            write_tsv(bank_file,banks,["genome_id","genbank","genbank_sha256"])
            images=read(R/"assets/reads_tools.json")
            image=images.get("reads-bakta",{}).get("container","")
            if any(b["genbank"] for b in banks) and not (re.fullmatch(r"sha256:[a-f0-9]{64}",image) or re.search(r"@sha256:[a-f0-9]{64}$",image)):
                raise ValueError("Crosswalk requires a digest-pinned reads-bakta image")
            core_state=out/"pipeline_info/run_status.json"
            if args.resume and out.exists() and any(out.iterdir()) and not core_state.exists():
                raise ValueError("Partial core output has no saved run status; cannot safely resume")
            before=read(core_state).get("attempt") if core_state.exists() else None
            command=[str(R/"bin/run-core"),"--input",str(genome_file),"--outdir",str(out),
                     "--genbanks",str(bank_file),"--metadata",str(metadata)]
            if reference:command+=["--reference-dir",str(reference)]
            if args.resume and core_state.exists():command+=["--resume"]
            record.update(command=command,status="waiting_resources",stage="core")
            dump(state,record)
            HEAVY_LOCK.parent.mkdir(parents=True,exist_ok=True)
            # run_core.py holds only its own output lock; it does not acquire this
            # shared heavy lock. One parent lock covers core plus crosswalk work.
            with (contextlib.nullcontext() if spark_runtime.enabled() else HEAVY_LOCK.open("a")) as heavy:
                if heavy is not None:fcntl.flock(heavy,fcntl.LOCK_EX)
                else:spark_runtime.allocation()
                record["status"]="running";dump(state,record)
                run_logged(command,attempt/"core_console.log",R)
                genes,session=core_receipt(out,before);core_verified=True
                record.update(stage="crosswalk",core_status="completed",core_session_id=session,
                              core_manifest_sha256=sha(out/"core_data_manifest.json"))
                dump(state,record)
                statuses=crosswalks(out,banks,genes,image,attempt,args.resume)
            failed=[s for s in statuses if s["status"]=="failed"]
            record.update(status="partial_failure" if failed else "success",stage="complete",
                          annotation_crosswalk=str(out/"annotation_crosswalk"),crosswalk_statuses=statuses,
                          exit_code=1 if failed else 0)
            if failed:record["error"]="One or more annotation crosswalks failed; core outputs preserved"
        except Exception as error:
            record.update(status="partial_failure" if core_verified else "failed",error=str(error),exit_code=1)
        record["completed_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
        dump(state,record);dump(attempt/"bridge_manifest.json",record)
        print(json.dumps(dict(status=record["status"],outdir=str(out),receipt=str(state),error=record.get("error"))))
        return record["exit_code"]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assembly-result",required=True);p.add_argument("--outdir",required=True)
    p.add_argument("--resume",action="store_true");p.add_argument("--metadata",required=True)
    p.add_argument("--reference-dir")
    return execute(p.parse_args())


if __name__=="__main__":raise SystemExit(main())
