#!/usr/bin/env python3
"""Auditable adapters to B2 long.tsv; absent/unknown measurements never become zero."""
import argparse
import json
from pathlib import Path
from association import load_spec, read_tsv, sha, write_tsv
from viral_catalog import POLICY, eligibility

MISSING={"","na","nan","null","none","unknown","not_assessed"}
TRUE={"true","eligible","pass","1","yes"}

def known(value):
    return str(value).strip().lower() not in MISSING

def unique(rows,key,label):
    result={row[key]:row for row in rows}
    if len(result)!=len(rows):raise ValueError("Duplicate "+key+" in "+label)
    return result

def candidate_decision(hit,mapping,host_qc,scope):
    """Legacy core is explicitly the bacterial-genome route; never infer host from taxonomy."""
    reasons=[]
    if scope=="host_resolved_prophage":
        if hit["candidate_type"]!="provirus_locus":
            return "excluded","outside_provirus_locus_estimand"
        if not known(hit.get("source_host_genome")) or hit["source_host_genome"]!=hit["genome_id"]:
            reasons.append("bacterial_source_host_missing_or_mismatched")
        if not known(hit.get("integration_evidence")):
            reasons.append("integration_evidence_unknown")
    if hit.get("caller")!="geNomad":reasons.append("viral_identity_not_supported_by_reviewed_caller")
    if hit["candidate_type"] not in {"provirus_locus","full_contig_virus"}:
        reasons.append("boundary_extent_unknown")
    if not mapping or str(mapping.get("slice_verified")).lower()!="true":
        reasons.append("upstream_source_slice_not_verified")
    elif mapping.get("coordinate_system")!="0-based-half-open":
        reasons.append("coordinate_system_not_supported")
    else:
        for key in ("candidate_id","genome_id","contig_id","start0","end0","strand","length","sequence_sha256","boundary_version"):
            if str(mapping.get(key,""))!=str(hit.get(key,"")):
                raise ValueError("Master/coordinate mapping mismatch: "+hit["candidate_id"]+"/"+key)
        start,end=int(hit["start0"]),int(hit["end0"])
        if not 0<=start<end or end-start!=int(hit["length"]) or hit["strand"] not in {"+","-"}:
            raise ValueError("Invalid verified candidate coordinates")
    if not host_qc or host_qc.get("checkm2_status")!="completed" or host_qc.get("eligibility")!="eligible":
        reasons.append("source_genome_QC_not_eligible")
    if str(hit.get("main_catalog_eligible","unknown")).lower() not in TRUE:
        reasons.append("legacy_main_catalog_eligibility_missing_or_false")
    if hit.get("checkv_status")!="completed":reasons.append("CheckV_not_completed")
    if reasons:return "unknown",";".join(reasons)
    quality,reason=eligibility(dict(hit,virus_scope=hit.get("virus_scope","unknown")))
    return ("excluded" if quality=="ineligible" else quality),reason

def genome(core,spec,output):
    receipt=json.loads((core/"pipeline_info/workflow_status.json").read_text())
    if receipt.get("success") is not True:raise ValueError("Core workflow must have a successful receipt")
    catalog=json.loads((core/"catalog_manifest.json").read_text())
    if catalog.get("catalog_version")!=spec["catalog_version"]:raise ValueError("Declared catalog version differs from actual core catalog")
    signed=core/"core_data_manifest.json"
    if signed.exists():
        hashes=json.loads(signed.read_text()).get("outputs",{})
        for name in ("sample_status_core.tsv","genome_qc.tsv","vOTU_members.tsv","prophage_master.tsv","coordinate_mapping.tsv","taxonomy.tsv","catalog_manifest.json"):
            if not hashes.get(name) or sha(core/name)!=hashes[name]:raise ValueError("Core signed output missing/changed: "+name)
    if spec["catalog_scope"]=="virome_discovery":raise ValueError("Bacterial-genome adapter cannot create a virome-discovery estimand")
    qc,_=read_tsv(core/"genome_qc.tsv",{"genome_id","eligibility","checkm2_status"})
    statuses,_=read_tsv(core/"sample_status_core.tsv",{"genome_id","A0_A1_status","n_candidates"})
    bystatus=unique(statuses,"genome_id","sample status")
    if set(bystatus)!={r["genome_id"] for r in qc}:raise ValueError("Sample/QC denominator genomes differ")
    members,_=read_tsv(core/"vOTU_members.tsv",{"genome_id","candidate_id","votu_id","membership_status"})
    master,_=read_tsv(core/"prophage_master.tsv",{"candidate_id","genome_id","candidate_type"})
    coords,_=read_tsv(core/"coordinate_mapping.tsv",{"candidate_id","slice_verified","coordinate_system"})
    taxonomy,_=read_tsv(core/"taxonomy.tsv",{"genome_id","taxonomy_status","reference_supported_species"})
    byqc=unique(qc,"genome_id","QC");bytax=unique(taxonomy,"genome_id","taxonomy")
    bycandidate=unique(master,"candidate_id","master");bycoord=unique(coords,"candidate_id","coordinates")
    if any(row["genome_id"] not in byqc for row in master):raise ValueError("Candidate genome absent from QC")
    if any(row["candidate_id"] not in bycandidate or row["genome_id"]!=bycandidate[row["candidate_id"]]["genome_id"] for row in members):
        raise ValueError("vOTU member foreign key mismatch")
    unique(members,"candidate_id","vOTU membership")
    known_votu={row["votu_id"] for row in members if known(row["votu_id"])}
    allowed=known_votu|({"total_candidate_count"} if spec["outcome"]=="count" else set())
    if not set(spec["feature_ids"]).issubset(allowed):raise ValueError("Feature absent from catalog; unmeasured feature cannot become zero")
    decisions={};audit=[]
    for hit in master:
        decision,reason=candidate_decision(hit,bycoord.get(hit["candidate_id"]),byqc.get(hit["genome_id"]),spec["catalog_scope"])
        decisions[hit["candidate_id"]]=(decision,reason)
        audit.append(dict(candidate_id=hit["candidate_id"],genome_id=hit["genome_id"],candidate_type=hit["candidate_type"],
            source_kind="legacy_bacterial_genome",catalog_scope=spec["catalog_scope"],decision=decision,reason=reason,
            slice_evidence="upstream_verified_coordinates_and_hash_crosscheck" if bycoord.get(hit["candidate_id"],{}).get("slice_verified")=="true" else "unknown",
            source_slice_recomputed_this_adapter=False,viral_eligibility_policy=POLICY["version"]))
    eligible_features={r["votu_id"] for r in members if decisions[r["candidate_id"]][0]=="eligible" and r["membership_status"] in {"representative","member"}}
    rows=[];unit_audit=[]
    for gid,item in sorted(byqc.items()):
        tx=bytax.get(gid,{})
        unit_reasons=[]
        sample=bystatus[gid]
        if sample["A0_A1_status"]!="completed":unit_reasons.append("caller_not_completed:"+sample["A0_A1_status"])
        elif int(sample["n_candidates"])!=sum(r["genome_id"]==gid for r in master):raise ValueError("Caller sample/master candidate counts differ")
        if item["eligibility"]!="eligible" or item["checkm2_status"]!="completed":
            unit_reasons.append("host_QC_not_eligible")
        # Species is analysis eligibility for this fixed-species template, not global catalog membership.
        if tx.get("taxonomy_status")!="reference_supported" or not known(tx.get("reference_supported_species")):
            unit_reasons.append("analysis_species_unresolved")
        elif tx["reference_supported_species"]!=spec["target_species"]:
            unit_reasons.append("outside_target_species")
        gm=[r for r in members if r["genome_id"]==gid]
        assigned={r["candidate_id"] for r in gm if r["membership_status"] in {"representative","member"}}
        unassigned=[h["candidate_id"] for h in master if h["genome_id"]==gid and decisions[h["candidate_id"]][0]=="eligible" and h["candidate_id"] not in assigned]
        if unassigned:unit_reasons.append("eligible_candidate_without_resolved_vOTU:"+",".join(sorted(unassigned)))
        if any(r["membership_status"]=="ambiguous_multiple_representatives" for r in gm):
            unit_reasons.append("ambiguous_vOTU_membership")
        unit_audit.append(dict(unit_id=gid,analysis_status="eligible" if not unit_reasons else "unknown",
            reason=";".join(unit_reasons),reference_supported_species=tx.get("reference_supported_species","unknown"),
            target_species=spec["target_species"],catalog_membership_independent_of_target_species=True))
        for fid in spec["feature_ids"]:
            reasons=list(unit_reasons)
            applicable=[r for r in gm if fid=="total_candidate_count" or r["votu_id"]==fid]
            unknown=[r for r in applicable if decisions[r["candidate_id"]][0]=="unknown"]
            if unknown:reasons.append("candidate_qualification_unknown:"+",".join(sorted(r["candidate_id"] for r in unknown)))
            if fid!="total_candidate_count" and fid not in eligible_features:reasons.append("feature_has_no_eligible_catalog_member")
            hits={r["candidate_id"] for r in applicable if r["membership_status"] in {"representative","member"} and decisions[r["candidate_id"]][0]=="eligible"}
            status="unknown" if reasons else "assessed"
            value=len(hits) if spec["outcome"]=="count" else int(bool(hits))
            rows.append(dict(unit_id=gid,feature_id=fid,value=value if status=="assessed" else "",status=status,
                reason=";".join(reasons) or "assessed_development_catalog_policy"))
    write_tsv(output,["unit_id","feature_id","value","status","reason"],rows)
    write_tsv(output.with_suffix(output.suffix+".candidates.tsv"),
        ["candidate_id","genome_id","candidate_type","source_kind","catalog_scope","decision","reason","slice_evidence","source_slice_recomputed_this_adapter","viral_eligibility_policy"],audit)
    write_tsv(output.with_suffix(output.suffix+".units.tsv"),
        ["unit_id","analysis_status","reason","reference_supported_species","target_species","catalog_membership_independent_of_target_species"],unit_audit)
    return [core/p for p in ("sample_status_core.tsv","genome_qc.tsv","vOTU_members.tsv","prophage_master.tsv","coordinate_mapping.tsv","taxonomy.tsv","catalog_manifest.json","pipeline_info/workflow_status.json")]+([signed] if signed.exists() else [])

def b1(source,spec,output):
    rows,fields=read_tsv(source,{"unit_id","feature_id","value","status"})
    if not set(r["feature_id"] for r in rows).issubset(spec["feature_ids"]):raise ValueError("B1 contains features outside declared hypothesis family")
    write_tsv(output,["unit_id","feature_id","value","status"]+(["reason"] if "reason" in fields else []),rows)
    return [source]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kind",choices=["genome","b1"],required=True)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--spec",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():p.error("Refusing to overwrite adapter output")
    spec=load_spec(args.spec)
    if args.kind=="genome" and spec["input_kind"]!="genome":p.error("Genome adapter requires genome spec")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    sources=genome(args.source,spec,args.output) if args.kind=="genome" else b1(args.source,spec,args.output)
    outputs=[args.output]+list(args.output.parent.glob(args.output.name+".*.tsv"))
    manifest=dict(adapter=args.kind,scientific_status="development_exploratory",
        input_hashes={str(f.resolve()):sha(f) for f in sources+[args.spec]},
        output_sha256=sha(args.output),output_hashes={f.name:sha(f) for f in outputs},
        catalog_version=spec["catalog_version"],catalog_scope=spec["catalog_scope"],
        absence_semantics="not_detected_by_assay_in_assessable_unit",
        viral_eligibility_policy=POLICY if args.kind=="genome" else None,
        adapter_sha256=sha(Path(__file__)),shared_eligibility_code_sha256=sha(Path(__file__).with_name("viral_catalog.py")))
    args.output.with_suffix(args.output.suffix+".manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")

if __name__=="__main__":main()
