#!/usr/bin/env python3
"""Validate study metadata and write a fresh audit report. No run/resume/submit action."""
import argparse
import csv
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from core_preflight import fresh_result_dir, check_platform
from study_contract import load_contract, evaluate, local_path


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def check(manifest,outdir):
    check_platform();out=fresh_result_dir(outdir)
    report=dict(mode='metadata_check_only',schema_valid=False,contract_valid=False,execution_authorized=False,nextflow_started=False,tasks_submitted=0,tools_executed=0,network_requests=0,real_data_validated=False,stage_readiness='not_ready')
    tables=None;cfg=None
    try:
        cfg,tables,evidence,protected=load_contract(manifest)
        for p in protected:
            if p==out or p.is_relative_to(out) or out.is_relative_to(p):raise ValueError('Output overlaps declared input/identity path: '+str(p))
        report=evaluate(cfg,tables);report['input_manifest']=str(Path(manifest).resolve());report['input_evidence']=evidence
    except (ValueError,OSError,KeyError,TypeError) as exc:
        report['errors']=[dict(code='invalid_contract',message=str(exc))]
    # Only new direct results child, never an input path; load/overlap errors do not write.
    if not report.get('schema_valid'):
        return report,2
    out.mkdir()  # atomic nonexistence check; no exist_ok and no overwrite
    write_json(out/'report.json',report);write_json(out/'parameter_plan.json',report['parameter_plan'])
    if report['contract_valid']:
        write_json(out/'study_tables.preserved.json',tables)
        with (out/'feature_coordinates.draft.tsv').open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['feature_id','genome_id','contig_id','start0','end0','strand','source_coordinate_system'])
            for r in tables['features']:
                if r['start'] is not None:w.writerow([r['feature_id'],r['genome_id'],r['contig_id'],r['start']-(1 if r['coordinate_system']=='one_based_closed' else 0),r['end'],r['strand'],r['coordinate_system']])
        # Existing core supports only genome_id/fasta and assembly metadata. The full
        # study source tables and every hash remain in input_evidence; no clinical
        # information is represented as consumed by the old core workflow.
        selected=[g for g in tables['genomes'] if g['master_member'] and g['role']=='discovery' and g['fasta_path']]
        with (out/'core_input.draft.tsv').open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['genome_id','fasta'])
            for g in selected:w.writerow([g['genome_id'],str(local_path(g['fasta_path'],Path(manifest).resolve().parent))])
        with (out/'core_metadata.draft.tsv').open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['genome_id','ncbi_assembly_accession'])
            for g in selected:w.writerow([g['genome_id'],g['assembly_accession'] or ''])
        write_json(out/'adaptation_manifest.json',dict(status='draft_not_executable_research_plan',source_tables=report['input_evidence'],preserved_columns='all original fields retained in hashed source tables; only two limited core adapters emitted',rows_adapted=len(selected),master_rows_without_path=[g['genome_id'] for g in tables['genomes'] if g['master_member'] and not g['fasta_path']],unsupported='This preflight emits draft assembly adapters only. Dedicated host_source/community_analysis consumers require their own frozen model contracts and runtime evidence; no run command generated',parameter_overrides_activated=False))
    (out/'README.md').write_text('# Study metadata check only\n\nExit 0 means a structurally valid draft with consistent declared relationships, not readiness to analyze.\n\nStage: '+report['stage_readiness']+'; actual biological inputs and evidence truth were not assessed.\nSee report.json and parameter_plan.json. *.draft.tsv are limited field adapters, not approved inputs. No analysis command was generated.\n')
    return report,0 if report['contract_valid'] else 2


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args(argv)
    try:report,code=check(a.manifest,a.outdir)
    except (ValueError,OSError) as exc:print(str(exc),file=sys.stderr);return 2
    print(json.dumps(dict(contract_valid=report['contract_valid'],stage_readiness=report['stage_readiness'],execution_authorized=False,errors=report.get('errors',[])),ensure_ascii=False))
    return code

if __name__=='__main__':raise SystemExit(main())
