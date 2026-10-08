"""Quantification-group report; no expansion of one group's measurements to vOTUs."""
from pathlib import Path
import html,shutil
from cohort_common import read_json,write_json,write_table,sha256,require,table
from cohort_measure import MEASUREMENT_FIELDS,DETECTION_FIELDS,HOST_FIELDS
def aggregate_groups(root,out,units,measurements,detections,hosts,manifests):
    ref=read_json(root/'reference/reference_manifest.json')
    mapping={r['reference_id']:r for r in ref['references'] if r['reference_kind']=='target'}
    for row in measurements+detections:
        target=mapping[row['reference_id']]
        require(row.get('quantification_group')==target['quantification_group'] and row.get('votu_id') in ('','NA'),'Group identity mislabeled as vOTU')
        for k in ('biological_votu_ids','reference_host_id','host_assignment_status'):
            require(row.get(k)==str(target.get(k,'')),'Group provenance mismatch: '+k)
    write_table(out/'measurement_long.tsv',MEASUREMENT_FIELDS,measurements)
    write_table(out/'detection_long.tsv',DETECTION_FIELDS,detections)
    write_table(out/'host_core.tsv',HOST_FIELDS,hosts)
    write_table(out/'host_backbone.tsv',HOST_FIELDS,hosts)
    metadata=[{**u.get('metadata',{}),**{k:u.get(k) for k in ['unit_id','biological_sample_id','subject_id','library_id','material_type','layout']}} for u in units.values()]
    fields=sorted(set().union(*(set(r) for r in metadata))) if metadata else ['unit_id','biological_sample_id','subject_id']
    write_table(out/'unit_metadata.tsv',fields,metadata)
    shutil.copyfile(root/'reference/quantification_group_members.tsv',out/'quantification_group_members.tsv')
    hostindex={(r['unit_id'],r['host_id']):r for r in hosts};mi={(r['unit_id'],r['reference_id']):r for r in measurements};signal=[]
    for d in detections:
        m=mi[(d['unit_id'],d['reference_id'])];h=hostindex.get((d['unit_id'],d['reference_host_id']),{})
        reason=''
        if d['host_assignment_status']!='single_reference_host_not_sample_phasing':reason='multiple_or_unknown_reference_hosts_no_sample_link'
        elif not h or h.get('status')!='assessed' or h.get('host_adjustment_eligibility') not in ('eligible','provisional'):reason='host_core_unassessable'
        elif float(h.get('mean_backbone_depth') or 0)<=0:reason='host_core_not_detected'
        elif d['inference_status']!='assessed':reason=d['inference_reason']
        if str(m['sample_in_discovery_set']).lower() in ('true','unknown'):reason='discovery_overlap_or_unknown'
        signal.append(dict(unit_id=d['unit_id'],biological_sample_id=d['biological_sample_id'],subject_id=d['subject_id'],
            reference_host_id=d['reference_host_id'],quantification_group=d['quantification_group'],detection_profile_id=d['detection_profile_id'],
            status='assessed' if not reason else 'unassessable',reason=reason,value=d['value'] if not reason else '',
            fragment_rpkm_all_clean=m['fragment_rpkm_all_clean'],host_core_depth=h.get('mean_backbone_depth',''),
            host_core_breadth=h.get('backbone_breadth',''),clean_fragments_denominator=m['clean_fragments_denominator'],
            catalog_version=m['catalog_version'],mapping_reference_version=m['mapping_reference_version'],
            detection_status='not_detected_at_this_depth' if d['detection_status']=='not_detected' else d['detection_status'],
            lineage='unknown',sample_host_phasing='not_established'))
    sf=['unit_id','biological_sample_id','subject_id','reference_host_id','quantification_group','detection_profile_id','status','reason','value','fragment_rpkm_all_clean','host_core_depth','host_core_breadth','clean_fragments_denominator','catalog_version','mapping_reference_version','detection_status','lineage','sample_host_phasing']
    write_table(out/'community_signals.tsv',sf,signal)
    write_table(out/'unit_status.tsv',['unit_id','measurement'],[dict(unit_id=u,measurement='completed') for u in units])
    (out/'OUTPUTS.md').write_text('Research community measurement uses one sequence-equivalence quantification_group per reference, never one independent observation per mapped vOTU. host_core_depth is competitive depth on filtered canonical core markers. It is not MetaPhlAn relative abundance. Detected signal does not establish physical host linkage, induction or copy number. Ambiguous-only and missing host evidence stay unassessable. Science thresholds remain uncalibrated.\n')
    (out/'report_cohort.html').write_text('<!doctype html><meta charset="utf-8"><title>Research community measurement</title><h1>Quantification-group measurement</h1><p>Engineering output; scientific calibration pending.</p><p>'+str(len(units))+' units; '+str(len(mapping))+' groups.</p><a href="community_signals.tsv">Statistical interface</a> · <a href="OUTPUTS.md">Definitions</a>')
    write_json(out/'cohort_manifest.json',dict(schema_version=1,stage='B1_research_community',status='completed',measurement_contract='research_quantification_group_v1',reference_manifest_sha256=sha256(root/'reference/reference_manifest.json'),
        mapping_reference_version=ref['mapping_reference_version'],catalog_version=ref['catalog_version'],unit_manifests=manifests,
        group_measurements=len(measurements),identity_numeric_detection_validation='passed',biological_votu_expansion=False,
        outputs={p.name:sha256(p) for p in out.iterdir() if p.is_file()},scientific_calibration=False))
    return signal
