import pathlib,json,csv,collections,math,hashlib,itertools
import numpy as np
from scipy.stats import fisher_exact,binom
R=pathlib.Path('/srv/CRC-PHIRE/analysis1');A=R/'14_postpilot_adjustment';P=R/'12_reports/pilot10'
def table(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def save(p,d):p.write_text(json.dumps(d,indent=2))
def bh(p):
    p=np.array(p);order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty(len(p));out[order]=np.minimum(q,1);return out
def main():
    out=A/'04_measurement_audit';m=table(P/'donor_votu_measurement.tsv');d=table(P/'exploratory/votu_differential.tsv');occ=table(P/'prophage_evidence/occurrence_evidence.tsv');metadata=table(R/'01_manifests/reference_panel_metadata_original.tsv')
    assert len(m)==4730 and len(d)==473 and len(metadata)==258
    fp=[float(r['fisher_p']) for r in d if r['fisher_p'] not in ['NA','']];ap=[float(r['abundance_p']) for r in d if r['abundance_p'] not in ['NA','']]
    changes=[];counts=[]
    for threshold in [.75,.30]:
        for r in m:
            assigned=int(float(r['assigned_fragment_count']));amb=int(float(r['ambiguous_fragment_count']));b=float(r['breadth_full'])
            if r['measurement_status']!='assessed':status='unassessable';value='NA'
            elif b>=threshold and assigned>=1:status='detected';value=1
            elif amb and assigned==0:status='not_resolvable';value='NA'
            elif amb:status='not_detected_at_resolved_signal';value='NA'
            else:status='not_detected';value=0
            if threshold==.75:assert status==r['detection_status'],(status,r['detection_status'])
            changes.append({'sample':r['biological_sample_id'],'group':r['group'],'votu_id':r['votu_id'],'breadth_threshold':threshold,'breadth_full':b,'assigned_fragments':assigned,'ambiguous_target_mentions':amb,'status':status,'value':value,'baseline_status':r['detection_status']})
        for group in ['CRC','control']:
            c=collections.Counter(r['status'] for r in changes if r['breadth_threshold']==threshold and r['group']==group)
            for k,v in c.items():counts.append({'threshold':threshold,'group':group,'status':k,'sample_votu_cells':v})
    write(out/'threshold_states.tsv',changes);write(out/'threshold_summary.tsv',counts)
    meta=[r for r in occ if r['source_kind']=='metagenome_assembly' and r['candidate_type']=='provirus_locus'];pro=[r for r in occ if r['candidate_type']=='provirus_locus']
    floor=2/math.comb(10,5);counter=bh([floor]*130+[1]*685)
    power=[]
    for n in [5,25,50]:
        for alpha in [.05,5e-5]:
            p=0.
            for x in range(n+1):
                for y in range(n+1):
                    if fisher_exact([[x,n-x],[y,n-y]]).pvalue<alpha:p+=binom.pmf(x,n,.224)*binom.pmf(y,n,.026)
            power.append({'n_per_group':n,'alpha':alpha,'assumed_case_prevalence':.224,'assumed_control_prevalence':.026,'exact_power':p,'scope':'illustrative_species_prevalence_NOT_prophage_power'})
    write(A/'06_statistics_design/exact_power_correction.tsv',power)
    facts={'rows_measurement':len(m),'votus':len(d),'fisher_tests':len(fp),'abundance_tests':len(ap),'fisher_nominal':sum(p<.05 for p in fp),'abundance_nominal':sum(p<.05 for p in ap),'global_BH_min':float(bh(fp+ap).min()),'endpoint_BH_abundance_min':float(bh(ap).min()),'p_floor_5v5':floor,'BH_counterexample_n_tests':815,'BH_counterexample_n_floor_p':130,'BH_counterexample_min_q':float(counter.min()),'measurement_status_counts':dict(collections.Counter(r['detection_status'] for r in m)),'reference_taxa':dict(collections.Counter(r['adopted_taxon'] for r in metadata)),'eligible_isolate_comparison':sum(r['eligible_for_disease_comparison']=='true' for r in metadata),'unknown_clinical':sum(r['clinical_label_or_unknown']=='unknown' for r in metadata),'predicted_provirus_total':len(pro),'meta_provirus':len(meta),'two_flanks_total':sum(int(r['left_flank_bp'])>0 and int(r['right_flank_bp'])>0 for r in pro),'two_flanks_meta':sum(int(r['left_flank_bp'])>0 and int(r['right_flank_bp'])>0 for r in meta),'two_1kb_flanks_meta':sum(int(r['left_flank_bp'])>=1000 and int(r['right_flank_bp'])>=1000 for r in meta),'two_hostmarkers_meta':sum(r['two_sided_host_marker_evidence']=='true' for r in meta)}
    detected=collections.Counter(r['votu_id'] for r in m if r['detection_status']=='detected');facts['n_detected_donor_distribution']=dict(collections.Counter(detected[r['votu_id']] for r in d));save(out/'facts.json',facts)
    # Count actual ambiguity fragments once per unique readset, not once per target.
    ref=table(R/'09_measurement/pilot10/prepared/reference/reference_map.tsv');types={r['reference_id']:r['reference_kind'] for r in ref};votu={r['reference_id']:r['votu_id'] for r in ref if r['reference_kind']=='target'}
    units=json.loads((R/'01_manifests/pilot10_clean_reads.json').read_text())['units'];summary=[];targetrows=[];paircounts=collections.Counter();filters=[]
    for u in units:
        directory=R/'09_measurement/pilot10/units'/('unit_'+u['unit_id']);c=collections.Counter();tc=collections.Counter();manifest=json.loads((directory/'measurement_manifest.json').read_text())
        for k,v in manifest['fragment_filter_counts'].items():filters.append({'sample':u['biological_sample_id'],'unit':u['unit_id'],'layout':u['layout'],'filter':k,'fragments':v})
        with open(directory/'ambiguity.tsv') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                ids=r['reference_ids'].split(',');ts=[x for x in ids if types[x]=='target'];hs=[x for x in ids if types[x]=='host']
                kind='host_only' if not ts else ('target_host' if hs else ('multiple_targets' if len(ts)>1 else 'same_target_locations_or_saturation'))
                c[kind]+=1
                for target in ts:tc[(target,kind)]+=1
                if len(ts)>1:
                    for x,y in itertools.combinations(sorted(ts),2):paircounts[(x,y)]+=1
        assert sum(c.values())==manifest['fragment_filter_counts'].get('ambiguous_or_reporting_saturated',0)
        for k,v in c.items():summary.append({'sample':u['biological_sample_id'],'unit':u['unit_id'],'layout':u['layout'],'group':u['metadata']['disease_group'],'ambiguity_class':k,'unique_within_readset_fragments':v})
        for (target,k),v in tc.items():targetrows.append({'sample':u['biological_sample_id'],'unit':u['unit_id'],'reference_id':target,'votu_id':votu[target],'ambiguity_class':k,'fragment_mentions':v})
    write(out/'ambiguity_fragment_classes.tsv',summary);write(out/'ambiguity_target_mentions.tsv',targetrows);write(out/'read_filter_counts.tsv',filters)
    if paircounts:write(out/'target_ambiguity_pairs.tsv',[{'target_1':x,'target_2':y,'votu_1':votu[x],'votu_2':votu[y],'shared_fragment_mentions':n} for (x,y),n in paircounts.most_common()])
    # Reference-origin linkage is a context label, not patient host identification.
    accession={r['accession']:r['adopted_taxon'] for r in metadata};links=collections.defaultdict(set)
    for r in occ:
        if r['votu_id'] and r['assembly_id'] in accession:links[r['votu_id']].add(accession[r['assembly_id']])
    write(out/'votu_reference_taxa.tsv',[{'votu_id':r['votu_id'],'reference_taxa':';'.join(sorted(links[r['votu_id']])),'has_metagenome_predicted_provirus':r['has_metagenome_predicted_provirus'],'has_reference_predicted_provirus':r['has_reference_predicted_provirus']} for r in d])
    save(out/'completed.json',{'status':'completed','fact_summary':facts,'ambiguity_classes':dict(sum((collections.Counter({r['ambiguity_class']:r['unique_within_readset_fragments']}) for r in summary),collections.Counter())),'threshold_policy':'retain_0.75_primary;0.30_sensitivity_only_until_specificity_supported','known_negative_control_assessed':False,'new_read_allocation':False,'analysis_scope':'measurement_diagnostics_not_new_association_discovery'})
if __name__=='__main__':main()
