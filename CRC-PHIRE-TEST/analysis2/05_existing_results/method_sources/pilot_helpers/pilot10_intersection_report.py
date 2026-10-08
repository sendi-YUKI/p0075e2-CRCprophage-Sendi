"""Direct G-M sequence evidence, functional links and descriptive 5 vs 5 report."""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import csv,hashlib,html,json,os,statistics
R=Path('/srv/CRC-PHIRE/analysis1')
def read(p):return json.loads(p.read_text())
def table(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,rows,fields):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def pair_metrics(g,m,lengths,directions):
    if g==m:return dict(ani_percent=100.0,g_coverage_percent=100.0,m_coverage_percent=100.0,match='identical_oriented_sequence_sha256')
    a=directions.get((g,m));b=directions.get((m,g))
    if not a or not b:return None
    ani=min(float(a['ani_percent']),float(b['ani_percent']))
    gc=min(float(a['query_covered_bp']),float(b['subject_covered_bp']))*100/lengths[g]
    mc=min(float(a['subject_covered_bp']),float(b['query_covered_bp']))*100/lengths[m]
    short=gc if lengths[g]<lengths[m] else mc if lengths[m]<lengths[g] else min(gc,mc)
    if ani<95 or short<85:return None
    return dict(ani_percent=ani,g_coverage_percent=gc,m_coverage_percent=mc,
        match='broad_bidirectional_match' if gc>=85 and mc>=85 else 'partial_or_fragment_match')
def group_descriptions(donor,groups):
    """Retain every measured vOTU; missing/unevaluable is never a negative donor."""
    assert len(groups)==10 and Counter(groups.values())=={'CRC':5,'control':5}
    bytarget=defaultdict(list)
    for r in donor:
        assert groups[r['biological_sample_id']]==r['group']
        bytarget[r['reference_id']].append(r)
    output=[]
    for target,rows in sorted(bytarget.items()):
        assert len(rows)==10 and {r['biological_sample_id'] for r in rows}==set(groups)
        first=rows[0]
        for field in ['catalog_version','mapping_reference_version','votu_id','candidate_id','target_length']:
            assert len({r[field] for r in rows})==1,(target,field)
        summary={k:first[k] for k in ['catalog_version','mapping_reference_version','votu_id','candidate_id','target_length']}
        summary['reference_id']=target
        for group,prefix in [('CRC','crc'),('control','control')]:
            selected=[r for r in rows if r['group']==group]
            measured=[r for r in selected if r['measurement_status']=='assessed']
            evaluable=[r for r in selected if r['detection_status'] in {'detected','not_detected'}]
            detected=sum(r['detection_status']=='detected' for r in evaluable)
            summary.update({prefix+'_donors':5,prefix+'_measurement_assessed':len(measured),
                prefix+'_detection_evaluable':len(evaluable),prefix+'_detected':detected,
                prefix+'_detection_fraction':detected/len(evaluable) if evaluable else None})
            for measure in ['fragment_rpkm_all_clean','mean_depth_full','breadth_full']:
                summary[prefix+'_median_'+measure]=statistics.median(float(r[measure]) for r in measured) if measured else None
        a=summary['crc_detection_fraction'];b=summary['control_detection_fraction']
        summary.update(crc_minus_control_detection_fraction=a-b if a is not None and b is not None else None,
            interpretation='descriptive_same_discovery_samples;not_host_adjusted;no_p_or_q_values')
        output.append(summary)
    return output
def main():
    assert os.environ.get('SLURM_JOB_ID')
    root=R/'06_viral_catalog/pilot10_combined';cat=root/'catalog';out=R/'12_reports/pilot10';out.mkdir(exist_ok=True)
    accepted=read(R/'13_handoff/pilot10_selective_publication_acceptance.json')
    assert accepted['status']=='completed_catalog_selective_publication'
    assert sha(root/'pipeline_info/terminal_output_manifest.json')==accepted['terminal_manifest_sha256']
    for rel,digest in read(root/'pipeline_info/terminal_output_manifest.json').items():assert sha(root/rel)==digest
    measured=read(out/'donor_measurement_acceptance.json');assert measured['status']=='completed' and measured['donors']==10
    for rel,digest in measured['sha256'].items():assert sha(out/rel)==digest
    occ=table(cat/'source_occurrences.tsv');masters=table(cat/'viral_sequence_master.tsv')
    assert all(o['slice_verification']=='verified_this_run' for o in occ)
    seq={r['viral_sequence_id']:r for r in masters};lengths={k:int(r['length']) for k,r in seq.items()}
    G=defaultdict(list);M=defaultdict(list)
    for o in occ:
        (G if o['source_kind']=='legacy_bacterial_genome' else M)[o['viral_sequence_id']].append(o)
    assert sum(map(len,G.values()))==812
    frozen=read(R/'01_manifests/pilot10_frozen_assemblies.json');frozen_ids=set(frozen['groups'])
    assert len(frozen_ids)==10 and Counter(frozen['groups'].values())=={'CRC':5,'control':5}
    assert {o['biological_sample_id'] for v in M.values() for o in v}<=frozen_ids
    assert {u['biological_sample_id'] for u in read(cat/'discovery_manifest.json')['units']}==frozen_ids
    directions={(r['query_id'],r['subject_id']):r for r in table(cat/'votu_directional_ani.tsv')}
    possible={(g,m) for g,m in directions if g in G and m in M}
    possible|={(g,m) for m,g in directions if g in G and m in M}
    possible|={(s,s) for s in set(G)&set(M)}
    metadata={r['accession']:r for r in table(R/'01_manifests/reference_panel_metadata_original.tsv')}
    links=[];sequence_links=[]
    for g,m in sorted(possible):
        metrics=pair_metrics(g,m,lengths,directions)
        if metrics is None:continue
        common=dict(g_sequence_id=g,m_sequence_id=m,**metrics,g_length=lengths[g],m_length=lengths[m],
            both_length_ge5000=lengths[g]>=5000 and lengths[m]>=5000,
            g_quality_eligibility=seq[g]['viral_catalog_eligibility'],m_quality_eligibility=seq[m]['viral_catalog_eligibility'])
        sequence_links.append(common)
        for go in G[g]:
            meta=metadata.get(go['assembly_id'],{})
            for mo in M[m]:
                links.append(dict(common,g_occurrence_id=go['source_occurrence_id'],m_occurrence_id=mo['source_occurrence_id'],
                    reference_genome=go['assembly_id'],reference_strain=meta.get('canonical_strain_name','unknown'),
                    reference_taxon=meta.get('adopted_taxon','unknown'),patient_sample=mo['biological_sample_id'],
                    g_candidate_id=go['candidate_id'],m_candidate_id=mo['candidate_id'],
                    g_candidate_type=go['candidate_type'],m_candidate_type=mo['candidate_type'],
                    g_integration_evidence=go['integration_evidence'],m_integration_evidence=mo['integration_evidence'],
                    g_contig=go['contig_id'],g_start0=go['start0'],g_end0=go['end0'],
                    m_contig=mo['contig_id'],m_start0=mo['start0'],m_end0=mo['end0'],
                    patient_host_claim='reference_host_link_only_not_confirmed_patient_host',
                    cargo_identity_claim='not_inferred_from_sequence_cluster_or_global_match'))
    fields=['g_sequence_id','m_sequence_id','ani_percent','g_coverage_percent','m_coverage_percent','match','g_length','m_length',
        'both_length_ge5000','g_quality_eligibility','m_quality_eligibility']
    write(out/'G_M_sequence_matches.tsv',sequence_links,fields)
    write(out/'G_M_occurrence_host_links.tsv',links,fields+['g_occurrence_id','m_occurrence_id','reference_genome','reference_strain','reference_taxon','patient_sample','g_candidate_id','m_candidate_id','g_candidate_type','m_candidate_type','g_integration_evidence','m_integration_evidence','g_contig','g_start0','g_end0','m_contig','m_start0','m_end0','patient_host_claim','cargo_identity_claim'])
    matchedG={r['g_sequence_id'] for r in sequence_links if r['match']!='partial_or_fragment_match'}
    matchedM={r['m_sequence_id'] for r in sequence_links if r['match']!='partial_or_fragment_match'}
    dispositions=[dict(viral_sequence_id=s,source_side=side,direct_broad_match=s in matched,
        quality_eligibility=seq[s]['viral_catalog_eligibility'],length=lengths[s]) for side,values,matched in [('G',G,matchedG),('M',M,matchedM)] for s in sorted(values)]
    write(out/'G_M_sequence_disposition.tsv',dispositions,['viral_sequence_id','source_side','direct_broad_match','quality_eligibility','length'])
    funcs=root/'functions';gene_to_seq={r['gene_id']:r['viral_sequence_id'] for r in table(root/'genes/gene_table.tsv')}
    feature_rows=[]
    for r in table(funcs/'function_hits.tsv'):
        # Raw PHROGs output also contains rejected / non-best domains.
        if r.get('accepted')!='true' or r.get('is_best')!='true':continue
        sid=gene_to_seq[r['gene_id']]
        if sid in matchedG|matchedM:feature_rows.append(dict(viral_sequence_id=sid,gene_id=r['gene_id'],type='PHROG_homology',feature=r.get('phrog',''),evidence_file=str(funcs/'function_hits.tsv'),interpretation='homology_not_experimental_function'))
    for r in table(funcs/'systems.tsv'):
        sid=r.get('viral_sequence_id') or r['contig_id']
        if sid in matchedG|matchedM:feature_rows.append(dict(viral_sequence_id=sid,gene_id='',type=r['activity'],feature=r['subtype'],evidence_file=str(funcs/'systems.tsv'),interpretation='viral_sequence_context_only_not_whole_host_defense'))
    ko=root/'viral_kofam/kofam/viral_member_ko_counts.tsv'
    for r in table(ko):
        if r['viral_sequence_id'] in matchedG|matchedM:feature_rows.append(dict(viral_sequence_id=r['viral_sequence_id'],gene_id=r['gene_ids'],type='KO_homology',feature=r['ko'],evidence_file=str(ko),interpretation='candidate_cargo_requires_context_review_not_confirmed_AMG'))
    write(out/'intersection_function_evidence.tsv',feature_rows,['viral_sequence_id','gene_id','type','feature','evidence_file','interpretation'])
    donor=table(out/'donor_votu_measurement.tsv')
    descriptive=group_descriptions(donor,frozen['groups'])
    descriptive_fields=['catalog_version','mapping_reference_version','votu_id','candidate_id','target_length','reference_id']
    for prefix in ['crc','control']:
        descriptive_fields += [prefix+'_'+k for k in ['donors','measurement_assessed','detection_evaluable','detected','detection_fraction','median_fragment_rpkm_all_clean','median_mean_depth_full','median_breadth_full']]
    descriptive_fields += ['crc_minus_control_detection_fraction','interpretation']
    write(out/'pilot10_descriptive_group_comparison.tsv',descriptive,descriptive_fields)
    report=dict(status='completed_balanced_engineering_pilot',created_utc=datetime.now(timezone.utc).isoformat(),
        reference_genomes=258,donors=10,groups=dict(Counter(frozen['groups'].values())),reference_raw_candidates=812,
        reference_raw_provirus_loci=541,metagenome_raw_occurrences=sum(map(len,M.values())),
        metagenome_occurrence_types=dict(Counter(o['candidate_type'] for v in M.values() for o in v)),
        combined_unique_sequences=len(seq),engineering_quality_votus=read(cat/'catalog_manifest.json')['representatives'],
        direct_sequence_matches=dict(Counter(r['match'] for r in sequence_links)),
        broad_matched_reference_sequences=len(matchedG),broad_matched_metagenome_sequences=len(matchedM),
        occurrence_host_links=len(links),donor_broad_detections=dict(Counter(r['biological_sample_id'] for r in donor if r['detection_status']=='detected')),
        disease_association='descriptive_5_vs_5_only_no_significance_test_single_cohort_discovery_samples',
        descriptive_group_features=len(descriptive),
        patient_host_assignment='not_proven_by_reference_intersection',
        scope='Balanced ten-donor technical and descriptive pilot; not completion of nine-cohort study or independent validation')
    paragraphs=[
        '# CRC-PHIRE：258份参考＋10人先导补充结果',
        '本报告只在十人联合目录、功能注释、十人统一定量及逐文件验收全部完成后生成；原六人结果单独保留。',
        '## 实际完成范围',
        '参考端258份基因组；患者端为Feng队列5例CRC、5名研究对照。保留原先导六人，按既有清单顺序补入SID530368、SID532796、SID530295、SID31512四名对照；SID31285另行保留、不纳入。四名新增对照在本轮分析结果产生前固定，没有根据噬菌体结果选人。研究对照可有高血压或脂肪肝等共病，不等同于全部健康人。',
        f"参考端812条原始候选中有541个预测provirus位点。患者端得到{report['metagenome_raw_occurrences']}个原始病毒实例；原始病毒实例不能全部称为prophage。联合目录包含{len(seq)}条按序列SHA去重的序列，工程质量集合形成{report['engineering_quality_votus']}个vOTU。",
        '## 两端能否连接',
        f"直接匹配分类：{json.dumps(report['direct_sequence_matches'],ensure_ascii=False)}。有广泛直接匹配的参考端序列{len(matchedG)}条、患者端序列{len(matchedM)}条。",
        '匹配使用既定95% ANI及85%覆盖起点；广泛匹配要求两条序列都覆盖至少85%，片段匹配单列。对两个方向的非重叠比对证据取保守值。相同方向序列SHA完全相同另列。',
        'G_M_occurrence_host_links.tsv给出菌株/基因组、患者样本、两端候选编号和坐标。这个连接证明患者序列与参考菌株中的病毒序列相似；不能据此确认它在患者体内属于该菌，也不能推断功能基因完全一致。',
        '## 功能与定量怎么阅读',
        '四名新增对照同样参加组装和病毒发现，十人共同建立新版目录；十人全部向同一新版参考重新定量，未把六人旧目录的数值与四人新目录的数值直接拼接。258份参考和已有六人的主要预测结果经校验复用。',
        '功能从各条病毒序列自身预测基因后统一进行PHROGs、DefenseFinder/AntiDefenseFinder和KOfam搜索。intersection_function_evidence.tsv汇总匹配对象的证据链接。KO、cargo或反防御同源命中仍需检查基因上下文、边界和系统完整性，不是机制验证。',
        '每人PE、SE和fastp孤儿reads分别比对，再按碱基深度相加计算同一个人的覆盖广度；没有把readsets当作独立患者。竞争参考加入258份细菌基因组并遮蔽目录中的病毒位点，减少明显的细菌背景误配，但不是完整肠道背景库。',
        '广泛检出采用整条病毒参考至少75%碱基在过滤后达到至少1×覆盖；保留连续深度、广度和歧义信号。未达阈值不等于不存在。258份单独菌株的竞争比对深度不能直接当作已验证的物种丰度或宿主校正结果。',
        '## 这次可以与不可以下的结论',
        '可以评价：现有数据是否支持两端序列连接、功能证据能否追踪到具体序列、按人定量能否贯通。若未得到广泛交集，保留片段匹配、目录覆盖和宿主检出背景，不能为了取出交集降低阈值。',
        'pilot10_descriptive_group_comparison.tsv按同一vOTU列出两组的检出人数、可评价人数、检出比例及丰度、深度、广度中位数。未评价状态不算阴性，分母单独列出；所有已测量vOTU均保留。这里只描述这十人的观察差异，没有进行显著性检验、P值筛选或疾病预测。',
        '十人同时参与发现与测量，仍是样本内先导。5∶5增加了对照侧信息，但样本少、单队列、年龄和共病等尚不能充分控制，不证明CRC关联、致病或保护作用，也不是独立验证。两组都是4女1男，不代表已经完成年龄或其他临床因素匹配；原始临床元数据的核查状态仍须保留。正式研究与analysis2交接范围需独立推进。',
        '## 文件入口',
        '- G_M_sequence_matches.tsv：直接序列匹配与双向覆盖。\n- G_M_occurrence_host_links.tsv：参考菌株与患者病毒实例对应。\n- G_M_sequence_disposition.tsv：匹配与未匹配序列均保留。\n- intersection_function_evidence.tsv：匹配对象的功能证据。\n- donor_votu_measurement.tsv：十人的定量、覆盖和检出状态。\n- pilot10_descriptive_group_comparison.tsv：5∶5的描述性比较。\n- pilot10_final_acceptance.json：本报告数值与文件校验。',
    ]
    md='\n\n'.join(paragraphs)+'\n';(out/'十人先导补充结果说明.md').write_text(md,encoding='utf-8')
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CRC-PHIRE 十人先导补充</title><style>body{max-width:1000px;margin:40px auto;padding:0 24px;color:#183142;font:17px/1.9 system-ui}article{white-space:pre-wrap}h1{font-size:28px}</style><h1>CRC-PHIRE 十人先导补充</h1><article>'+html.escape(md)+'</article></html>'
    (out/'十人先导补充结果说明.html').write_text(page,encoding='utf-8')
    report['output_sha256']={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='pilot10_final_acceptance.json'}
    (out/'pilot10_final_acceptance.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
