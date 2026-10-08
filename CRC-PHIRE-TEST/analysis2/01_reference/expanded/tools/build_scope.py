"""Write explicit evidence-led search scope. Caps are metadata screening budgets, not power targets."""
import json
from pathlib import Path
B=Path(__file__).resolve().parents[1]
sources={
'CRC2025':'https://www.nature.com/articles/s41591-025-03693-9',
'CRC2019':'https://www.nature.com/articles/s41591-019-0405-7',
'QMP2024':'https://www.nature.com/articles/s41591-024-02963-2',
'VIROME2023':'https://pmc.ncbi.nlm.nih.gov/articles/PMC10334131/',
'FNA2024':'https://www.nature.com/articles/s41586-024-07182-w',
'FNTAX2025':'https://journals.asm.org/doi/10.1128/mbio.00941-25',
'BF2026':'https://www.nature.com/articles/s43856-026-01403-1',
'PROPHAGE2025':'https://www.nature.com/articles/s41467-025-66733-5',
'GPD2021':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7895897/',
'BUTYRICICOCCUS2020':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7577080/',
'AKK2021':'https://www.nature.com/articles/s41598-021-82465-0',
'FPER2026':'https://www.nature.com/articles/s41467-026-74591-y',
'HIBC2025':'https://www.nature.com/articles/s41467-025-59229-9',
'PGENOTOX2022':'https://www.nature.com/articles/s41586-022-04444-3',
}
# query, role, direction, cap, citations, interpretation. Taxonomy returned by API retained independently.
specs=[
('Fusobacterium animalis','deep','CRC_enriched_lineage_dependent',180,'FNA2024;FNTAX2025;CRC2025','保留旧名与新种/谱系对应；不把旧Fna C1/C2当永恒分类'),
('Fusobacterium nucleatum','deep','CRC_enriched_lineage_dependent',120,'CRC2025;FNTAX2025','查询包含子类群时全局去重复；宿主/口腔与肠道来源分开'),
('Fusobacterium polymorphum','lineage_context','CRC_enriched_lineage_dependent',70,'CRC2025;FNTAX2025','Fn复合群的谱系覆盖'),
('Fusobacterium vincentii','lineage_context','CRC_enriched_lineage_dependent',70,'CRC2025;FNTAX2025','Fn复合群的谱系覆盖'),
('Fusobacterium watanabei','taxonomy_context','unresolved_taxonomic_context',40,'FNTAX2025','新分类参考；不外推疾病方向'),
('Fusobacterium paranimalis','taxonomy_context','unresolved_taxonomic_context',40,'FNTAX2025','保留分类争议和类型株锚点'),
('Fusobacterium periodonticum','focused','CRC_association_with_recent_mechanistic_support',60,'FPER2026','新增研究线索，不等同多队列普适结论'),
('Parvimonas micra','deep','CRC_enriched',180,'CRC2025;QMP2024','多来源种内prophage覆盖；CRC来源不是参考纳入硬条件'),
('Peptostreptococcus anaerobius','deep','CRC_enriched',150,'QMP2024;CRC2019','重复关联和机制关注宿主'),
('Peptostreptococcus stomatis','deep','CRC_enriched',120,'CRC2025;CRC2019','不与P. anaerobius合并；公开培养genome可能较少'),
('Bacteroides fragilis','deep','strain_dependent_mixed',200,'BF2026;CRC2025','ETBF/NTBF依基因证据；不把整个物种等同致病'),
('Bacteroides hominis','lineage_context','unresolved_taxonomic_context',80,'BF2026','历史division II单列；当前taxonomy与旧标签核查'),
('Gemella morbillorum','focused','CRC_enriched',100,'CRC2025;CRC2019','口腔/肠道/临床体位分开'),
('Solobacterium moorei','focused','CRC_enriched',80,'CRC2019;CRC2025','补足反复出现的口腔相关CRC宿主'),
('Porphyromonas asaccharolytica','focused','CRC_enriched',80,'QMP2024;CRC2019','定量混杂调整后仍有支持'),
('Porphyromonas somerae','focused','CRC_association_context_dependent',60,'CRC2019;CRC2025','不将属级方向直接推广所有种'),
('Prevotella intermedia','focused','CRC_enriched',80,'QMP2024;CRC2019','口腔来源可作参考，不能当健康人对照'),
('Anaerococcus vaginalis','focused','CRC_enriched',60,'QMP2024','补足混杂调整后支持的宿主'),
('Dialister pneumosintes','focused','CRC_enriched',60,'QMP2024;CRC2025','补足混杂调整后支持的宿主'),
('Hungatella hathewayi','focused','CRC_enriched',100,'CRC2025;VIROME2023','种内与旧Clostridium标签交叉核查'),
('Clostridium symbiosum','focused','CRC_association_context_dependent',80,'CRC2019;CRC2025','数据库新旧分类分别保留'),
('Streptococcus gallolyticus','focused','CRC_relevant_clinical_association',100,'CRC2019','临床宿主资源；血培养不等同粪便病例对照'),
('Faecalibacterium','deep_depleted_background','control_associated_subset_direction_not_genus_wide',200,'VIROME2023;CRC2019;HIBC2025','属级检索覆盖F. prausnitzii复合群的新种；逐种解释'),
('Roseburia','deep_depleted_background','control_associated_subset_direction_not_genus_wide',160,'VIROME2023;CRC2019','属级检索补种内/种间覆盖，不宣称全属保护'),
('Agathobacter rectalis','depleted_background','control_associated',150,'VIROME2023','旧Eubacterium rectale名称对应；耗竭不等同保护'),
('Anaerostipes hadrus','ecological_context','context_dependent',100,'CRC2025','肠道功能/谱系背景；不预设CRC方向'),
('Coprococcus comes','ecological_context','direction_not_prespecified',60,'GPD2021','常见肠道发酵生态背景，未据此宣称CRC耗竭'),
('Bifidobacterium longum','ecological_context','control_association_context_dependent',100,'CRC2019','种/亚种和年龄来源需区分，不预设保护'),
('Bifidobacterium adolescentis','ecological_context','direction_not_prespecified',80,'GPD2021','补充非CRC富集的肠道宿主背景'),
('Akkermansia muciniphila','ecological_context','context_dependent',100,'AKK2021;CRC2025','CRC中可富集，不能自动列为保护菌'),
('Bacteroides uniformis','ecological_context','direction_not_prespecified',80,'GPD2021','常见肠道Bacteroides背景，避免只看B. fragilis'),
('Bacteroides thetaiotaomicron','ecological_context','direction_not_prespecified',80,'GPD2021','可实验宿主及功能背景，不预设疾病方向'),
('Phocaeicola vulgatus','ecological_context','context_dependent',80,'CRC2019','旧Bacteroides vulgatus；高丰度近缘背景'),
('Phocaeicola dorei','ecological_context','context_dependent',60,'CRC2019','旧Bacteroides dorei；与vulgatus分开'),
('Escherichia coli','mechanistic_context','strain_gene_dependent',150,'CRC2025;PGENOTOX2022','仅RefSeq完整组装盘点作为有界结构背景；不是E. coli全物种穷尽检索'),
('Butyricicoccus pullicaecorum','experimental_context','potential_protective_preclinical_only',30,'BUTYRICICOCCUS2020','动物实验依据；非人类型株只作实验/分类背景，单列不伪装患者宿主'),
]
taxa=[]
for n,role,direction,cap,refs,note in specs:
    q={'key':n.lower().replace(' ','_'),'query':n,'role':role,'disease_direction':direction,
       'candidate_budget_cap':cap,'source_ids':refs.split(';'),'interpretation':note,
       'scope_rank':'genus' if ' ' not in n else 'named_taxon','api_filters':{}}
    if n=='Escherichia coli':q['api_filters']={'filters.assembly_source':'refseq','filters.assembly_level':'complete_genome'}
    taxa.append(q)
cfg={'date':'2026-10-08','stage':'metadata_preselection_not_final_sequence_QC','taxa':taxa,'sources':sources,
     'selection_rules':{'source':'NCBI Datasets v2 current taxonomy; exclude API-labelled MAGs',
     'metadata_QC':'completeness >=90 if reported; contamination <=5 if reported; missing values retained as conditional',
     'fragmentation':'<=200 contigs preferred; 201-500 reserve unless otherwise justified; >500 excluded from main preselection',
     'priority':'human gut then other human then unknown; nonhuman normally exclude except type/experimental anchors separate',
     'duplicate':'paired GCA/GCF one preferred assembly; same BioSample one best assembly; named-strain aliases audited separately',
     'representativity':'rank quality within BioProject/country/body-site strata then round-robin; preserve type anchors; no claimed sequence phylogenetic dereplication',
     'old258':'no automatic retention privilege; current metadata and same rules; individually reconciled',
     'missing_checkm':'not equivalent to passing QC; Jinlong evaluates from sequence',
     'no_new_FASTA_here':True,'no_disease_pvalue_selection':True}}
B.mkdir(parents=True,exist_ok=True)
(B/'selection_scope.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf-8')
print('search_units',len(taxa),'budget_ceiling',sum(t['candidate_budget_cap'] for t in taxa))
