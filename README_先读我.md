# CRC-PHIRE：给 Jinlong 的正式研究交接

CRC-PHIRE = **Colorectal Cancer Prophage–Host Interactions and Repertoire Exploration**。

本交接包更新于2026-10-08。目标是从多队列CRC与合格对照粪便宏基因组发现可重复的prophage信号，以细菌参考基因组解释宿主、整合结构和功能，最后选择1–2个可实验对象。

**Jinlong直接开展正式研究，不再另做先导。** 已有Feng队列5 CRC＋5对照pilot及后诊断供了解经验；是否复用任何输入、中间结果或代码由Jinlong决定。工具、数据库、参数、并行和服务器路径按Jinlong的成熟pipeline选定并记录。

## 建议阅读顺序

1. [全局研究计划](00_project/CRC-PHIRE_全局研究计划.md)：背景、问题、分工、完整分析和冻结时机。
2. [Jinlong执行说明](00_project/CRC-PHIRE_Jinlong执行说明.md)：上机后按什么顺序做。
3. [参考预选说明与数量](01_reference/expanded/参考集合选择说明.md)和[可下载候选清单](01_reference/expanded/download_manifest.tsv)。
4. [已有pilot方法与结果](05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)：真实数字、图、代码和后续用途。
5. [方法与参数问题表](03_execution/Jinlong_方法与参数决策表.md)：重点回答影响科学结论的选择，并把全部实际参数导出。
6. [九队列说明](02_cohorts/九队列正式研究说明.md)及[交付与冻结要求](04_deliverables/交付清单与表字段.md)。

核心MD旁有同名HTML，下载仓库后可在浏览器阅读；GitHub网页直接阅读MD。

## 目前拿到了什么

- 公共检索返回21,608条组装记录，身份整理为15,064个元数据组；1,931份预选组装，覆盖36个检索单元；按NCBI当前标签分为54个具名种标签及22份未定种记录。
- 预选中177个身份与旧258对应，1,754个为新增参考身份。旧258逐份保留、替换、备用或待查见[对账表](01_reference/expanded/old258_reconciliation.tsv)。
- 这些是**可以开始下载和统一QC的候选**，不是最终通过序列QC的菌株数量；命名别名和近克隆仍需Jinlong核查。
- 已有258 FASTA保留本地，GitHub不存这约1.03 GB文件；可按原accession下载或另行传输。候选清单中新增FASTA尚未批量下载。
- 已有pilot：258参考、10人、541＋382个预测provirus实例、473个定量vOTU；有G—M序列连接和2个复制活动候选。小样本关联未通过全局FDR。

## 分工与两个冻结点

本地Part 0已整理设计、公共元数据预选、pilot证据与执行记录模板。Jinlong在Part 1并行处理参考FASTA与九队列患者上游。**先冻结正式参考accession版本和数量，再冻结由G与开发患者M构成的联合病毒目录，之后进入Part 2联合定量、关联和保留队列验证。** 不必等待参考最终数量才开始患者QC和组装。

公开包不含服务器登录信息、患者FASTQ、BAM、数据库或原始大缓存。详见[公开交付范围](90_package/GitHub交付范围.md)。历史材料仅是证据，当前执行要求以本入口及新计划为准。
