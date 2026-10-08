# 给 Jinlong 的方法与参数决策表

请按自己的成熟pipeline决定工具、数据库和数值参数。这里不是要求逐项征求研究者审批，也不要求新先导；目的是让研究者最后清楚“实际做了什么、参数在哪里改、改后影响哪些步骤”。主要选择写理由，其余完整参数从真实配置/命令导出即可。

同名TSV提供可填写列：旧pilot值、建议、Jinlong决定、真实运行值、单位/分母、配置位置、软件版本、理由、日期与运行证据。**当前Jinlong值全部待决定，运行值全部未运行。** 不能直接把建议列复制成已执行值。涉及临床定义、开发保留角色或研究问题的变更与研究者讨论；常规工具/参数由Jinlong选择。

## 先解释75%与几个容易混淆的“覆盖”

pilot的75%：分母为不含宿主侧翼的病毒参考全长，分子为通过identity、query比对范围、MAPQ和baseQ及归属检查后，深度至少1的位置数。PE重叠不重复计覆盖，同供者readsets用覆盖并集。至少1归属fragment，没有额外平均深度门槛。95% reads identity、90% aligned query、vOTU 85%较短序列AF、CheckV完整度各是不同指标。

降低到30%得到的是“可信部分覆盖”敏感性端点；不能说病毒准确率变成30%。Jinlong可采用其他规则，依据自己的成熟方法及测量经验，在看正式疾病结果前决定并写清分母。无需把pilot数字当成必须沿用的标准。

## D01｜同一供者的多次采样、多文库和多个run，哪些属于同一标本，哪些应该合并？

- 为什么要决定：错误合并或将run当人会改变独立样本数。
- pilot已用：10人/40 readsets/50 clean FASTQ；供者层测量。
- 建议起点：donor/specimen/library/run逐级映射；技术run合并；重复采样预定义。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`identity.merge_policy`、生效run/commit。

## D02｜哪些算CRC和合格对照，腺瘤、0期、术后及有冲突标签者怎么处理？

- 为什么要决定：改变临床含义需研究者共同决定。
- pilot已用：Feng5 CRC/5研究对照；不是全部健康。
- 建议起点：按原文临床定义核查；冲突单列；0期做明确主/敏感性决策。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`clinical.inclusion`、生效run/commit。

## D03｜怎样检查文件完整，失败重试、版本不可用与备用来源怎样记？

- 为什么要决定：只检查路径存在不能验收；失败不是生物学排除。
- pilot已用：既有下载按大小/MD5及收据验收。
- 建议起点：官方大小/MD5＋下载后SHA256；原accession版本锁定。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`download.validation`、生效run/commit。

## D04｜哪些genome可进入正式参考；公开质量缺失和旧CheckM警告怎样处理？

- 为什么要决定：工具/库可替换；残缺/污染会改变prophage检出。
- pilot已用：旧258历史QC与新公开候选元数据分开。
- 建议起点：统一序列QC；常规≥90%完整/≤5%污染；结构锚点警告重查。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reference.qc`、生效run/commit。

## D05｜采用什么分类版本，如何处理旧名、GCA/GCF、同株多组装与近克隆？

- 为什么要决定：最终种数/株数由此确定；不默认配对FASTA相同。
- pilot已用：旧4分析标签；新候选仅显式配对/BioSample身份组。
- 建议起点：固定taxonomy/ANI；别名核查；近克隆避免按95%ANI误删种内变异。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reference.taxonomy_dedup`、生效run/commit。

## D06｜是否需要补参考宿主/谱系缺口，何时锁定最终accession和数量？

- 为什么要决定：禁止按开发P或保留效应追加有利参考。
- pilot已用：pilot258；本次新候选是预选而非最终。
- 建议起点：Part1按预设缺口补齐；联合目录构建前冻结正式参考。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reference.freeze`、生效run/commit。

## D07｜用什么质控工具、接头、质量和长度标准，PE/SE及singleton如何保留？

- 为什么要决定：更严格长度会影响低丰度/小基因组检出。
- pilot已用：fastp0.24.0；min length50、qualified Phred20、unqualified40%；不把默认值当逐批实际资源。
- 建议起点：按成熟pipeline选定并导出完整命令和read数变化。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reads.qc`、生效run/commit。

## D08｜用哪个人参考版本和比对规则，怎样记录去除量？

- 为什么要决定：更广人参考或更严格规则改变微生物可用reads。
- pilot已用：GRCh38-GENCODE-v47；Bowtie2 very-sensitive-local；identity.95/query.90/MAPQ0/-k10；任一mate可信则去对。
- 建议起点：固定参考/哈希、工具版本和参数；保留汇总与清洁reads身份。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reads.dehumanization`、生效run/commit。

## D09｜默认逐样本还是共组装，最短contig、k-mer、线程和内存如何设？

- 为什么要决定：共组装需追踪个体来源与混合株；影响BAM/复制活动。
- pilot已用：MEGAHIT1.2.9逐样本，min-contig-len500；曾分批调整资源。
- 建议起点：优先逐样本保留患者上下文；资源按服务器实测。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`assembly.strategy`、生效run/commit。

## D10｜是否复用pilot输入或中间结果，如何记录兼容性？

- 为什么要决定：无需新先导；重算/复用都不伪增患者。
- pilot已用：10人和258已完成；新服务器不要求复现。
- 建议起点：可完全重算或按输入/软件/数据库/参数哈希选择复用。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reuse.policy`、生效run/commit。

## D11｜用哪个prophage caller或组合，冲突边界如何处理？

- 为什么要决定：全caller取交集可能损失敏感性；单caller需证据分层。
- pilot已用：geNomad/CheckV，保存原contig及区段。
- 建议起点：成熟主caller＋重点结构复核；保留所有边界来源。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`viral.calling`、生效run/commit。

## D12｜病毒最低长度、CheckV完整度/污染、低或未知质量如何纳入？

- 为什么要决定：短片段可误匹配；质量未知不等于假病毒。
- pilot已用：原定量473个vOTU；长度起点5kb。
- 建议起点：明确主/补充/结构资格；保留小片段供解释。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`viral.eligibility`、生效run/commit。

## D13｜患者中什么证据可支持整合及具体宿主，怎样记录多宿主？

- 为什么要决定：参考匹配只证明参考关联，不能直接证明患者宿主。
- pilot已用：原位坐标/侧翼marker/跨界reads与参考连接分别列。
- 建议起点：原contig细菌上下文＋分类优先；参考、CRISPR、MAG分别标来源。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`host.evidence`、生效run/commit。

## D14｜哪些对象确实需要分箱或公共MAG辅助？

- 为什么要决定：全量MAG成本更高；错误bin会误连宿主。
- pilot已用：prophage可在未分箱原contig上发现。
- 建议起点：按宿主证据缺口定向分箱；不设全量MAG前置。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`host.mag_policy`、生效run/commit。

## D15｜怎样分别定义近全长与局部匹配，双方覆盖分母是什么？

- 为什么要决定：短片段高identity不能冒称完整prophage同一。
- pilot已用：ANI≥95%；双方≥85%为广泛，短序列≥85%可为片段。
- 建议起点：分别报告G和M覆盖与方向；保留多对多关系。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`matching.G_M`、生效run/commit。

## D16｜采用什么聚类算法、ANI/覆盖和代表选择，怎么控制片段桥接？

- 为什么要决定：改变代表/测量单元会触发重新建索引和定量。
- pilot已用：95%ANI＋85%较短序列覆盖。
- 建议起点：固定算法/成员/代表/质量顺序；检查短片段连接。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`catalog.clustering`、生效run/commit。

## D17｜病毒与宿主竞争背景怎么建，identity、query覆盖、MAPQ和baseQ设多少？

- 为什么要决定：过滤更严可能丢共享真实reads，更松易交叉比对。
- pilot已用：Bowtie2 end-to-end very-sensitive；.95/.90/MAPQ20/baseQ20/-k20。
- 建议起点：沿用或采用成熟pipeline；宿主prophage区遮蔽；参数逐项记录。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/input_manifests/pilot10_measurement_profile.json)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`mapping.filters`、生效run/commit。

## D18｜多重比对、共享区、成对reads和重复reads怎样计数，什么情况写NA？

- 为什么要决定：不同分配算法改变丰度与检出，不能混矩阵。
- pilot已用：alternative locations审计；PE重叠碱基并集；duplicates保留。
- 建议起点：明确unique/ambiguous/not_evaluable/zero；保持片段与read单位。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/input_manifests/pilot10_measurement_profile.json)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`mapping.assignment`、生效run/commit。

## D19｜一条病毒至少多少比例的位置被合格reads覆盖才算检出，分母和最低深度是什么？

- 为什么要决定：不是75%准确率、75×深度、CheckV完整度或vOTU85%AF。
- pilot已用：病毒全长≥75%位置达到≥1×；至少1归属fragment；均深门槛0。
- 建议起点：75%可作起点，由Jinlong决定；明确不含宿主侧翼。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/input_manifests/pilot10_detection_profile.json)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`detection.breadth_full_min`、生效run/commit。

## D20｜是否另报30%等部分覆盖，如何与主检出分开？

- 为什么要决定：较低覆盖提高敏感性也更依赖局部特异性。
- pilot已用：0.75:457人×vOTU；0.30:695，增加238测量单元。
- 建议起点：主终点固定；部分覆盖作为独立敏感性，不用显著性挑阈值。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/input_manifests/pilot10_detection_profile.json)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`detection.partial_breadth_min`、生效run/commit。

## D21｜相似病毒难以区分时，保留NA、合并测量单元还是特异区段？

- 为什么要决定：合并可提高可测性但损失细粒度；需冻结后统一重算。
- pilot已用：歧义主要来自病毒间共享；11目标有限技术校准。
- 建议起点：依据可归属性选择；测量单元和目录vOTU分开标记。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`mapping.resolution`、生效run/commit。

## D22｜用哪种计数/深度/归一化单位，和检出终点有何区别？

- 为什么要决定：组成效应、测序深度和低覆盖影响不同。
- pilot已用：fragment_rpkm_all_clean；低于检出门槛的非零值未强制置零。
- 建议起点：选一种主丰度定义并保留有效分母；检出单独报告。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`abundance.normalization`、生效run/commit。

## D23｜用什么工具和数据库，标签缺失或低丰度如何记录？

- 为什么要决定：数据库分辨率可能不足以区分参考谱系。
- pilot已用：MetaPhlAn4.2.4/vJan25；B.hominis精确标签缺失为NA。
- 建议起点：固定版本/分类crosswalk；标签不支持与未报告信号分开。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`host.profile`、生效run/commit。

## D24｜哪些宿主—病毒可进入调整，怎样保证前后样本可比？

- 为什么要决定：简单病毒/宿主比值会放大低宿主噪声；调整不等于因果。
- pilot已用：已有10人profile与背景表，未完成可靠疾病宿主调整。
- 建议起点：可信连接＋共同支持/事件足够；同供者子集对照调整前后。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`statistics.host_adjustment`、生效run/commit。

## D25｜数据库版本与同源阈值是什么，如何处理未知和成员差异？

- 为什么要决定：不能把代表功能广播全部成员；同源不等于实验功能。
- pilot已用：PHROGs v4/38880；pyhmmer E/domE≤1e−5、model coverage≥.4。
- 建议起点：按实际成员序列注释；保留best与多命中来源。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/src_pilot6_measurement/bin/core_functions.py)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`functions.phrogs`、生效run/commit。

## D26｜怎样判完整系统与组件，如何区分病毒内部/边界/宿主？

- 为什么要决定：边界误差可把宿主系统算进prophage。
- pilot已用：DefenseFinder模型3.1.0；系统与gene坐标保留。
- 建议起点：选成熟模型并固定版本；关键系统复核完整性和位置。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/src_pilot6_measurement/bin/core_functions.py)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`functions.defense`、生效run/commit。

## D27｜采用哪些反防御模型，怎样表述和复核候选？

- 为什么要决定：不能从命中直接推断感染谱或已证明抑制机制。
- pilot已用：DefenseFinder anti-defense模型命中。
- 建议起点：模型同源候选分级，关键对象结合结构/系统背景。
- 核对入口：[历史方法/配置](../05_existing_results/method_sources/src_pilot6_measurement/bin/core_functions.py)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`functions.anti_defense`、生效run/commit。

## D28｜哪些只称代谢同源，哪些可作为重点AMG/cargo候选？

- 为什么要决定：不必全量VFDB/ARG；扩展须有明确科学问题。
- pilot已用：KOfam/KO同源已完成；未把所有KO当AMG。
- 建议起点：病毒上下文、污染/边界、完整蛋白和功能证据复核。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`functions.amg_review`、生效run/commit。

## D29｜哪些宿主值得种内树和carriage比较，如何控制研究/体位偏倚？

- 为什么要决定：参考株比例不是患者人群比例；不强做病例分离株比较。
- pilot已用：旧参考少数宿主；疾病来源比较资格不足。
- 建议起点：按覆盖/质量做种内核心树；控制相关性与组装可测性。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`reference.phylogeny`、生效run/commit。

## D30｜使用PropagAtE还是其他方法，BAM/identity/比值/效应量怎么设？

- 为什么要决定：改BAM或分母可能改变候选；后续实验支持具体活动。
- pilot已用：PropagAtE1.1.0；identity.97、ratio2、Cohen d.7。
- 建议起点：成熟方法＋原组装上下文BAM；保留工具默认与额外规则。
- 核对入口：[历史方法/配置](../05_existing_results/postpilot_diagnostics/02_pipeline/replication_activity.py)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`replication.method`、生效run/commit。

## D31｜宿主侧翼、边界、深度、重复区和短候选怎样筛？

- 为什么要决定：382可评价不等于382都合格进入主目录。
- pilot已用：病毒/宿主均深≥1、广度≥.5；宿主≥1kb；endmask150；局部±5kb。
- 建议起点：明确结构/覆盖可评价与病毒身份两层资格。
- 核对入口：[历史方法/配置](../05_existing_results/postpilot_diagnostics/02_pipeline/replication_activity.py)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`replication.eligibility`、生效run/commit。

## D32｜哪些对象信息足够进入正式检验，过滤能否与P值无关？

- 为什么要决定：固定前置规则；不是通用法定门槛或阳性保证。
- pilot已用：397/473仅0或1人广泛检出。
- 建议起点：可参考5%检出、20阳/20阴、至少2队列有信息；按规模定。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`statistics.information_filter`、生效run/commit。

## D33｜检出、丰度和功能各用什么模型，稀疏/完全分离怎么办？

- 为什么要决定：模型复杂度受事件数制约；不可估计保留原因。
- pilot已用：pilot Fisher＋精确秩置换；未做多变量模型。
- 建议起点：主检出回归/稀疏合适后端；丰度补充；功能只分析可解释单位。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`statistics.models`、生效run/commit。

## D34｜调整哪些变量，如何处理队列异质性与Thomas同研究来源？

- 为什么要决定：缺失协变量与共享批次需记录；不能简单混样本忽略研究。
- pilot已用：单队列pilot不支持跨队列效应。
- 建议起点：年龄/性别/深度等预定；队列内建模、RE汇总与研究级敏感性。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`statistics.covariates_meta`、生效run/commit。

## D35｜主/补充对象与终点如何分家族，报告哪一套P/Q？

- 为什么要决定：不同家族FDR不能冒称全目录统一控制。
- pilot已用：815项联合BH均Q1；endpoint Q另列。
- 建议起点：F1含M成员主、G-only补；检出/丰度主次事先定义。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`statistics.multiplicity`、生效run/commit。

## D36｜保留队列前锁定哪些目录、候选、方向和校正规则？

- 为什么要决定：保留序列不得完善开发目录；效应看过不能仍称完全独立。
- pilot已用：pilot发现和测量同10人；无独立验证。
- 建议起点：固定两保留集；候选后验证，可考虑Holm FWER.05。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`validation.freeze`、生效run/commit。

## D37｜如何兼顾重复效应、宿主结构、功能与菌株可获得性？

- 为什么要决定：优先关键实验；完整机制作为后续连续课题。
- pilot已用：已有两个复制候选均为control；未最终选实验株。
- 建议起点：选择1–2个问题清楚且可培养对象；临床方向不单一限定。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`candidate.prioritization`、生效run/commit。

## D38｜以后要重画或补注释，哪些文件必须留下，哪些临时输出可删？

- 为什么要决定：删除HMM原始明细后重查需重算；不能说有备份。
- pilot已用：已清理部分DefenseFinder HMM临时目录，系统/基因结果保留。
- 建议起点：序列/坐标/注释/配置/BAM索引和图底表保留；大文件索引。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`retention.policy`、生效run/commit。

## D39｜研究者怎样从图找到数据、代码、命令、参数修改入口和重跑范围？

- 为什么要决定：无需人工复写每条CLI；未执行模板不能冒充实际记录。
- pilot已用：已打包pilot代码/结果/运行配置和哈希。
- 建议起点：workflow自动导出＋小文件Git；记录figure→table→script→run。
- 核对入口：[历史方法/配置](../05_existing_results/CRC-PHIRE_pilot方法结果与正式研究启示.md)。
- Jinlong填写：选择、理由、工具/数据库版本、数值与单位、实际配置路径及键`provenance.capture`、生效run/commit。
