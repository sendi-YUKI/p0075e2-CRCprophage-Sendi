> 历史诊断阶段方案，仅作结果解释；当前执行以[新全局计划](../../../00_project/CRC-PHIRE_全局研究计划.md)为准。Jinlong不再做新先导。

# Pipeline调整清单与下一批操作

本轮在现有Nextflow＋Slurm流程上补诊断模块。原258参考、十人QC/组装、病毒目录和定量结果继续复用。下表区分真实完成的改动与正式扩量前仍需接通的部分；“软件已经安装”不等于该研究接口已经验收。

## 1. 本轮已做的必要改动

| 模块 | 复用内容 | 本轮调整与输出 | 验收与兼容方式 |
|---|---|---|---|
| 输入身份 | `pilot10_clean_reads.json`、`pilot10_frozen_assemblies.json`及下载/归档收据 | 10名供者、40 readsets、50 FASTQ，逐人profile输入与原组装clean输入完全一致 | `donor_input_identity.tsv`；大型输入按既有SHA重新核验一次，避免每个阶段重复扫描 |
| 宿主测量 | 现成MetaPhlAn 4.2.4、vJan25数据库、`cohort_tools.py`解析器 | 数据库清单格式适配；按供者合并不重复的SE/PE/孤儿reads；逐人profile、mapout、命令和耗时 | 原数据库不改；1132预检失败保留，1133为恢复作业；最终以十份收据为准 |
| 测量诊断 | 原4730行测量及40份歧义表/深度 | 按fragment区分host-only、病毒间、病毒—host及同目标多位置/饱和；0.75/0.30完整状态重算 | 0.75重算逐行重现原状态；不重新比对患者，不把展开计数当唯一reads |
| 有限技术校准 | 原竞争参考、mask、Bowtie2索引、原归属函数 | 技术类别选11个不同目标（12个类别实例）；SE/PE、已知阳性与近缘/宿主背景混合；并恢复同一组对象真实覆盖曲线 | 模拟独立标识；不按P值选对象；保存全部命令与source data；只支持所测技术范围 |
| 复制活动 | 十人原组装FASTA、全组装readback BAM、原位点坐标 | 抽取所需完整contig及其原始比对；PropagAtE固定版本；输出382个位点状态、覆盖曲线和有限敏感性 | 序列切片SHA、坐标长度及独立覆盖均值与工具结果比较；官方输出和项目可评价状态分别保留 |
| 统计复算 | 原探索性统计表 | 重算BH、5＋5下界反例和原功效表 | 原815个检验结果不修改；纠错表说明后果 |
| 下一批清单 | 九队列inventory、Feng107候选和现有下载收据 | 逐run记录已核验字节、剩余字节、已完成十人、待临床冻结状态 | 尚未提交扩量；不把候选人数写成最终合格人数 |
| 输出与可恢复性 | 原路径约定 | 所有新小文件在`14_postpilot_adjustment`；大工作文件在NVMe `work/postpilot_adjustment` | 原报告/收据不覆盖，有限目录同步有SHA回执，真实图需视觉检查 |

本轮宿主计算已经完成：1133与1139分别实际处理5名供者，1139的末份结果由1141原子移交，1133对已完成结果校验后跳过，未重复映射同一供者。1140完成10份profile及输入身份验收、382个位点BAM坐标核查和宿主图；1142生成诊断报告与方案HTML。实际代码见`../02_pipeline/`，收据见`../09_handoff/host_convergence_acceptance.json`和各模块`completed.json`。这些完成状态不代表下表的正式扩量接口已经全部实现。

## 2. 正式扩量前的最小调整

| 优先级 | 待做内容 | 具体输入/输出 | 如何验收 |
|---|---|---|---|
| P0 | 批量clinical/subject/sample/run清单冻结 | 原作者/策展/ENA清单→纳入排除表、去重表、最终runs、SHA冻结JSON | 样本/分组冲突逐项有结论；同供者唯一；开发/保留用途一致 |
| P0 | 上游driver参数化 | 替换原`run_remaining14.py`固定14人、16核/40GiB、全量归档等硬编码；参数指定新manifest/outdir/resources | 复用1个完成样本不重跑；1个新开发样本从下载到病毒预测完成；Nextflow与Slurm资源一致 |
| P0 | 存储与发布更新 | NAS输入先本地暂存；完成的小结果独立发布，大工具输出保留单份并记录依赖 | 删除某个明确可再生临时文件前，消费者均不依赖它；不对旧raw软链接直接清空work |
| P0 | 纳入/测量家族与版本 | C及成员来源→M_discovery-containing、G-only、实例到测量ID交叉表 | G∩M只计一次；保留队列序列不混入C；更换目录只重算受影响阶段 |
| P0 | 正式二元结局和meta后端 | 冻结检测表、临床表→队列log OR/SE/CI、meta、FDR、失败原因 | 小型已知分离案例和有限真实对象验证；稀疏/不可估计不伪造有限效应；报告完整分母 |
| P0 | 宿主子集接口 | phage-host证据+分类crosswalk+MetaPhlAn→可测子集、同样本调整前后比较 | 物种/属/SGB分辨率一致；多宿主和无共同支持明确留出；零与NA不混淆 |
| P1 | 功能成员级比较 | 完整比较集合内一致注释、长度/质量/宿主/供者→预定义功能结果 | 代表注释不无条件传播；碎片缺失不作真阴性 |
| P1 | 复制活动正式接入 | 本轮适配器＋结构资格清单→逐位点结果、稳定子集、队列描述 | 把覆盖可评价与目录/宿主证据分开；新NAS输入先本地暂存；不要求每位患者有活跃位点 |
| P1 | 验证与候选冻结 | 目录/候选/方向/模型版本→两保留队列结果 | 候选冻结早于保留疾病效应；Holm范围完整；不支持/不可评价也输出 |
| 按需 | 关键宿主分类、定向MAG、谱系及补参考 | 实际候选提出的问题→更精确宿主/结构/菌株资源 | 每项有明确对象与停止条件；不启动全量附加分析 |

正式统计的具体信息过滤和稀疏模型后端仍标`conditional`，见`../06_statistics_design/analysis_rules.json`。科学问题、主终点、发现/验证边界已明确；条件项需要的是执行前的有限实现验收，不是重新讨论整个课题。

## 3. 下一批实际入口

现有九队列清单均有候选来源记录，但尚未完成全部供者级临床冻结。Feng候选107人、214 runs、321 FASTQ，合计490,654,560,109字节；其中71,880,930,564字节有既有下载收据并通过本轮大小/修改时间复核，候选池剩余418,773,629,545字节（约418.8 GB）。这是候选池预算，最终下载量随临床纳排减少；不需要重新下载已完成十人的数据。

本轮生成：

- `../01_manifests/next_cohort_inventory.tsv`：七开发、两保留及每个来源待完成事项。
- `../01_manifests/next_feng_candidate_samples.tsv`：107名候选身份、组别、现有进度和临床状态。
- `../01_manifests/next_feng_candidate_runs.tsv`：真实run、URL、官方MD5/字节、已下载/待下载量。
- `../01_manifests/next_batch_inventory.json`：规模与当前`not clinically frozen`事实。

先解决原作者临床来源中的具体冲突，并输出选定开发批次的`final_runs.tsv`和`clinical_freeze.json`。后者至少包含：

```json
{
  "status": "clinically_frozen",
  "batch_id": "development_batch01",
  "run_manifest_sha256": "实际final_runs.tsv的SHA256",
  "eligible_sample_ids": ["实际合格样本ID"],
  "donor_deduplication_checked": true,
  "clinical_sources": ["实际原始临床表及核查记录"]
}
```

在Spark上先只检查并生成适配副本：

```bash
spark-status
A=/srv/CRC-PHIRE/analysis1/14_postpilot_adjustment
python3 "$A/02_pipeline/next_download.py" \
  --clinical-freeze "$A/01_manifests/clinical_freeze.json" \
  --runs "$A/01_manifests/final_runs.tsv"
```

真实清单验收通过且准备推进下一批时，提交：

```bash
sbatch --export=ALL,CLINICAL_FREEZE="$A/01_manifests/clinical_freeze.json",RUN_MANIFEST="$A/01_manifests/final_runs.tsv" \
  "$A/02_pipeline/next_download.sbatch"
```

该入口复用已有`prefetch_engineering16.py`的断点续传与MD5/SHA逻辑，只更换输入manifest和新的状态输出目录，保存父代码/适配代码SHA。当前缺少最终临床冻结，所以**入口已准备，批次状态为not_ready；本轮未提交下载作业**。上游组装driver仍需按上表参数化，不能把原1090脚本直接作为全量入口。

跨队列技术检查优先选两个其他开发来源的少量平衡样本，例如每来源4 CRC＋4对照，按合格身份和深度层选择，避免按是否已有phage命中选择。确定样本后纳入整个开发计算范围，不把这16人变成新的永久样本量上限；保留队列继续保持独立。

## 4. 资源、耗时和效率

本轮原位点复用BAM的实际计算约十余分钟。1133（6核/44GiB）总墙钟2小时7分22秒，1139（14核/64GiB）并行运行57分12秒，两作业时间重叠，不能相加当作用户等待时长。逐人宿主计算耗时和真实线程见`../03_host_profile/host_profile_qc.tsv`；1141移交1秒，1140验收/汇总/宿主图7秒，1142报告4秒。技术校准是数万条已知来源模拟reads和原覆盖数组，远小于重跑真实mapping。

1139退出后，hui的新任务已使用12核/64GiB；轻量移交1141从4GiB调整为1GiB后立即完成。后续1140/1142利用1133释放的资源自动完成，没有修改其他用户作业。用户要求提高1140/1142优先级的请求被Slurm权限拒绝，优先级未改变；实际结果交付未因此等待hui完成。资源调整和权限结果均保留在`../09_handoff/`。

新一批所需墙钟时间由下载、QC/组装、病毒预测和共享排队决定。旧十人组装已经说明每人的上游可达小时级；本轮快速位点小试不能外推为全量宏基因组也能十几分钟完成。正式估算采用：各深度层样本数×实测CPU-hours/可持续CPU，加下载与I/O瓶颈及有限重试。九队列最终人数和总字节未冻结前，提供数周量级规划，具体完成日期暂待清单和跨队列吞吐。

效率优先顺序：复用完成输入/阶段；批量参数化和合理并发；目录冻结后一次正式测量；小结果独立发布、大raw单份保留；只为具体测量问题做局部修复。不追加与先导决策无关的工具、全量MAG或全面重注释。

## 5. 服务器新说明

后续先读[ai.md](http://SERVER_ADDRESS:3001/ai.md)和`spark-status`，使用`/opt/bio/share/templates`，新工具先查共享安装、固定版本，NAS输入先暂存到NVMe。旧参数不能取代实时配额；不调整他人排队或使用sudo。内存保护器冻结时按管理员指引处理。

当前共享工具更完整，可减少后续安装工作；不需要因此重跑已验收的旧结果。正式迁移共享路径先验证版本/内容及小样本兼容性，记录差异。
