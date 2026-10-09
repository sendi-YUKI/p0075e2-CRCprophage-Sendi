# p0075e2s-020-cohort

| 项 | 值 |
|---|---|
| 目的 | 九队列患者端上游:A 段建立供者—标本—文库—run 身份表与合格 run 清单;B 段在身份冻结后下载、QC、去人源、逐样本组装与病毒/prophage 调用,产出患者端实例 M |
| pm 任务 | 本仓库不用 pm;方法决定对应 [决策表](../../03_execution/Jinlong_方法与参数决策表.md) D01–D03、D07–D12 |
| 输入来源 | [九队列说明](../../02_cohorts/九队列正式研究说明.md)、[来源清单](../../02_cohorts/nine_cohort_sources.tsv)、`02_cohorts/sources/`、`02_cohorts/context/` @ `67c7568`,随 git 到达,本 bundle 不用 `in/` |
| HPC 计算目录 | 建议 `/ictstr01/project/genomics/jru/p0075e2s/020-cohort/`,开工时写入实际路径 |
| 提交命令 | |
| job ID | |
| 状态 | A 段待运行 |

本 bundle 分两段。A 段完成后停下,由用户与 Sendi 裁决临床身份;用户把冻结的 run 清单提交为 `scripts/runs_frozen.tsv` 并把状态改成「身份已冻结」后,再开始 B 段。

## A 段:身份表与 run 清单(D01、D02)

先读 [九队列说明](../../02_cohorts/九队列正式研究说明.md) 全文,它列出了每个队列已知的冲突与连接线索。

1. 按 [clinical_final 模板](../../04_deliverables/templates/clinical_final.template.tsv) 的字段,为九个队列建立 `donor_id → specimen_id → library_id → run_accession` 表。
   - 临床分组以原文、补充表与作者资料为准,策展表与 ENA/SRA 用来定位归档身份。
   - 每个供者给出 `clinical_group`(CRC / control / adenoma / other / ambiguous)、`role`(discovery / holdout)与证据出处。
2. 逐条处理九队列说明「待解决事项」表里的对象。证据足以裁决的写明依据;证据不足的标 `ambiguous`,写进 `identity_issues.tsv`(对象、冲突内容、各来源说法、建议处理)。
3. 只保留 shotgun(WGS)run,技术 run 归到同一标本;同一供者的多份标本单列 `duplicate_group`。
4. 由 `02_cohorts/sources/*_ena.tsv` 生成候选 run 清单:`run_accession`、`fastq_https`、`fastq_md5`、`fastq_bytes`、`library_layout`,并按队列汇总人数与下载量。八个 BioProject 的全部 WGS run 合计约 3,857 个、7.4 TB,共享项目的 run 不要重复计入。
5. 结果放 `out/A/`,`gcln up` 后把状态改成「身份表待裁决」并 push,然后停下。

A 段只读写元数据。耗时的 ENA 查询循环用 `sbatch -c 1 --mem 4G` 提交,登录节点只跑秒级查询。

## B 段:患者上游(D03、D07–D12)

只处理 `scripts/runs_frozen.tsv` 中的 run。开发队列走完整流程;holdout 队列(WIRBEL_DE_2019、YANG_CN_2021)只做下载、QC 与去人源,它们后续用冻结目录定量,组装与调用都不跑。

下表是沿用 Sendi 试跑与 vpipe 现成做法的建议,用户冻结身份时在同一次提交里确认或修改。

| 步骤 | 建议 | vpipe 入口 |
|---|---|---|
| 下载与校验 | ENA https 直下,按 `fastq_md5` 校验,失败重试后单列 | `download.slurm ena_fastq <run> <outdir>`,array over runs |
| reads QC(D07) | fastp,PE 保持配对,singleton 单独保存 | `assembly.slurm fastp` |
| 去人源(D08) | 用 HMGU 现有人类参考(`$DB_BBMAP_REF_HUMAN`),identity 0.9;试跑用的是 Bowtie2,两者按用户确认择一 | `assembly.slurm decontam <id> R1 R2 0.9 human` |
| 组装(D09) | 逐样本 MEGAHIT,同一标本的技术 run 先合并 | `assembly.slurm megahit` |
| 病毒与 prophage 调用(D11、D12) | geNomad + CheckV,保留 caller 得分、边界、CheckV 质量与原 contig | `--analysis genomad_sg`、`--analysis checkv_sg` |
| 回贴 BAM | clean reads 回贴本样本 contig,保留 BAM 与索引,供边界复核与复制活动分析 | `abundance.slurm bowtie2_idx` → `abundance.slurm bowtie2`(`-k 1 --very-sensitive`,产出排序 BAM) |

资源按 RUNBOOK「env: hpc」第 2 节:`cpu_p` 不写 `--mem` 时每 CPU 只有 100 MB;一天以上的驱动作业用 `--signal` 自我重投;整批是否完成以产物计数为准。

B 段 `out/` 只放 spark 要用的结果:逐样本的 QC 与去人源读数统计、组装统计、病毒区段表、`viral_contigs.fna.gz`、CheckV 表、`software_versions.tsv`、`large_file_index.tsv` 和 `RUNLOG.tsv`。FASTQ、组装全集与 BAM 留在 HPC 计算目录,登记进 `large_file_index.tsv`。

## 停下来问用户

- A 段结束时(必停)。
- B 段中,某队列下载失败或去人源后读数过低的样本超过 5%。
- 发现同一供者跨队列重复、或临床分组与冻结表矛盾。
- 需要安装新软件、改 vpipe 本体,或预计计算目录超过 15 TB。

## 结果

| 文件 | 内容 | 生成命令 |
|---|---|---|
| `out/A/clinical_identity_draft.tsv` | 九队列身份表草稿,字段同 clinical_final 模板 | |
| `out/A/identity_issues.tsv` | 待裁决对象、证据与建议 | |
| `out/A/candidate_runs.tsv` | 候选 WGS run、URL、MD5、字节 | |
| `out/A/cohort_summary.tsv` | 各队列 CRC / control / 其他人数与下载量 | |
| `out/B/…` | B 段开始后补齐 | |
