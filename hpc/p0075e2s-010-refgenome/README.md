# p0075e2s-010-refgenome

| 项 | 值 |
|---|---|
| 目的 | 1,931 份预选细菌参考组装:下载与校验、序列 QC、分类与同株核查、prophage 调用,产出参考端 prophage 实例 G 的证据表 |
| pm 任务 | 本仓库不用 pm;方法决定对应 [决策表](../../03_execution/Jinlong_方法与参数决策表.md) D03、D04、D05、D11、D12 |
| 输入来源 | [download_accessions.txt](../../01_reference/expanded/download_accessions.txt)、[preselected.tsv](../../01_reference/expanded/preselected.tsv) @ `7c901e8`,随 git 到达,本 bundle 不用 `in/` |
| HPC 计算目录 | 建议 `/ictstr01/project/genomics/jru/p0075e2s/010-refgenome/`,开工时写入实际路径 |
| 提交命令 | |
| job ID | |
| 状态 | 待运行 |

## 开工前读

1. [参考集合选择说明](../../01_reference/expanded/参考集合选择说明.md) 第 4、6 节:QC 起点、条件预选的 81 份、固定版本下载。
2. vpipe 仓库 `docs/RUNBOOK.md` 的「env: hpc」节;`bin/p0075e2.slurm --help` 打印本步沿用的参考基因组链(下载、CheckM2、GTDB-Tk、skani、合并后的 prophage 调用、PHROGs)及各步的坑。

## 做什么

本步对全部 1,931 份组装计算证据,去留在 spark 侧冻结正式参考时决定。95% ANI 去冗余会删掉种内 prophage 变异,所以本步不跑 dRep,只标记同株别名与近克隆。

1. **下载与校验**(D03)。用 `datasets` 按 `download_accessions.txt` 的 `accession.version` 下载(`download.slurm assembly <list> <outdir>`,或参考说明第 6 节的 dehydrated + rehydrate)。
   - 写 `download_status.tsv`:`accession`、`status`、`bytes`、`official_md5`(取自数据包的 `md5sum.txt`)、`sha256`、`failure_reason`。
   - 指定版本取不到时记为失败,不换成 latest。下载失败与 QC 排除分开记录。
   - 基因组文件统一命名 `<accession>.fa` 放进 `genomes/`(checkm2 写死 `.fa`)。
2. **序列 QC 与分类**(D04、D05)。
   - `sbatch -c 24 --mem 120G assembly.slurm checkm2 p0075e2s_ref genomes .fa`
   - `p0075e2.slurm n1_gtdbtk genomes .fa`(加 `q_` 前缀规避 GTDB 参考 ID 撞名)
   - `seqkit stats -a -T genomes/*.fa` 取 contig 数与 N50
   - `sbatch -c 32 --mem 64G p0075e2.slurm n1_clonal_clusters genomes .fa 99.9`,输出 skani 边表与近克隆分组
3. **prophage 调用**(D11、D12)。沿用 N1 链:
   - `p0075e2.slurm n1_merge_genomes genomes .fa` 合并成一个 FASTA,并生成 `genome_id_map.tsv`
   - 先用 3 个基因组跑通、核对 header 保留 `gNNNN_`,再全量 `nextflow run <vpipe>/main.nf --fasta merged_genomes.fa --analysis prophage -profile slurm -resume`(PhiSpy + PhageBoost + geNomad,保留各工具的 BED 与 consensus)
   - 对 prophage 区段跑 CheckV(`--analysis checkv_sg` 或 `anno_contig.slurm checkv`)
   - `p0075e2.slurm n1_prophage_proteins` → `anno_prot.slurm phrogs` → `n1_merge_phrogs`
4. **汇总进 `out/`**,上传前核对文件清单与下面「结果」一致;更新本 README 状态并 push。

坐标一律写成 0-based 半开区间,并在表头注明。实际命令、sbatch 脚本与配置放进 `scripts/`,随 git 提交。

## 停下来问用户

- 下载失败超过 2%,或集中在某个宿主。
- HMGU 缺少所需容器或数据库,需要安装新软件或改 vpipe 本体。
- 重点宿主(Fusobacterium、P. micra、P. anaerobius、B. fragilis)的 GTDB 分类与 NCBI 标签大面积冲突,或 CheckM2 污染普遍偏高。

其他工具参数、资源与可重试的失败按 RUNBOOK 处理,写进 RUNLOG。

## 结果

`out/` 只放 spark 要用的表与序列;全基因组 FASTA、Bakta 全量注释和中间产物留在 HPC 计算目录,在 `large_file_index.tsv` 登记。

| 文件 | 内容 | 生成命令 |
|---|---|---|
| `download_status.tsv` | 逐 accession 下载状态、字节、MD5、SHA256 | |
| `assembly_stats.tsv` | seqkit 统计(contig 数、N50、总长) | |
| `checkm2_quality_report.tsv` | CheckM2 完整度与污染 | |
| `gtdbtk.bac120.summary.tsv` | GTDB-Tk 分类(含 GTDB release) | |
| `skani_edges.tsv`、`clonal_clusters.tsv` | 两两 ANI/AF 与近克隆分组 | |
| `genome_id_map.tsv` | `gNNNN` 与 accession 的对应 | |
| `prophage_regions.tsv` | 各 caller 与 consensus 的区段、坐标、得分 | |
| `prophages.fna.gz`、`prophage_flanks.fna.gz` | prophage 区段与两侧各 10 kb 侧翼 | |
| `checkv_quality_summary.tsv` | CheckV 质量 | |
| `prophage_phrogs.tsv` | 逐蛋白 PHROGs 注释 | |
| `software_versions.tsv` | 工具、容器、数据库及版本 | |
| `large_file_index.tsv` | 留在 HPC 的大文件路径、字节、SHA256 | |
| `RUNLOG.tsv` | `started`、`jobid`、`workdir`、`command`、`state` | |
