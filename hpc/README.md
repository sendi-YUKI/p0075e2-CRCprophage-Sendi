# hpc/ — spark 与 HMGU HPC 之间的大文件中转

本目录在 spark 与 HMGU HPC 之间搬运任务文件:代码与小文件走 git,大文件走
`gcln`(Google Drive `My Drive/transfer`)。两侧的人与 AI 都按本文件操作。

## git 追踪什么

`.gitignore` 默认忽略 `hpc/` 下的全部内容,只放行三类:

- `hpc/README.md`(本文件)
- `hpc/<bundle>/README.md` 与 `hpc/rawdata/README.md`
- `hpc/<bundle>/scripts/` 下的全部内容:sbatch 脚本、samplesheet、小配置

其余文件只经 `gcln` 传输。拿不准某个文件归哪边,用 `git check-ignore -v <路径>` 查。

## 布局

一项 HPC 任务一个 bundle,原始数据单放 `rawdata/`:

```
hpc/
├── rawdata/                 原始数据,一个交付批次一个目录(按需建立)
│   ├── README.md            git:各批次的内容、来源与时限
│   └── d0xxx/               批次号 = transfer 下的名字
└── p0075e2s-<NNN>-<slug>/     bundle 名 = transfer 下的名字
    ├── README.md            git:目的、输入来源、命令、状态
    ├── scripts/             git:sbatch 脚本、samplesheet、小配置
    ├── in/                  gcln:spark → HPC 的输入
    └── out/                 gcln:HPC → spark 的结果
```

- bundle 名以项目代号 `p0075e2s-` 开头。`transfer` 是所有项目共用的一层扁平目录,
  同名 bundle 会被 `gcln` 合并,前缀让各项目的数据互不串扰。
- `<NNN>` 取消费结果的分析步骤编号,即该步骤在 `analyses/data/` 下的目录编号。
- `rawdata/` 下的目录名是 labdata / colabdata 登记的批次号(`d0xxx` / `d9xxx`),
  批次号本身全局唯一,直接作传输名。批次内保持厂商交付的原始结构与文件名,只读;
  多个 bundle 共用同一批原始数据。

## 往返

上传整个 bundle,下载只取对向的那一半。`gcln` 按 checksum 跳过已有内容,
所以重复上传很便宜;`scripts/` 两侧都以 git 为准。

spark → HPC:

```bash
# spark
cd <repo>/hpc && gcln up p0075e2s-<NNN>-<slug>
git push                                   # README 与 scripts/
# HPC
git pull
cd <repo>/hpc/p0075e2s-<NNN>-<slug> && gcln dl p0075e2s-<NNN>-<slug>/in
```

HPC → spark:

```bash
# HPC(计算节点上也能跑,可写进 sbatch 末尾)
cd <repo>/hpc && gcln up p0075e2s-<NNN>-<slug>
git push                                   # 更新后的 bundle README
# spark
git pull
cd <repo>/hpc/p0075e2s-<NNN>-<slug> && gcln dl p0075e2s-<NNN>-<slug>/out
```

原始数据从 NAS 归档区直接上传,HPC 侧下到 `rawdata/`:

```bash
# spark
cd /mnt/nasdata/hpc-hmgu/data_archive && gcln up d0xxx
# HPC
cd <repo>/hpc/rawdata && gcln dl d0xxx
```

`gcln` 的行为(实测见 [vpipe RUNBOOK](file:///home/allen/vpipe/docs/RUNBOOK.md)「gcln」一节):

- 相同内容跳过,内容变化的同名文件原位替换,目标端多出的文件保留。
  本地删掉的文件在 Drive 上仍在,用 `gcln rm <名字>` 撤回(只作用于 transfer,移入 Drive 回收站)。
- 每次调用固定开销约 5.5 s,每个小文件上传约 0.45 s。上千个小文件先打成 tar 再放进 `in/` 或 `out/`。
- Google Drive 每个账号每天最多上传 750 GB,超过这个量的 bundle 分天传。

## 规矩

- **`in/` 由 spark 侧分析步骤产出。** 从 `analyses/data/<step>/` 复制或打包进来,
  bundle README 写明来源步骤与当时的 `git rev-parse --short HEAD`。
- **`out/` 只放 spark 要用的结果。** 中间产物留在 HPC 的计算目录。
  HPC 侧的仓库克隆与计算目录放在 `/ictstr01` 项目区,大文件不进 `$HOME`。
- **分析步骤只读 `analyses/data/`。** 结果回到 spark 后,经用户确认落到
  `analyses/data/00-raw/hpc/<bundle>/`(配 README 与 `RUNLOG.tsv`),下游步骤从这里读取,
  并按[代码参数与结果归档规范](../03_execution/代码参数与结果归档规范.md)第 7 节登记进大文件索引。
  落位细则见 [HPC 计算产物回迁](file:///home/allen/.configure/rules/refs/hpc-repatriation.md)。
- **bundle README 记运行状态。** 谁推进了一步,谁更新状态并 commit:目的、pm 任务 ID、
  输入来源、HPC 计算目录、提交命令与 job ID、当前状态。
- **收尾由用户拍板。** 结果落位并核对后,`gcln rm <bundle>` 撤掉 Drive 上的 bundle,再清两侧的 `in/`、`out/`。
- **spark 侧的备份。** `~/github/rujinlong` 整体进每日 NAS 备份,`hpc/` 的内容一并备份。
  原始数据的原件已在 NAS 归档区时,spark 侧不再往 `rawdata/` 复制。
- **HMGU 的环境事实**(分区、vpipe、conda、资源档位、踩过的坑)写进
  [vpipe RUNBOOK](file:///home/allen/vpipe/docs/RUNBOOK.md) 的 `env: hpc` 节,各项目共用。

## 调度与动态迭代

- 目标是在账号权限、有效站点要求与授权预算内,缩短排队、计算及传输总时间,提高成功任务吞吐量。
- 默认显式 `--nice=0`,partition 与 QoS 按实时资源和任务需求选择;正数 nice 只在用户明确指定或有效站点强制规定要求时使用。
  核对 wrapper、模板及 Nextflow 最终生成的提交参数。
- walltime 按预计运行时间加合理余量,CPU、内存与 GPU 按程序需求申请;并发、分块和小任务合并按容量与预算调整。
  站点快照中的资源数字和示例参数在使用前核对。终端迭代、notebook 用 interactive;独立短测试可用站点允许的 batch 队列。
- 排队接近预计运行时间时查 PENDING 原因与有效请求;在已有授权内按证据调整队列或资源,重投前核对原 job 状态。
  preemptible 按中断概率、重跑成本与恢复能力评估。
- **HPC AI 已获授权根据实际 Slurm jobs 动态迭代调度策略。** 在现有 REPORT(使用 bundle README 或 HANDOFF 的项目写状态节)记录
  job ID、输入规模、硬件、有效参数、Submit/Eligible/Start/End、CPU/GPU 利用、峰值内存、I/O/传输与失败重跑证据;
  缺失指标记 unknown。策略注明适用范围,单次观察先作为候选,后续真实任务或已授权的小测试验证;比较总完成时间、
  成功吞吐量与重跑成本,出现退化时恢复已验证参数。当前验收完成即交付。
- 任务特有参数留 bundle;可复用规则更新本文件或当前可写的共享 hpc skill,commit + push 并在报告引用 commit。
  spark 接收后回写 canonical skill 与模板、同步副本并录入 pm;私有站点材料留在内部。
  迭代沿用用户优先级偏好、有效站点要求、账号权限与预算;扩大授权或改动受保护数据时由用户决定。

## bundle 索引

每建一个 bundle 加一行;收尾后删去该行,记录留在 git 历史与 pm。

| bundle | 用途 | pm 任务 | 状态 |
|---|---|---|---|
| [p0075e2s-010-refgenome](p0075e2s-010-refgenome/README.md) | 1,931 份参考组装的下载、QC、分类与 prophage 调用(G) | 决策表 D03–D05、D11、D12 | 待运行 |
| [p0075e2s-020-cohort](p0075e2s-020-cohort/README.md) | 九队列身份表(A 段),冻结后的患者上游(B 段,M) | 决策表 D01–D03、D07–D12 | A 段待运行 |
