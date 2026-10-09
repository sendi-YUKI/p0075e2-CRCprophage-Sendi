# p0075e2s · CRC-PHIRE 正式研究

本仓库是 Sendi 交接的 CRC-PHIRE 研究包,Jinlong(用户)在 `dev-ru` 分支上推进正式计算。研究问题、分工与顺序见 [README](README.md) 和 [Jinlong 执行说明](00_project/CRC-PHIRE_Jinlong执行说明.md)。

## 目录

| 路径 | 内容 |
|---|---|
| `00_project/`–`06_literature/`、`90_package/` | Sendi 的设计、输入清单、模板与试跑结果;作为输入读取,修改前先问用户 |
| `hpc/` | 与 HMGU HPC 交换大文件的中转区,在 HPC 跑任务前先读 [hpc/README.md](hpc/README.md) |
| `analyses/data/` | spark 侧的分析数据与 HPC 回迁结果,git 只追踪 `README.md` |

## 在 HMGU HPC 上

1. `git pull`,读 [hpc/README.md](hpc/README.md),再读「bundle 索引」中未收尾的 bundle README,按其中「做什么」执行。两个 bundle 互不依赖,可并行。
2. 工具调用以 vpipe 为准:vpipe 仓库 `docs/RUNBOOK.md` 的「env: hpc」节,计算一律 `sbatch`。
3. 每推进一步,更新 bundle README 的状态、计算目录、提交命令与 job ID,commit 并 push 到 `dev-ru`。

## 规矩

- 仓库在 GitHub 公开。凭据、token、邮箱与非公开的个人数据不进 git。
- [决策表](03_execution/Jinlong_方法与参数决策表.md) D01–D39 由用户决定。bundle README 写明的建议可以直接执行;改变临床定义、队列角色或主终点的事项停下来问用户。
- holdout 队列(WIRBEL_DE_2019、YANG_CN_2021)不参与开发目录构建、参数挑选与候选排序。
- 大文件只经 `gcln` 传输;git 只放 README 与 bundle 的 `scripts/`。
