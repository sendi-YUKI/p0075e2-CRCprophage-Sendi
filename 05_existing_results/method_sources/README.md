# 方法源码、版本差异与真实运行证据

结果见[综合报告](../CRC-PHIRE_Sendi试跑方法结果与正式研究参考.md)，完整图注、测量和统计说明见[图注与测量统计方法](../图注与测量统计方法.md)。这里保存实际方法的可追溯记录；Jinlong可使用自己的pipeline并记录实际方法与参数，无需再跑pilot或复现旧数值。

| 路径 | 内容 |
|---|---|
| `src_pilot6_measurement/` | 可直接阅读的主源码快照；十人联合目录与测量继承此树，名字pilot6不表示只分析六人 |
| `source_overlays/src_reads_v1/` | reads身份、QC、去人源和组装版本中与主树不同或主树缺少的文件 |
| `source_overlays/src_campaign_v2/` | 参考及早期模块版本中与主树不同或主树缺少的文件 |
| [source_snapshot_manifest.json](source_snapshot_manifest.json) | 三个原版本的全部路径、对应主树／差异文件、字节数和SHA256；相同文件仅存一份 |
| `pilot_helpers/`、`pilot10_runtime_bin/` | 十人统计、证据整合、绘图及实际辅助逻辑 |
| `input_manifests/` | 供者和输入身份、组装来源、测量与检出profile |
| `run_records/` | 逐阶段实际参数、代码、软件／容器、数据库及完成记录 |
| [诊断代码](../postpilot_diagnostics/02_pipeline/) | 宿主、测量校准、复制活动与诊断图代码 |

差异目录需要按manifest恢复后才能组成对应完整版本；不能把主树的`viral_catalog.py`、`cohort_measure.py`等直接当成所有历史版本的同一实现。源码存在某个可选模块也不等于该模块已运行。运行范围由完成收据及结果表确定，各批CPU／内存以实际命令和参数为准。

使用[compact_pilot_sources.py](../../90_package/tools/compact_pilot_sources.py)可校验存储内容，或在空目录中恢复三个完整源码树。从`analysis2/`目录执行：

```bash
python 90_package/tools/compact_pilot_sources.py verify --archive 05_existing_results/method_sources
python 90_package/tools/compact_pilot_sources.py restore --archive 05_existing_results/method_sources --output restored_pilot_sources
```

恢复会核对全部文件SHA256和目录结构，并拒绝覆盖非空输出目录。本地源码按原字节保存；公开源码中的内部路径使用逻辑别名，manifest校验该公开副本的实际字节。两种副本分别维护哈希，不把公开副本哈希称为原服务器源码哈希。重画或移植时需设置数据根、输入及环境；源码记录用于审阅方法，不能直接代表新服务器配置。

软件工具版本与数据库版本分开保存。第三方代码保留原版权与许可证；文献和数据库使用各自许可，不重新授予第三方内容许可。
