# 已有结果使用说明

先读[综合报告](CRC-PHIRE_pilot方法结果与正式研究启示.md)，了解258份参考、Feng十人及后续诊断的已完成结果和11张图；查图的选择规则、NA／零值、统计和复制活动参数时，读[图注与测量统计方法](图注与测量统计方法.md)。

Jinlong直接推进正式研究，可使用自己的pipeline，自行选择并记录参数；没有新增pilot或复现旧数值的要求。正式研究以[全局研究计划](../00_project/CRC-PHIRE_全局研究计划.md)为准。

| 目录／文件 | 用途 |
|---|---|
| `reference258/` | 参考基因组、预测prophage与质量、坐标结果 |
| `pilot10/` | 十人联合目录连接、供者测量、探索统计、功能／边界证据及七张图 |
| `postpilot_diagnostics/` | 宿主profile、歧义和技术校准、原位点复制活动、实际代码与四张图 |
| [method_sources/README.md](method_sources/README.md) | 输入与运行记录；可读主源码、版本差异与哈希映射，可恢复三个原版本 |
| [pilot_verified_metrics.tsv](pilot_verified_metrics.tsv) | 综合报告数值与底表的对账入口 |
| [source_data_index.tsv](source_data_index.tsv) | 九份相同图源数据的原路径、唯一底表及SHA256映射 |

图注直接链接科学底表；各图`source_data/`保留独立显示表、别名和设置。数据量单位、缺失原因与分母见方法文档，不把实例、序列、vOTU和患者数混用。
