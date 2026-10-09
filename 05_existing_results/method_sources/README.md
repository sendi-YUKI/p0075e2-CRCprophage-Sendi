# 历史方法源码和真实运行证据

src_reads_v1：reads身份/QC/去人源/组装改造；src_campaign_v2：参考及早期模块；src_pilot6_measurement：十人合并目录/测量所继承的源码树（目录名pilot6不代表只算六人）。pilot_helpers包含实际十人报告、统计、证据和绘图逻辑；run_records是实际运行配置/数据库/工具manifest；input_manifests给出十人身份与测量profile。postpilot_diagnostics/02_pipeline另存后诊断代码。

这些是历史源码快照，不代表每个文件里的可选模块都已经执行。实际是否运行以run_records、验收和结果表为准。不能根据旧默认配置推断每个批次CPU/内存。当前Jinlong使用自己的pipeline，保留等价记录即可。

本地快照保留原路径。公开副本替换内部路径/服务器信息并记录前后哈希，所以公开源码是便于审阅的移植副本，需设置数据根后使用，不宣称原样可重提Spark任务。软件/第三方模块保留原版权与许可证；文献及数据库许可分别引用，不给第三方内容重新授予许可。
