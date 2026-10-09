# GitHub交付范围

目标为用户指定的公开仓库sendi-YUKI/p0075e2-CRCprophage-Sendi，保留既有内容。公开内容包括当前研究计划、候选与下载清单、公共来源/质量摘要、参数问题和记录模板、pilot关键表图/图源数据、方法代码/配置的安全副本。

GitHub交接内容直接放在仓库根目录，使用00_project、01_reference等目录，不再套CRC-PHIRE-TEST/analysis2。研究者本地仍保留analysis2目录。后续发布脚本沿用这一根目录布局。

不上传：258份FASTA、患者FASTQ/组装/BAM、数据库/索引、原始公开API缓存、完整原作者附件/采集缓存、本地旧文档归档、运维/调度日志、私钥、token、内部地址、Spark实时状态或其他人的作业信息。已有258可按reference258_manifest中的accession重新获取，或由研究者单独传FASTA并核对哈希；新增参考由Jinlong下载。

公开方法快照中的路径改为逻辑别名，原始本地文件不改；文件前后哈希、是否脱敏和体积记录在public_file_manifest.tsv。历史验收hash对应原始结果，不能拿它验证已脱敏的路径文本，公开文件完整性请用公开manifest。数值数据及科学参数不改变。

历史报告若引用未公开的大文件，按large_files_local_index和新large_file_index模板定位/重建。内部连接配置通过私下安全渠道填写，不从GitHub获得密码或服务器地址。候选与队列元数据为公开来源整理，保留来源及引用；未知疾病标签不补为健康。

GitHub网页以MD为入口。HTML供下载本地阅读，不自动启用Pages。此包不包含新服务器已运行结果或最终QC后参考清单。
