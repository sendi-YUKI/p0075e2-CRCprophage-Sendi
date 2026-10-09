# 文件完整性与维护

这里存放整个交接包的文件校验清单，以及恢复历史方法快照和重复图源路径的小工具。研究内容从[项目首页](../README.md)进入。

`public_file_manifest.tsv`逐文件记录公开副本的大小、SHA256、原本地SHA256及是否发生公开路径替换。自身不计入自身清单。`current_validation_report.json`记录本次整理的实际检查。

`05_existing_results/method_sources`采用一份完整源码与两份版本差异存储；三个历史快照都可恢复，不能把差异目录单独当完整pipeline运行：

```sh
python 90_package/tools/compact_pilot_sources.py verify --archive 05_existing_results/method_sources
python 90_package/tools/compact_pilot_sources.py restore --archive 05_existing_results/method_sources --output restored_sources
```

输出目录必须为空或尚不存在。公开代码中的服务器路径已改为逻辑路径；在新环境执行前按实际数据位置配置。工具、数据库和第三方代码沿用各自许可证。

九份逐字相同的图源副本已用`05_existing_results/source_data_index.tsv`连接到底表；其余独特图源数据保留。旧图源manifest是绘图当时的历史记录，按原路径核验前可恢复这些副本：

```sh
python 90_package/tools/restore_figure_source_aliases.py
python 90_package/tools/restore_figure_source_aliases.py --restore
```

当前Markdown是文档唯一维护源，同名HTML由Markdown生成，不分别改写。公开包不含原始测序大文件、FASTA、原始API缓存、旧文档归档或服务器凭据。需要恢复旧方案时查看Git历史，不把旧方案作为当前入口。
