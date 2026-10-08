# 参考下载与后续QC

这是给Jinlong执行的入口，本地没有运行新FASTA下载。使用NCBI官方Datasets工具（固定实际版本）按accession.version获取，避免通过物种名在以后无意换成新清单。

```bash
datasets download genome accession --inputfile download_accessions.txt --include genome --dehydrated --filename reference_candidates.zip
unzip reference_candidates.zip -d reference_candidates
datasets rehydrate --directory reference_candidates
```

以上为待在新环境运行的官方CLI模板。工作目录、并行、缓存及网络重试由Jinlong选择；逐accession命令也已放入download_manifest.tsv。来源：[Datasets下载文档](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)。不加入自动latest替换逻辑；若所需版本不可用，记录问题并决定替换。

逐文件保存下载状态、实际字节数、官方校验文件及SHA256。现有258若要复用，先与reference258_manifest.tsv的本地哈希核对；替换GCA/GCF不能假定文件相同。下载失败不能当QC排除，也不能当该菌株不存在。

统一QC后核对分类（固定NCBI/GTDB与ANI）、污染、contig/N50、同株别名和明显重复；保留短读碎片对prophage调用的影响。建议对正式主要参考优先使用完整度≥90%、污染≤5%，但最终组合由Jinlong决定并记录。重查预选中的条件结构锚点和缺失指标。为每条记录提供final_included/excluded/replaced/review及理由。

正式冻结至少包括：参考ID、taxon、strain/别名、accession.version、BioSample/Project、实际来源、classification版本、QC工具/库版本、完整度/污染/碎片化、FASTA路径别名、SHA256、参考角色、纳排与替换理由。全部数量从这张表计算，不手工填一个总数。
