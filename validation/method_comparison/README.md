# PCOS-CAM 方法比较

本目录记录 GSE268919 小鼠全卵巢探索性比较。脚本从 PCOS-CAM 根目录运行。四个位置参数依次为含 RNA data 层和 sample_id、seurat_clusters 元数据的 Seurat RDS、项目的 refined_cluster_annotation.csv、ScType 官方 sctype_score_.R 路径及输出目录。

## 比较口径

- 与前一次 A/B 试算沿用相同的本地参照、纳入条件和父类对照表。原试算用人工挑选的 marker 平均表达面板，这次从单细胞 RNA 层重新取值，两个实验的分数不能直接拼接。
- 固定随机种子 17021。每个样本内每个簇最多抽 40 个细胞，保留更小的组；各方法共享抽样细胞、既有簇和小鼠–人同源表。指标按簇等权，每类 F1 在参照出现的父类上等权平均。
- A 使用完整细分字典后映射父类，B 使用预设阈值（2 个命中、0.10 相对分差）。ScType 和 UCell 使用 B 的父类正标记集；ScType 使用官方实现，按簇求和后采用官方低分拒判规则，无外加负标记。UCell 在每细胞全转录组上计算秩分数（maxRank 1500），对簇内分数取均值并选择最高正分，精确并列或零分保留未判定。这个标签决策是本次比较的约定。各自原生知识库的端到端表现需独立测试。
- 用已有人工注释中置信度高且预先可以明确映射至父类的簇进行评价。每种方法的拒判和并列均计入总分母。参考标签来自本地表达标记分析，可能共享知识来源。单队列结果仅供检错和比较方法行为。

## 复现信息

ScType 官方实现版本 630e15cf1e51f2612eda4ad0406dfb17503fa8c9。UCell 2.17.1，提交 f778e091f19e28d3d05bde134663fd6b2f938267。若需要指定临时 R 包库，可设置 PCOS_BENCH_R_LIB。数据来源 GSE268919。输出包括逐簇预测、实际抽样细胞、各父类召回率与 F1、汇总指标及 PDF/PNG 图。

[ScType 源码](https://github.com/IanevskiAleksandr/sc-type/blob/630e15cf1e51f2612eda4ad0406dfb17503fa8c9/R/sctype_score_.R) · [UCell 源码](https://github.com/carmonalab/UCell/tree/f778e091f19e28d3d05bde134663fd6b2f938267) · [GSE268919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268919)

此目录用于比较共用 PCOS-CAM 标记资源时的评分策略。各工具使用自身资源的主比较见 [native_comparison](../native_comparison/README.md)。
