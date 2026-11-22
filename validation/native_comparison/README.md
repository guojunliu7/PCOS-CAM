## 定稿双面板图的保留文件

从 PCOS-CAM 项目根目录运行以下命令。依赖 R 和 ggplot2。

```sh
cd validation/native_comparison
Rscript compare_granulosa_candidates.R
```

以下路径均相对于本 README 所在目录。清理旧结果时请保留这些文件。

| 用途 | 文件 |
|---|---|
| 绘图与统计代码 | `compare_granulosa_candidates.R` |
| 颗粒细胞候选判定代码 | `../../hierarchy.R` |
| 逐簇预测、大类映射、三算法逐簇核对输入 | `../hierarchical_comparison/full_predictions.csv`、`../../data/label_hierarchy.csv`、`all_cluster_audit.csv` |
| 输出图 | `granulosa_candidate_comparison.png`、`granulosa_candidate_comparison.pdf` |
| 输出数据 | `granulosa_candidate_comparison.csv`、`granulosa_six_cluster_review.csv` |

# End-to-end comparison on GSE268919

The primary comparison pairs PCOS-CAM with File_1.xlsx, ScType with its bundled ScTypeDB_full.xlsx, and SingleR with the mouse reference provided by celldex (MouseRNAseqData, fine labels, cluster mode, pruned labels). SingleR requires the investigator to choose a reference. UCell alone has no ovarian annotation reference and is not assigned an end-to-end rank.

The cell IDs, cluster partitions, and 10 reviewed clusters are exactly those of the previous benchmark. PCOS-CAM A and B raw labels are reused from that run. No PCOS-CAM marker sets or orthology tables are passed to external tools. ScType uses its supplied positive and negative gene sets, the Seurat RNA scale.data, its automatic tissue selector and built-in rejection rule. The upstream detector at commit 630e15cf1e51f2612eda4ad0406dfb17503fa8c9 used a literal layer name in Seurat v5; the run script fixes that accessor and pins the source assets locally. The highest-scoring tissue was Hippocampus. The native database has 16 tissues and no Ovary entry. SingleR 2.6.0 and celldex 1.14.0 used the RNA log-normalized expression directly.

The main count is the number of distinct narrower output labels among clusters whose predicted broad family matches the review family. Synonyms across methods are not counted together; within each method repeated identical labels count once. Mesenchymal and stromal are compatible; CD45-negative immune-like and microglial calls are counted as immune at the broad level under this liberal hierarchy. Mixed theca/stroma versus a single stromal reference is shown separately. A second count requires the reference detail or its recorded markers to support the identity at the claimed resolution. For example NK versus T/NK is only partial; a novel A/B subgroup or a hippocampal microglia label is not established by an ovarian immune reference. Distinct supported counts indicate local crosswalk agreement, not independent biological accuracy. All predicted labels, including unreviewed clusters, are in all_cluster_audit.csv.
