# PCOS-CAM

用于多囊卵巢综合征（PCOS）研究，基于卵巢标记字典的细胞簇注释。内置 data/File_1.xlsx，包含 17 篇研究、153 个标签和 1,501 个标签–基因组合。核心加权评分来自 [FPCAM](https://github.com/guojunliu7/FPCAM/tree/e90cf57d61543f655a27131acdc3fa8539456d22)，采用 Apache-2.0 许可证。

本仓库包含应用源码、标记字典、示例表达表和验证脚本。原始单细胞数据、大型分析对象及生成的图表和结果文件未纳入版本控制；验证脚本所需的外部数据来源见各验证目录的 README。

在本目录的 R 会话中启动：

```r
install.packages(c("shiny", "readxl", "ggplot2", "openxlsx")) # 首次安装
shiny::runApp(".")
```

也可以直接在 R 中注释：

```r
source("Global.R")
source("hierarchy.R")
d <- read_pcos_dictionary("data/File_1.xlsx")
x <- read_pcos_expression("www/PCOS_GSE240688_cluster_expression.csv")
h <- read_pcos_hierarchy("data/File_1.xlsx", d, "data/label_hierarchy.csv")
r <- pcos_cam(x, d, h, "Human")
r$annotations
# Seurat 对象中的簇编号须与输入表一致
# obj$PCOS_CAM <- pcos_cluster_labels(obj$seurat_clusters, r)
```

输入为各簇标记基因的非负表达量。支持成对的 cluster_0_Markergene、cluster_0_Express 列，首列为 gene 的宽表，以及 cluster、gene、expression 长表。成对输入保留 cluster_0 等原编号。每簇每基因应有一个表达值。使用精确基因符号，小鼠输入通过工作簿 Mouse_to_human 页转换。

每个命中基因的贡献为其表达量除以该簇所有输入基因的表达总量，再乘以字典 percent 权重。未映射和未命中字典的输入仍计入表达分母。得分是各标签命中贡献之和，输入基因集合会影响得分尺度。零得分标记为 Unassigned，最高分并列标记为 Ambiguous 并列出全部候选。分差 margin 为第一、第二名之差；delta1 和 delta2 沿用上游的中位数差定义。权重和得分均不表示注释正确的概率。

下载工作簿包含逐簇结果、完整得分矩阵、逐基因贡献和未映射的小鼠基因。热图显示各簇前三个正分候选的并集；UMAP 使用用户提供的已有坐标和本次簇注释。独立 PCOS 队列的注释准确率尚待验证。

## PCOS-CAM 分层注释

153 个原始标签通过显式映射归到 13 个父类。Theca / stroma 和 CD45- Immune-like 各自成组。父类使用 Evidence 中所有成功映射的候选基因；同一投票单位在父类内对同一基因只计一次。依次按票数、各单位最佳排名均值和基因符号选择前十个基因，在父类内归一化权重。细分标签继续使用原 17 篇字典权重。

两级分别比较候选，广义标签也可在第二级胜出。默认每级至少命中两个正表达标记，并要求相对分差（第一名减第二名，再除以第一名）达到 0.10。这两个可调的初始阈值尚待独立队列校准。大类零命中或命中数不足返回 Unassigned，大类并列或分差不足返回 Ambiguous。父类通过、细分证据不足时返回父类。diagnostics 保留各级候选、次选、命中基因数、命中字典权重之和及首位基因贡献占比；并列第一时没有唯一标签可供计算覆盖度。

稳定性检查将首位贡献基因从同级所有候选中移除，再比较得分；多个基因并列贡献第一时，分别移除并报告结果。stable 为 FALSE 表示至少一次移除会改变唯一首选。当前版本将其作为复核标志。热图和分差图展示父类竞争。下载包含全字典评分、父类字典及两级逐基因贡献。

当正式注释仅为 `Granulosa`，且全字典唯一首选标签属于该大类下的亚型时，`candidate_subtype` 记录该候选，`annotation` 保留正式注释。统计时单列候选的增量，同一亚型只计一次。对外图表、横轴、图例和图注统一使用方法名 **PCOS-CAM**。
