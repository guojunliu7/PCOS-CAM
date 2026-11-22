# PCOS-CAM

PCOS-CAM annotates cell clusters for polycystic ovary syndrome (PCOS) research using an ovarian marker dictionary. The included `data/File_1.xlsx` covers 17 studies, 153 labels, and 1,501 label–gene combinations. Its core weighted scoring method is derived from [FPCAM](https://github.com/guojunliu7/FPCAM/tree/e90cf57d61543f655a27131acdc3fa8539456d22), which is licensed under Apache-2.0.

This repository contains the application source code, marker dictionary, example expression table, and validation scripts. Raw single-cell data, large analysis objects, and generated figures and result files are not tracked. The README files in the validation directories describe the external data sources required by those scripts.

Run the Shiny app from an R session in this directory:

```r
install.packages(c("shiny", "readxl", "ggplot2", "openxlsx")) # First-time setup
shiny::runApp(".")
```

You can also annotate clusters directly in R:

```r
source("Global.R")
source("hierarchy.R")
d <- read_pcos_dictionary("data/File_1.xlsx")
x <- read_pcos_expression("www/PCOS_GSE240688_cluster_expression.csv")
h <- read_pcos_hierarchy("data/File_1.xlsx", d, "data/label_hierarchy.csv")
r <- pcos_cam(x, d, h, "Human")
r$annotations
# Cluster IDs in a Seurat object must match those in the input table.
# obj$PCOS_CAM <- pcos_cluster_labels(obj$seurat_clusters, r)
```

The input contains non-negative expression values for marker genes in each cluster. Supported layouts are paired `cluster_0_Markergene` and `cluster_0_Express` columns; a wide table with `gene` as its first column; and a long table with `cluster`, `gene`, and `expression` columns. Paired input retains the original cluster IDs, such as `cluster_0`. Each cluster–gene pair should have one expression value. Use exact gene symbols. Mouse genes are converted through the `Mouse_to_human` sheet in the workbook.

For each matched gene, its contribution is its expression divided by the total expression of all input genes in that cluster, multiplied by the dictionary's `percent` weight. Unmapped input genes and genes absent from the dictionary still contribute to the expression denominator. A label's score is the sum of its matched gene contributions, so the set of input genes affects the score scale. A zero score yields `Unassigned`. A tie for the highest score yields `Ambiguous` and lists all tied candidates. The `margin` is the difference between the highest and second-highest scores; `delta1` and `delta2` retain the median-difference definitions of the upstream method. Neither the weights nor the scores are probabilities that an annotation is correct.

The downloadable workbook includes cluster-level results, the full score matrix, per-gene contributions, and unmapped mouse genes. The heatmap shows the union of the top three positively scored candidates for each cluster. The UMAP display uses coordinates supplied by the user and the annotations from the current run. Annotation accuracy on an independent PCOS cohort remains to be established.

## Hierarchical annotation in PCOS-CAM

An explicit mapping assigns the 153 original labels to 13 parent classes. `Theca / stroma` and `CD45- Immune-like` are separate classes. For each parent class, the method considers all successfully mapped candidate genes in the `Evidence` sheet and counts a gene only once per voting unit within that class. It selects the top ten genes by vote count, then by the mean of each unit's best rank, and finally by gene symbol. Weights are normalized within each parent class. Fine-grained labels retain the weights from the original 17-study dictionary.

Candidates are compared separately at the parent and fine-grained levels. A broad label can also win at the second level. By default, each level requires at least two matched markers with positive expression and a relative score margin of at least 0.10, calculated as `(top score - second score) / top score`. These adjustable initial thresholds have not yet been calibrated on an independent cohort. A parent class with no matches or too few matched markers yields `Unassigned`; a tie or insufficient margin yields `Ambiguous`. When the parent class passes but fine-grained evidence is insufficient, the parent class is returned. `diagnostics` retains the leading and runner-up candidates at each level, the number of matched genes, the sum of their dictionary weights, and the fraction of the leading score contributed by its top gene. When first place is tied, there is no unique label for calculating coverage.

For the stability check, the highest-contributing gene is removed from all candidates at the same level before rescoring. If several genes tie for the highest contribution, each is removed in a separate check. `stable = FALSE` means at least one removal changes the unique top candidate. The current version uses this as a review flag. The heatmap and margin plot show competition among parent classes. Downloads include the full-dictionary scores, parent-class dictionary, and per-gene contributions at both levels.

When the final annotation is only `Granulosa` and the unique top label from the full dictionary is a subtype within that class, `candidate_subtype` records that candidate while `annotation` remains `Granulosa`. Summaries report the candidate's added contribution separately and count each subtype only once. Public-facing plots, axes, legends, and captions use the method name **PCOS-CAM** consistently.
