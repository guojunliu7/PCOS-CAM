# Figure 3 · ScType marker database comparison

Run from the PCOS-CAM root.

```sh
Rscript validation/sctype_dictionary/compare.R validation/sctype_dictionary/upstream ../../PCOS_JQF/GSE138518/GSE268919_scRNA_refined_annotation/GSE268919_Seurat_refined_annotation.rds
```

Both arms use the same 4,913 sampled cells, RNA scale.data, clusters, ScType scoring and cluster rejection threshold. The official database uses the tissue selected by the original automatic run (Hippocampus). All 26 default calls are checked against that run. The ovarian database supplies positive mouse marker sets through the dictionary's one-to-one ortholog mapping. ScType supplies its own marker weighting. PCOS-CAM vote weights and its candidate rule are not applied.

The left panel counts broad-family matches across the ten prespecified clusters. The right panel counts distinct fine labels within matching families. This measures annotation breadth conditional on family agreement. Individual subtype correctness requires separate validation.

Source: https://github.com/IanevskiAleksandr/sc-type/tree/630e15cf1e51f2612eda4ad0406dfb17503fa8c9 . Required source files and the database are saved in upstream. The analysis reads the existing sample manifest, native ScType calls, tissue scores and review audit. Outputs are cluster_predictions.csv, metrics.csv, and Figure_3.pdf/png in Figures_Tables.
