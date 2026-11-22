# Figure 4 · SingleR with ovarian marker restriction

Run from the PCOS-CAM root:

```sh
Rscript validation/singler_dictionary/compare.R ../../PCOS_JQF/GSE138518/GSE268919_scRNA_refined_annotation/GSE268919_Seurat_refined_annotation.rds
```

Both arms use SingleR 2.6.0, the fine labels of celldex MouseRNAseqData, the sampled cell manifest, cluster mode, and pruned calls. The default arm uses all shared genes for marker selection. The custom arm supplies the intersection of dictionary markers (one-to-one mouse orthologs), test genes, and reference genes via SingleR's restrict argument. The dictionary supplies candidate genes for marker selection. Reference profiles and output labels come from MouseRNAseqData in both arms.

The script verifies all 26 default calls against the saved native run. Both metrics use the same ten reviewed clusters as Figures 2 and 3. cluster_predictions.csv records label, family, and counting flags. marker_genes.csv lists the eligible genes. metrics.csv supplies the plotted counts. The PNG and PDF are saved in Figures_Tables.
