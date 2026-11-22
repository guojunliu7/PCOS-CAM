# PCOS-CAM user guide

Select the species of the input gene symbols and upload a cluster expression table. PCOS-CAM reads the ovarian marker dictionary bundled with the application. Scores update when the input or species changes.

The human example from the PCOS study GSE240688 is available from the sidebar. Its first column contains gene symbols and its remaining columns contain average expression by cluster. The example illustrates the input format. It does not provide independent validation labels.

To display an existing UMAP, upload a CSV with cluster, UMAP_1 and UMAP_2 columns. Use the same cluster identifiers as in the expression table. Download the results workbook for cluster annotations, candidate scores, marker contributions and unmapped mouse genes.

Mode B scores parent classes first and resolves labels within the leading parent. The default requirement of two matched markers and a relative margin of 0.10 at each stage is provisional. The diagnostics tab reports coverage and sensitivity to removing the largest contributing gene. Use mode A to inspect the original flat ranking.
