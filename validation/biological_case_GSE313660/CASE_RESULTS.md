# GSE313660 human follicular-fluid application: analysis record

Status: **exploratory draft; do not copy therapeutic claims into the manuscript**.

## Input and processing

- Source: [GSE313660](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313660), the human follicular-fluid/blood CITE-seq study of 4 PCOS and 4 control participants.
- Local input: `../../../RawData/GSE313660/GSE313660_Wegrzynowicz_PCOS_GEO.Rds`. The separate raw TSV count file does not carry usable gene row identifiers, so RNA and ADT raw count matrices were extracted from the author's RDS. Author cell-type labels were carried only for post hoc review.
- Follicular-fluid cells were independently normalized, clustered (Seurat, 2,500 variable genes, 30 PCs, resolution 0.5) and projected with UMAP. The frozen human PCOS-CAM dictionary and public `pcos_cam()` workflow were then run on per-cluster RNA expression, yielding 15,768 cells and 17 clusters.
- Reproduction entry points: `run_01_cluster.R`, `show_clustering_annotation.R`, `run_02_myeloid.R`, `run_03_myeloid_trajectory.R`, `run_04_gwas_screen.py`, `run_05_camera_pathways.R`, `run_06_dgidb_lookup.py`, and `plot_figure6_myeloid.R`. R packages needed by these scripts are installed locally.

## Figures and annotation audit

- `clustering_annotation_4panel.pdf` combines the independently derived UMAP, PCOS-CAM family calls, published labels, and PCOS-CAM full calls. The corresponding per-cluster table is `clustering_annotation_review.csv`.
- PCOS-CAM emitted 11 fine and 6 family-only calls. Several fine calls do not agree with published labels: clusters 4 and 5 are predominantly published *Theca*, clusters 7 and 9 are predominantly erythroblast, and cluster 13 is predominantly B cell. The macrophage-rich cluster 8 and granulosa-rich cluster 14 have broadly concordant calls. The published labels are an audit comparator, not a proven ground truth for every cell.
- `Figure_6_myeloid_application.pdf` shows RNA markers, selected surface-protein patterns, an exploratory Slingshot ordering, RNA trends along that ordering, and donor coverage. The trajectory uses 2,026 cells from clusters 6, 8 and 12. It is a cell-state ordering, not measured differentiation or a PCOS-specific direction.
- ADT features CD3-UCHT1 and CD14-M5E2 are nearly absent in this extracted follicular-fluid matrix and were not used for the protein display. CD11b, CD16 and HLA-DR are shown descriptively, not as an independent accuracy estimate.

## PCOS-versus-control and target-screening results

| Analysis | Result | Interpretation |
|---|---:|---|
| Macrophage-rich cluster 8 pseudobulk, donors with ≥20 cells | 3 PCOS + 3 controls | Small and unbalanced cell yields |
| RNA differential expression, edgeR quasi-likelihood | 0 genes at FDR <0.05 | No confirmed gene-level PCOS signal |
| TRRUST v2 TF activity, donor-level limma | 0 regulons at FDR <0.05 | No confirmed TF activity shift |
| Reactome, donor-level CAMERA | 0 pathways at FDR <0.05 | No confirmed pathway shift |
| GWAS Catalog PCOS associations | 147 associations; 17 distinct mapped genes represented in macrophage expression | Locus mappings do not establish causality |
| GWAS-mapped genes also passing macrophage DE FDR | 0 | No supported therapeutic target from the intersection |
| DGIdb interactions returned for 17 mapped genes | 391 raw records; 0 prioritized | A drug-gene record does not establish PCOS efficacy |

The separate gene-rank fgsea screen returned nominally striking pathway results, but it did not account for between-gene correlation and did not survive the donor-level CAMERA check. It is archived as `Reactome_fgsea_gene_rank_exploratory.csv` and should not be presented as a confirmed finding.

## Decision for the manuscript

The figures are available for inspection, but this cohort does **not** currently support a PCOS-specific TF, pathway, target or drug claim. In particular, the full PCOS-CAM fine-label map contains visible misannotations, and one PCOS donor dominates the follicular-cell compartment. A stronger, better-balanced human dataset would be needed before calling this a successful end-to-end disease mechanism and drug-screening case.
