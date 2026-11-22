# Figure 5 · SingleR ovarian reference and per-label marker comparison

The reference is GSE267729, a published mouse ovarian aging dataset separate from the GSE268919 query. The downloaded H5AD contains 53,119 cells and 11 author cell labels, including three granulosa labels. Dataset DOI https://doi.org/10.1073/pnas.2600323123. Versioned file https://datasets.cellxgene.cziscience.com/e6509f9c-0a3d-463f-bfa0-a6e368d78c37.h5ad.

Run from the PCOS-CAM root:

```sh
Rscript validation/singler_ovary/build_mapping.R
Rscript validation/singler_ovary/run.R ../../PCOS_JQF/GSE138518/GSE268919_scRNA_refined_annotation/GSE268919_Seurat_refined_annotation.rds
```

Both ovarian arms use identical reference cells and 11 author labels. Sample up to 50 cells per donor and label with seed 17021. SingleR (ovary) selects markers automatically. SingleR + PCOS−CAM uses a named genes list. For named granulosa states, pool markers from granulosa dictionary entries containing the author's defining gene (AMH, INHBB or NPPC). These are marker-based correspondences rather than proof of subtype equivalence. Other labels use the indicated lineage markers. M/B uses pooled monocyte, macrophage and B-cell markers and stays a mixed label. Marker lists use one-to-one mouse orthologs present in both matrices. SingleR uses gene identities; PCOS-CAM vote weights are outside this comparison.

SingleR is the general-reference baseline, using celldex MouseRNAseqData and automatic markers. All arms use identical query cells and cluster labels, SingleR cluster mode and pruned predictions. The 26 baseline calls are checked against Figure 4. The two metrics use ten prespecified reviewed clusters. Count each distinct finer label once among clusters matching the reference family. For the ovarian atlas the three granulosa states, M/B and NK/T count as finer than their parent family; each mixed immune population counts as one label. This describes output resolution. Subtype agreement needs independent finer reference labels.

The source mapping, actual mouse markers, sampled reference cells, predictions for each cluster and metrics accompany the analysis. PDF and PNG are in Figures_Tables.
