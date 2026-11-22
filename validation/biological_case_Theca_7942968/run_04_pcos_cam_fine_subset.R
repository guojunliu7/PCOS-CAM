# Patient-level exploratory contrast restricted to PCOS-CAM fine-labeled
# Theca & stroma 2 (Wu), Harmony cluster 5, with >=20 cells per donor.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(edgeR)
                                library(limma);library(msigdbr);library(decoupleR)})
out <- 'validation/biological_case_Theca_7942968'
query <- readRDS(file.path(out,'untreated_theca_harmony.rds'))
stopifnot(as.character(read.csv(file.path(out,'harmony_PCOS_CAM_predictions.csv'))$annotation[6])==
          'Theca & stroma 2 (Wu)')
cells <- colnames(query)[as.character(query$seurat_clusters)=='5']
n <- table(query$patient[match(cells,colnames(query))])
donors <- sort(names(n)[n>=20])
group <- factor(vapply(donors,function(id) unique(query$condition[query$patient==id]),
                       character(1)),levels=c('control','pcos'))
stopifnot(sum(group=='control')>=4,sum(group=='pcos')>=5)
count <- LayerData(query,assay='RNA',layer='counts')[,cells,drop=FALSE]
pb <- do.call(cbind,lapply(donors,function(id) {
  Matrix::rowSums(count[,query$patient[match(cells,colnames(query))]==id,drop=FALSE])
}))
colnames(pb) <- donors
write.csv(data.frame(patient=donors,condition=group,n_cells=as.integer(n[donors])),
          file.path(out,'fine5_pseudobulk_donors.csv'),row.names=FALSE)
saveRDS(pb,file.path(out,'fine5_pseudobulk_counts.rds'))
design <- model.matrix(~group)
y <- DGEList(counts=pb)
y <- y[filterByExpr(y,group=group),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust=TRUE)
test <- glmQLFTest(fit,coef=2)
de <- topTags(test,n=Inf,sort.by='PValue')$table
de$gene <- rownames(de)
write.csv(de,file.path(out,'fine5_PCOS_vs_control_DE.csv'),row.names=FALSE)
cat('PCOS-CAM fine5 cells:',length(cells),'donors:',length(donors),
    'DE FDR<0.05:',sum(de$FDR<.05),'\n')
v <- voom(y,design,plot=FALSE)
sets <- msigdbr(species='Homo sapiens',collection='C2',subcollection='CP:REACTOME')
paths <- split(sets$gene_symbol,sets$gs_name)
paths <- lapply(paths,function(z) intersect(unique(z),rownames(v)))
paths <- paths[lengths(paths)>=10 & lengths(paths)<=500]
pathway <- camera(v,index=paths,design=design,contrast=c(0,1),
                  inter.gene.cor=NA,sort=TRUE)
pathway$pathway <- rownames(pathway)
write.csv(pathway,file.path(out,'fine5_Reactome_camera.csv'),row.names=FALSE)
cat('CAMERA pathways FDR<0.05:',sum(pathway$FDR<.05),'\n')
