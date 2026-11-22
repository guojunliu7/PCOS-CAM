# Patient-level competitive pathway test, accounting for inter-gene correlation.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(edgeR);library(limma);library(msigdbr)})
out <- 'validation/biological_case_GSE313660'
counts <- readRDS(file.path(out,'macrophage_pseudobulk_counts.rds'))
meta <- read.csv(file.path(out,'macrophage_pseudobulk_donors.csv'))
meta <- meta[match(colnames(counts),meta$patient),]
stopifnot(identical(meta$patient,colnames(counts)))
group <- factor(meta$condition,levels=c('control','pcos'))
design <- model.matrix(~group)
y <- DGEList(counts=counts)
y <- y[filterByExpr(y,group=group),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
v <- voom(y,design,plot=FALSE)
sets <- msigdbr(species='Homo sapiens',collection='C2',subcollection='CP:REACTOME')
paths <- split(sets$gene_symbol,sets$gs_name)
paths <- lapply(paths,function(z) intersect(unique(z),rownames(v)))
paths <- paths[lengths(paths)>=10 & lengths(paths)<=500]
result <- camera(v,index=paths,design=design,contrast=c(0,1),
                 inter.gene.cor=NA,sort=TRUE)
result$pathway <- rownames(result)
write.csv(result,file.path(out,'Reactome_camera_patient_level.csv'),row.names=FALSE)
cat('Tested:',nrow(result),'FDR < 0.05:',sum(result$FDR<0.05,na.rm=TRUE),'\n')
print(head(result[,c('pathway','Direction','PValue','FDR')],10),row.names=FALSE)
