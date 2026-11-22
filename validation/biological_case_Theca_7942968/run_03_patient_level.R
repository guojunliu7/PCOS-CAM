# Donor-level PCOS vs control comparison of untreated cultured theca cells.
# Run from the PCOS-CAM project root after run_01_cluster.R.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(edgeR)
                                library(limma);library(msigdbr);library(decoupleR)})
out <- 'validation/biological_case_Theca_7942968'
query <- readRDS(file.path(out,'untreated_theca_reclustered.rds'))
donors <- sort(unique(query$patient))
group <- factor(vapply(donors,function(id) unique(query$condition[query$patient==id]),
                       character(1)),levels=c('control','pcos'))
count <- LayerData(query,assay='RNA',layer='counts')
pb <- do.call(cbind,lapply(donors,function(id) {
  Matrix::rowSums(count[,query$patient==id,drop=FALSE])
}))
colnames(pb) <- donors
write.csv(data.frame(patient=donors,condition=group,
                     n_cells=vapply(donors,function(id) sum(query$patient==id),integer(1))),
          file.path(out,'pseudobulk_donors.csv'),row.names=FALSE)
saveRDS(pb,file.path(out,'pseudobulk_counts.rds'))
design <- model.matrix(~group)
y <- DGEList(counts=pb)
y <- y[filterByExpr(y,group=group),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust=TRUE)
test <- glmQLFTest(fit,coef=2)
de <- topTags(test,n=Inf,sort.by='PValue')$table
de$gene <- rownames(de)
de <- de[,c('gene','logFC','logCPM','F','PValue','FDR')]
write.csv(de,file.path(out,'theca_PCOS_vs_control_DE.csv'),row.names=FALSE)
cat('Donor-level DE genes FDR<0.05:',sum(de$FDR<.05),'\n')

v <- voom(y,design,plot=FALSE)
sets <- msigdbr(species='Homo sapiens',collection='C2',subcollection='CP:REACTOME')
paths <- split(sets$gene_symbol,sets$gs_name)
paths <- lapply(paths,function(z) intersect(unique(z),rownames(v)))
paths <- paths[lengths(paths)>=10 & lengths(paths)<=500]
pathway <- camera(v,index=paths,design=design,contrast=c(0,1),
                  inter.gene.cor=NA,sort=TRUE)
pathway$pathway <- rownames(pathway)
write.csv(pathway,file.path(out,'Reactome_camera_patient_level.csv'),row.names=FALSE)
cat('Donor-level CAMERA pathways FDR<0.05:',sum(pathway$FDR<.05),'\n')

net <- read.delim('validation/biological_case_GSE313660/TRRUST_v2_human.tsv',
                  header=FALSE,col.names=c('source','target','mode','pmid'))
net <- net[net$mode %in% c('Activation','Repression'),]
net$mor <- ifelse(net$mode=='Activation',1,-1)
net <- unique(net[,c('source','target','mor')])
pair_signs <- aggregate(mor~source+target,net,function(z) length(unique(z)))
ambiguous <- paste(pair_signs$source[pair_signs$mor>1],
                   pair_signs$target[pair_signs$mor>1],sep='|')
net <- net[!paste(net$source,net$target,sep='|') %in% ambiguous,]
tf <- as.data.frame(decoupleR::run_ulm(edgeR::cpm(y,log=TRUE,prior.count=2),net,
                                       .source='source',.target='target',.mor='mor',
                                       minsize=5))
tf$patient <- tf$condition
tf$group <- group[match(tf$patient,donors)]
write.csv(tf,file.path(out,'TRRUST_TF_activity_by_donor.csv'),row.names=FALSE)
wide <- reshape(tf[,c('source','patient','score')],idvar='source',
                timevar='patient',direction='wide')
activity <- as.matrix(wide[,-1]);rownames(activity) <- wide$source
colnames(activity) <- sub('^score\\.','',colnames(activity))
activity <- activity[,donors,drop=FALSE]
activity <- activity[apply(activity,1,function(z) all(is.finite(z))),,drop=FALSE]
tf_fit <- eBayes(lmFit(activity,design))
tf_de <- topTable(tf_fit,coef=2,n=Inf,sort.by='P')
tf_de$TF <- rownames(tf_de)
write.csv(tf_de,file.path(out,'TRRUST_TF_activity_PCOS_vs_control.csv'),row.names=FALSE)
cat('Donor-level TF regulons FDR<0.05:',sum(tf_de$adj.P.Val<.05),'\n')
