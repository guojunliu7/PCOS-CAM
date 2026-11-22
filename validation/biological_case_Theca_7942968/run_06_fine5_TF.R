# TRRUST regulon activity in the PCOS-CAM fine-labeled theca/stroma cluster 5.
# Run from PCOS-CAM project root.
suppressPackageStartupMessages({library(edgeR);library(limma);library(decoupleR)})
out <- 'validation/biological_case_Theca_7942968'
count <- readRDS(file.path(out,'fine5_pseudobulk_counts.rds'))
meta <- read.csv(file.path(out,'fine5_pseudobulk_donors.csv'))
meta <- meta[match(colnames(count),meta$patient),]
group <- factor(meta$condition,levels=c('control','pcos'))
design <- model.matrix(~group)
y <- DGEList(counts=count)
y <- y[filterByExpr(y,group=group),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
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
tf$group <- group[match(tf$patient,meta$patient)]
write.csv(tf,file.path(out,'fine5_TRRUST_activity_by_donor.csv'),row.names=FALSE)
wide <- reshape(tf[,c('source','patient','score')],idvar='source',
                timevar='patient',direction='wide')
activity <- as.matrix(wide[,-1]);rownames(activity) <- wide$source
colnames(activity) <- sub('^score\\.','',colnames(activity))
activity <- activity[,meta$patient,drop=FALSE]
activity <- activity[apply(activity,1,function(z) all(is.finite(z))),,drop=FALSE]
fit <- eBayes(lmFit(activity,design))
result <- topTable(fit,coef=2,n=Inf,sort.by='P')
result$TF <- rownames(result)
write.csv(result,file.path(out,'fine5_TRRUST_TF_PCOS_vs_control.csv'),row.names=FALSE)
cat('Fine5 regulons FDR<0.05:',sum(result$adj.P.Val<.05),'\n')
print(head(result[,c('TF','logFC','P.Value','adj.P.Val')],10),row.names=FALSE)
