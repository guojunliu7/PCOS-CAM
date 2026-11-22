# Exploratory biological case: focus on the RNA/ADT-supported myeloid continuum.
# Run from the PCOS-CAM project root after run_01_cluster.R.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(edgeR)
  library(limma);library(fgsea);library(msigdbr);library(decoupleR)})
out <- 'validation/biological_case_GSE313660'
query <- readRDS(file.path(out,'follicular_reclustered.rds'))
cl <- as.character(query$seurat_clusters)
main <- c('6','8','12')
stopifnot(all(main %in% unique(cl)))

# RNA/ADT marker evidence, including negative lineage controls.
rna_markers <- c('CD3D','CD3E','MS4A1','NKG7','LST1','LYZ','S100A8','S100A9',
                 'FCN1','CD14','FCGR3A','C1QA','C1QB','APOE','CD68','HLA-DRA',
                 'HBB','FOXL2','CYP17A1')
adt_markers <- c('Hu.CD3-UCHT1','Hu.CD19','Hu.CD14-M5E2','Hu.CD16',
                 'Hu.CD11b','Hu.CD163','Hu.HLA.DR','Hu.CD45-HI30')
rna_markers <- intersect(rna_markers,rownames(query[['RNA']]))
adt_markers <- intersect(adt_markers,rownames(query[['ADT']]))
ids <- as.character(sort(as.integer(unique(cl))))
rna_data <- LayerData(query,assay='RNA',layer='data')[rna_markers,,drop=FALSE]
adt_data <- LayerData(query,assay='ADT',layer='data')[adt_markers,,drop=FALSE]
summarise_marker <- function(mat,assay) do.call(rbind,lapply(ids,function(id){
  z <- mat[,cl==id,drop=FALSE]
  data.frame(cluster=id,assay=assay,marker=rownames(z),
             mean=Matrix::rowMeans(z),fraction=Matrix::rowMeans(z>0))
}))
write.csv(rbind(summarise_marker(rna_data,'RNA'),summarise_marker(adt_data,'ADT')),
          file.path(out,'marker_RNA_ADT_by_cluster.csv'),row.names=FALSE)

# One macrophage-enriched cluster is used for a patient-level condition contrast.
# Retain donors with at least 20 cells; sparse donors are unsuitable for pseudobulk.
mac_cells <- colnames(query)[cl=='8']
tab <- table(query$patient[match(mac_cells,colnames(query))])
donors <- names(tab)[tab>=20]
donor_condition <- setNames(as.character(query$condition[match(donors,query$patient)]),donors)
stopifnot(sum(donor_condition=='control')>=3,sum(donor_condition=='pcos')>=3)
count_matrix <- LayerData(query,assay='RNA',layer='counts')[,mac_cells,drop=FALSE]
pb <- do.call(cbind,lapply(donors,function(donor){
  Matrix::rowSums(count_matrix[,query$patient[match(mac_cells,colnames(query))]==donor,drop=FALSE])
}))
colnames(pb) <- donors
write.csv(data.frame(patient=donors,condition=donor_condition[donors],
                     n_cells=as.integer(tab[donors])),
          file.path(out,'macrophage_pseudobulk_donors.csv'),row.names=FALSE)
saveRDS(pb,file.path(out,'macrophage_pseudobulk_counts.rds'))

group <- factor(donor_condition[donors],levels=c('control','pcos'))
design <- model.matrix(~group)
y <- DGEList(counts=pb)
y <- y[filterByExpr(y,group=group),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust=TRUE)
test <- glmQLFTest(fit,coef=2)
de <- topTags(test,n=Inf,sort.by='none')$table
de$gene <- rownames(de)
de <- de[order(de$PValue),c('gene','logFC','logCPM','F','PValue','FDR')]
write.csv(de,file.path(out,'macrophage_PCOS_vs_control_DE.csv'),row.names=FALSE)
cat('Macrophage donor cells:',paste(donors,tab[donors],collapse='; '),'\n')
cat('DE genes FDR<0.05:',sum(de$FDR<.05),'\n')

# Signed Reactome enrichment uses the full ranked gene list, not only a cutoff.
rank <- sign(de$logFC)*sqrt(pmax(de$F,0))
names(rank) <- de$gene
rank <- sort(rank[is.finite(rank)],decreasing=TRUE)
sets <- msigdbr(species='Homo sapiens',collection='C2',subcollection='CP:REACTOME')
pathways <- split(sets$gene_symbol,sets$gs_name)
pathway_result <- fgseaMultilevel(pathways,rank,minSize=10,maxSize=500,eps=0)
pathway_result$leadingEdge <- vapply(pathway_result$leadingEdge,
                                    function(z) paste(z,collapse=';'),character(1))
pathway_result <- pathway_result[order(pathway_result$padj),]
write.csv(as.data.frame(pathway_result),
          file.path(out,'Reactome_fgsea_gene_rank_exploratory.csv'),row.names=FALSE)
cat('Exploratory gene-rank fgsea pathways FDR<0.05 (not donor-level inference):',
    sum(pathway_result$padj<.05,na.rm=TRUE),'\n')

# TRRUST v2 human signed regulons. ULM scores are compared at donor level.
net <- read.delim(file.path(out,'TRRUST_v2_human.tsv'),header=FALSE,
                  col.names=c('source','target','mode','pmid'))
net <- net[net$mode %in% c('Activation','Repression'),]
net$mor <- ifelse(net$mode=='Activation',1,-1)
net <- unique(net[,c('source','target','mor')])
# A few TF-target pairs have both signs in TRRUST. Their direction is
# unresolved, so exclude them rather than passing duplicate edges to ULM.
pair_signs <- aggregate(mor~source+target,net,function(z) length(unique(z)))
ambiguous <- paste(pair_signs$source[pair_signs$mor>1],
                   pair_signs$target[pair_signs$mor>1],sep='|')
net <- net[!paste(net$source,net$target,sep='|') %in% ambiguous,]
logcpm <- edgeR::cpm(y,log=TRUE,prior.count=2)
tf <- decoupleR::run_ulm(logcpm,net,.source='source',.target='target',.mor='mor',minsize=5)
tf <- as.data.frame(tf)
tf$patient <- tf$condition
tf$group <- donor_condition[tf$patient]
write.csv(tf,file.path(out,'TRRUST_TF_activity_by_donor.csv'),row.names=FALSE)
wide <- reshape(tf[,c('source','patient','score')],idvar='source',timevar='patient',direction='wide')
activity <- as.matrix(wide[,-1]);rownames(activity) <- wide$source
colnames(activity) <- sub('^score\\.','',colnames(activity))
activity <- activity[,donors,drop=FALSE]
activity <- activity[apply(activity,1,function(z) all(is.finite(z))),,drop=FALSE]
tf_fit <- limma::eBayes(limma::lmFit(activity,design))
tf_de <- limma::topTable(tf_fit,coef=2,n=Inf,sort.by='P')
tf_de$TF <- rownames(tf_de)
write.csv(tf_de,file.path(out,'TRRUST_TF_activity_PCOS_vs_control.csv'),row.names=FALSE)
cat('TF regulons FDR<0.05:',sum(tf_de$adj.P.Val<.05,na.rm=TRUE),'\n')
