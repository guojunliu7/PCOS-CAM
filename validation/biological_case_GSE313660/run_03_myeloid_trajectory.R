# Exploratory state ordering within the myeloid compartment. It is not a
# time course or evidence of literal cell differentiation.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(Seurat);library(slingshot);library(ggplot2)})
out <- 'validation/biological_case_GSE313660'
obj <- readRDS(file.path(out,'follicular_reclustered.rds'))
clusters <- as.character(obj$seurat_clusters)
selected <- colnames(obj)[clusters %in% c('6','8','12')]
sub <- subset(obj,cells=selected)
DefaultAssay(sub) <- 'RNA'
set.seed(313660)
sub <- FindVariableFeatures(sub,nfeatures=2000,verbose=FALSE)
sub <- ScaleData(sub,features=VariableFeatures(sub),verbose=FALSE)
sub <- RunPCA(sub,features=VariableFeatures(sub),npcs=20,verbose=FALSE)
sub <- RunUMAP(sub,dims=1:15,seed.use=313660,verbose=FALSE)
loc <- Embeddings(sub,'umap')[selected,1:2,drop=FALSE]
cl <- factor(clusters[match(selected,colnames(obj))],levels=c('6','12','8'))
stopifnot(length(selected)>1000,all(table(cl)>100))
set.seed(313660)
fit <- slingshot(loc,clusterLabels=cl,start.clus='6',approx_points=200)
ptime <- slingPseudotime(fit)
weights <- slingCurveWeights(fit)
names_pt <- paste0('pseudotime_lineage_',seq_len(ncol(ptime)))
names_wt <- paste0('lineage_weight_',seq_len(ncol(weights)))
result <- data.frame(cell=selected,patient=obj$patient[match(selected,colnames(obj))],
                     condition=obj$condition[match(selected,colnames(obj))],
                     cluster=as.character(cl),UMAP_1=loc[,1],UMAP_2=loc[,2],
                     ptime,weights,check.names=FALSE)
names(result)[7:ncol(result)] <- c(names_pt,names_wt)
write.csv(result,file.path(out,'myeloid_slingshot_pseudotime.csv'),row.names=FALSE)
saveRDS(fit,file.path(out,'myeloid_slingshot_fit.rds'))

curves <- slingCurves(fit)
curve_data <- do.call(rbind,lapply(seq_along(curves),function(i) {
  z <- curves[[i]]$s[curves[[i]]$ord,,drop=FALSE]
  data.frame(UMAP_1=z[,1],UMAP_2=z[,2],lineage=paste0('Lineage ',i))
}))
pal <- c('6'='#0072B2','12'='#E69F00','8'='#CC79A7')
p1 <- ggplot(result,aes(UMAP_1,UMAP_2,colour=cluster))+
  geom_point(size=0.38,alpha=0.65)+
  geom_path(data=curve_data,aes(UMAP_1,UMAP_2,group=lineage),inherit.aes=FALSE,
            colour='#202E38',linewidth=1.3)+
  scale_colour_manual(values=pal,labels=c('6'='Classical monocyte-rich',
                                           '12'='FCGR3A+ monocyte-rich',
                                           '8'='Macrophage-rich'))+
  labs(title='A  Myeloid state continuum',subtitle='Slingshot rooted in cluster 6',
       x='UMAP 1',y='UMAP 2',colour=NULL)+theme_classic(base_size=17)
if (ncol(ptime)>=1) {
  chosen <- if (ncol(ptime)==1) 1 else which.max(colSums(!is.na(ptime)))
  result$chosen_pseudotime <- ptime[,chosen]
  result$chosen_pseudotime <- result$chosen_pseudotime/
    max(result$chosen_pseudotime,na.rm=TRUE)
  p2 <- ggplot(result,aes(UMAP_1,UMAP_2,colour=chosen_pseudotime))+
    geom_point(size=0.38,alpha=0.8,na.rm=TRUE)+
    scale_colour_viridis_c(option='C',na.value='#D9DDE0')+
    labs(title='B  One inferred state ordering',subtitle='Gray cells fall outside this fitted curve',
         x='UMAP 1',y='UMAP 2',colour='Scaled\npseudotime')+
    theme_classic(base_size=17)
} else stop('Slingshot returned no pseudotime lineage')
for (nm in c('p1','p2')) {
  z <- get(nm)+theme(plot.title=element_text(face='bold',size=19),
                      axis.title=element_text(face='bold'),
                      panel.border=element_rect(colour='#202E38',fill=NA,linewidth=0.7))
  assign(nm,z)
}
suppressPackageStartupMessages(library(patchwork))
fig <- p1+p2+plot_annotation(title='GSE313660 | Exploratory myeloid cell-state trajectory',
                             subtitle='Based on independently re-clustered follicular-fluid cells; pseudotime is an ordering, not a measured time axis')
ggsave(file.path(out,'myeloid_trajectory.png'),fig,width=14,height=6.6,dpi=250,bg='white')
ggsave(file.path(out,'myeloid_trajectory.pdf'),fig,width=14,height=6.6,bg='white')
cat('Cells:',length(selected),'lineages:',ncol(ptime),'\n')
