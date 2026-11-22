# Figure 6: descriptive myeloid application, with marker evidence and donor coverage.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2)
                                library(patchwork);library(slingshot)})
out <- 'validation/biological_case_GSE313660'
obj <- readRDS(file.path(out,'follicular_reclustered.rds'))
traj <- read.csv(file.path(out,'myeloid_slingshot_pseudotime.csv'))
fit <- readRDS(file.path(out,'myeloid_slingshot_fit.rds'))
clusters <- c('0','4','6','8','10','12','14')
meta <- obj@meta.data
meta$cluster <- as.character(meta$seurat_clusters)
rna_features <- c('CD3D','MS4A1','FCN1','CD14','FCGR3A','C1QA','HLA-DRA','FOXL2','CYP17A1')
adt_features <- c('Hu.CD11b','Hu.CD16','Hu.HLA.DR')
rna <- LayerData(obj,assay='RNA',layer='data')[rna_features,,drop=FALSE]
adt <- LayerData(obj,assay='ADT',layer='data')[adt_features,,drop=FALSE]
summarise <- function(mat,assay) do.call(rbind,lapply(clusters,function(id) {
  z <- mat[,meta$cluster==id,drop=FALSE]
  data.frame(cluster=id,marker=rownames(z),assay=assay,
             mean=Matrix::rowMeans(z),fraction=Matrix::rowMeans(z>0))
}))
rna_sum <- summarise(rna,'RNA');adt_sum <- summarise(adt,'ADT')
write.csv(rbind(rna_sum,adt_sum),file.path(out,'Figure6_marker_summary.csv'),row.names=FALSE)
rna_sum$marker <- factor(rna_sum$marker,levels=rev(rna_features))
adt_sum$marker <- factor(adt_sum$marker,levels=rev(adt_features))
rna_sum$cluster <- factor(rna_sum$cluster,levels=clusters)
adt_sum$cluster <- factor(adt_sum$cluster,levels=clusters)
rna_sum$relative <- ave(rna_sum$mean,rna_sum$marker,FUN=function(x) {
  if (sd(x)==0) rep(0,length(x)) else as.numeric(scale(x))
})
adt_sum$relative <- ave(adt_sum$mean,adt_sum$marker,FUN=function(x) {
  if (sd(x)==0) rep(0,length(x)) else as.numeric(scale(x))
})
theme_fig <- theme_classic(base_size=15)+theme(
  plot.title=element_text(face='bold',size=19),
  axis.title=element_text(face='bold'),
  axis.text=element_text(colour='#26323B',face='bold'),
  panel.border=element_rect(colour='#26323B',fill=NA,linewidth=.7))
pA <- ggplot(rna_sum,aes(cluster,marker,size=fraction,colour=relative))+
  geom_point()+scale_colour_gradient2(low='#3267AF',mid='#F4F5F6',high='#C34745',
                                      midpoint=0,limits=c(-2.5,2.5),oob=scales::squish)+
  scale_size(range=c(1,9),limits=c(0,1))+
  labs(title='A  RNA identity markers',x='Cluster',y=NULL,
       size='Fraction\ndetected',colour='Row-scaled\nmean')+theme_fig
pB <- ggplot(adt_sum,aes(cluster,marker,fill=relative))+
  geom_tile(colour='white',linewidth=.8)+
  scale_fill_gradient2(low='#3267AF',mid='#F4F5F6',high='#C34745',
                       midpoint=0,limits=c(-2.5,2.5),oob=scales::squish)+
  labs(title='B  Surface-protein patterns',x='Cluster',y=NULL,
       fill='Row-scaled\nCLR mean')+theme_fig+
  theme(panel.grid=element_blank())

curves <- slingCurves(fit)
curve_data <- do.call(rbind,lapply(seq_along(curves),function(i) {
  z <- curves[[i]]$s[curves[[i]]$ord,,drop=FALSE]
  data.frame(UMAP_1=z[,1],UMAP_2=z[,2],lineage=i)
}))
traj$cluster <- factor(traj$cluster,levels=c('6','12','8'))
pal <- c('6'='#0072B2','12'='#E69F00','8'='#CC79A7')
pC <- ggplot(traj,aes(UMAP_1,UMAP_2,colour=cluster))+
  geom_point(size=.45,alpha=.7)+
  geom_path(data=curve_data,aes(UMAP_1,UMAP_2,group=lineage),inherit.aes=FALSE,
            linewidth=1.2,colour='#26323B')+
  scale_colour_manual(values=pal)+
  labs(title='C  Myeloid state ordering',x='UMAP 1',y='UMAP 2',colour='Cluster')+theme_fig

genes <- c('FCGR3A','C1QA','HLA-DRA','IRF1')
cells <- traj$cell
expr <- LayerData(obj,assay='RNA',layer='data')[genes,cells,drop=FALSE]
pt <- traj$pseudotime_lineage_1
pt <- pt/max(pt,na.rm=TRUE)
breaks <- unique(quantile(pt,probs=seq(0,1,length.out=21),na.rm=TRUE))
bin <- cut(pt,breaks=breaks,include.lowest=TRUE,labels=FALSE)
trend <- do.call(rbind,lapply(seq_along(genes),function(i) {
  z <- tapply(as.numeric(expr[i,]),bin,mean,na.rm=TRUE)
  data.frame(gene=genes[i],bin=as.integer(names(z)),mean=as.numeric(z))
}))
trend$pseudotime <- (trend$bin-.5)/max(trend$bin)
trend$gene <- factor(trend$gene,levels=genes)
write.csv(trend,file.path(out,'Figure6_pseudotime_marker_trends.csv'),row.names=FALSE)
pD <- ggplot(trend,aes(pseudotime,mean,colour=gene))+
  geom_line(linewidth=.55,alpha=.35)+
  geom_smooth(method='loess',formula=y~x,se=FALSE,span=.55,linewidth=1.3)+
  facet_wrap(~gene,scales='free_y',ncol=2)+
  scale_colour_manual(values=c('FCGR3A'='#D69E00','C1QA'='#CC79A7',
                                'HLA-DRA'='#0072B2','IRF1'='#009E73'))+
  labs(title='D  RNA trends along ordering',x='Scaled pseudotime',
       y='Mean log-normalized RNA')+theme_fig+theme(legend.position='none')

donors <- c('Control A','Control B','Control C','Control D',
            'PCOS A','PCOS B','PCOS C','PCOS D')
count <- as.data.frame(table(patient=factor(traj$patient,levels=donors),
                             cluster=traj$cluster))
names(count)[3] <- 'n_cells'
write.csv(count,file.path(out,'Figure6_myeloid_patient_coverage.csv'),row.names=FALSE)
pE <- ggplot(count,aes(cluster,patient,fill=log10(n_cells+1)))+
  geom_tile(colour='white',linewidth=.8)+
  geom_text(aes(label=n_cells),size=3.6,fontface='bold')+
  scale_fill_gradient(low='#F5F8FA',high='#146C94')+
  labs(title='E  Donor coverage',x='Cluster',y=NULL,
       fill='log10\n(cells + 1)')+theme_fig+
  theme(panel.grid=element_blank())
top <- pA+pB+plot_layout(widths=c(1.35,.85))
bottom <- pC+pD+pE+plot_layout(widths=c(1,1,1))
figure <- top/bottom+plot_layout(heights=c(1,1.08))+
  plot_annotation(title='GSE313660 | Exploratory myeloid application',
                  subtitle='RNA/ADT patterns, inferred cell-state ordering and individual donor coverage; no treatment or causal inference')
ggsave(file.path(out,'Figure_6_myeloid_application.png'),figure,
       width=21,height=13.5,dpi=230,bg='white')
ggsave(file.path(out,'Figure_6_myeloid_application.pdf'),figure,
       width=21,height=13.5,bg='white')
