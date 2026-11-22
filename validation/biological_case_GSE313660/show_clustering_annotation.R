# Review-only visualization of the independent GSE313660 follicular-fluid re-clustering.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
out <- 'validation/biological_case_GSE313660'
obj <- readRDS(file.path(out,'follicular_reclustered.rds'))
pred <- read.csv(file.path(out,'PCOS_CAM_predictions.csv'),check.names=FALSE)
meta <- obj@meta.data
meta$cell <- rownames(meta)
meta$cluster <- as.character(meta$seurat_clusters)
meta$pcos_cam <- pred$annotation[match(meta$cluster,as.character(pred$cluster))]
meta$pcos_cam_family <- pred$parent_candidates[match(meta$cluster,as.character(pred$cluster))]
emb <- as.data.frame(Embeddings(obj,'umap'))
names(emb)[1:2] <- c('UMAP_1','UMAP_2')
df <- cbind(emb,meta[,c('cell','cluster','pcos_cam','pcos_cam_family',
                       'author_broad','author_fine','patient','condition')])
centres <- aggregate(cbind(UMAP_1,UMAP_2)~cluster,df,median)
most_common <- function(x) {
  t <- sort(table(x),decreasing=TRUE)
  paste0(names(t)[1],' (',sprintf('%.1f',100*as.numeric(t[1])/length(x)),'%)')
}
summary <- do.call(rbind,lapply(sort(unique(df$cluster)),function(id) {
  z <- df[df$cluster==id,]
  data.frame(cluster=as.integer(id),n_cells=nrow(z),
             PCOS_CAM=z$pcos_cam[1],PCOS_CAM_family=z$pcos_cam_family[1],
             annotation_level=pred$level[match(id,as.character(pred$cluster))],
             published_major_fine=most_common(z$author_fine),
             published_major_broad=most_common(z$author_broad),
             n_patients=length(unique(z$patient)))
}))
summary <- summary[order(summary$cluster),]
write.csv(summary,file.path(out,'clustering_annotation_review.csv'),row.names=FALSE)
write.csv(df[,c('cell','UMAP_1','UMAP_2','cluster','pcos_cam',
               'pcos_cam_family','author_broad','author_fine','patient','condition')],
          file.path(out,'clustering_annotation_cells.csv'),row.names=FALSE)

base <- theme_classic(base_size=16)+theme(
  plot.title=element_text(face='bold',size=20),
  axis.title=element_text(face='bold'),
  axis.text=element_text(face='bold',colour='#23313A'),
  legend.title=element_text(face='bold'),
  panel.border=element_rect(colour='#23313A',fill=NA,linewidth=0.8))
draw <- function(col,title,palette=NULL) {
  p <- ggplot(df,aes(UMAP_1,UMAP_2,colour=.data[[col]]))+
    geom_point(size=0.12,alpha=0.72)+
    geom_label(data=centres,aes(UMAP_1,UMAP_2,label=cluster),
               inherit.aes=FALSE,fontface='bold',size=3.7,
               linewidth=0.15,fill='white',alpha=0.90)+
    labs(title=title,colour=NULL,x='UMAP 1',y='UMAP 2')+base+
    theme(legend.text=element_text(size=11),legend.key.height=grid::unit(0.35,'cm'))
  if (!is.null(palette)) p <- p+scale_colour_manual(values=palette)
  p
}
cluster_palette <- setNames(grDevices::hcl.colors(length(unique(df$cluster)),
                                                     palette='Dark 3'),sort(unique(df$cluster)))
family_palette <- c('Immune'='#0072B2','Granulosa'='#009E73',
                    'Luteal'='#CC79A7','Theca / stroma'='#D55E00')
p1 <- draw('cluster','A  Re-clustered follicular-fluid cells',cluster_palette)
p2 <- draw('pcos_cam_family','B  PCOS-CAM broad-family calls',family_palette)
p3 <- draw('author_broad','C  Published study labels')
fig <- p1+p2+p3+plot_layout(widths=c(1,1,1))+
  plot_annotation(title='GSE313660 | 15,768 human follicular-fluid cells | 8 donors',
                  subtitle='Numbers indicate independently re-clustered groups; published labels are used only for review')
ggsave(file.path(out,'clustering_annotation_review.png'),fig,width=18,height=7.6,dpi=220,bg='white')
ggsave(file.path(out,'clustering_annotation_review.pdf'),fig,width=18,height=7.6,bg='white')
level_by_cluster <- pred$level[match(df$cluster,as.character(pred$cluster))]
df$display_call <- ifelse(level_by_cluster=='fine',df$pcos_cam,
                          paste0(df$pcos_cam,' [family only]'))
display_names <- unique(df$display_call)
display_colors <- setNames(grDevices::hcl.colors(length(display_names),
                                                   palette='Dark 3'),display_names)
fine_plot <- ggplot(df,aes(UMAP_1,UMAP_2,colour=display_call))+
  geom_point(size=0.25,alpha=0.75)+
  geom_label(data=centres,aes(UMAP_1,UMAP_2,label=cluster),
             inherit.aes=FALSE,fontface='bold',size=5.1,
             linewidth=0.25,fill='white',alpha=0.95)+
  scale_colour_manual(values=display_colors)+
  labs(title='PCOS-CAM full annotation | GSE313660',
       subtitle='11 clusters received fine labels; 6 reached family level only. Numbers identify re-clustered groups.',
       colour='Predicted label',x='UMAP 1',y='UMAP 2')+
  theme_classic(base_size=19)+
  theme(plot.title=element_text(face='bold',size=23),
        plot.subtitle=element_text(size=14),
        axis.title=element_text(face='bold'),
        axis.text=element_text(face='bold',colour='#23313A'),
        legend.title=element_text(face='bold',size=16),
        legend.text=element_text(size=13),
        legend.key.height=grid::unit(0.6,'cm'),
        panel.border=element_rect(colour='#23313A',fill=NA,linewidth=0.8))
ggsave(file.path(out,'PCOS_CAM_full_annotation.png'),fine_plot,
       width=12,height=8,dpi=250,bg='white')
ggsave(file.path(out,'PCOS_CAM_full_annotation.pdf'),fine_plot,
       width=12,height=8,bg='white')
combined <- (p1+theme(legend.position='none'))+
  p2+p3+
  (draw('display_call','D  PCOS-CAM full labels',display_colors)+
     theme(legend.text=element_text(size=10)))+
  plot_layout(ncol=2)+
  plot_annotation(title='GSE313660 | Human follicular-fluid application',
                  subtitle='15,768 cells from 8 donors; cluster numbers are shared across panels. [family only] marks calls without a fine label.')
ggsave(file.path(out,'clustering_annotation_4panel.png'),combined,
       width=20,height=14,dpi=250,bg='white')
ggsave(file.path(out,'clustering_annotation_4panel.pdf'),combined,
       width=20,height=14,bg='white')
cat('Saved UMAP review and full cluster/cell tables.\n')
