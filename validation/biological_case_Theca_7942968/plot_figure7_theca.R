# Human cultured-theca application. Run from the PCOS-CAM project root.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(ggrepel)
})
out <- 'validation/biological_case_Theca_7942968'
obj <- readRDS(file.path(out, 'untreated_theca_harmony.rds'))
pred <- read.csv(file.path(out, 'harmony_PCOS_CAM_predictions.csv'))
cells <- read.csv(file.path(out, 'harmony_cells_clusters.csv'))
de <- read.csv(file.path(out, 'fine5_PCOS_vs_control_DE.csv'))
ora <- read.csv(file.path(out, 'fine5_Reactome_DEG_ORA_up.csv'))
tf <- read.csv(file.path(out, 'fine5_TRRUST_TF_PCOS_vs_control.csv'))
donors <- read.csv(file.path(out, 'fine5_pseudobulk_donors.csv'))
stopifnot(nrow(cells) == ncol(obj), all(cells$cell %in% colnames(obj)))

emb <- as.data.frame(Embeddings(obj, 'umap'))
names(emb)[1:2] <- c('x', 'y')
emb$cell <- rownames(emb)
emb <- merge(emb, cells, by='cell', sort=FALSE)
emb$cluster <- factor(emb$cluster, levels=as.character(sort(as.integer(unique(emb$cluster)))))
emb$call_level <- pred$level[match(as.character(emb$cluster), as.character(pred$cluster))]
emb$call_level <- factor(emb$call_level,
                         levels=c('fine','parent','unresolved'),
                         labels=c('Fine subtype','Broad family','Unresolved'))
centres <- aggregate(cbind(x,y)~cluster,emb,median)

ink <- '#20323D'; green <- '#009E73'; orange <- '#E69F00'; blue <- '#0072B2';
purple <- '#9467BD'; grey <- '#A7B0B5'
base <- theme_classic(base_size=16)+theme(
  plot.title=element_text(face='bold',size=19,colour=ink),
  plot.subtitle=element_text(size=12,colour=ink),
  axis.title=element_text(face='bold',size=15,colour=ink),
  axis.text=element_text(face='bold',size=12,colour=ink),
  legend.title=element_text(face='bold',size=13,colour=ink),
  legend.text=element_text(size=11,colour=ink),
  panel.border=element_rect(colour=ink,fill=NA,linewidth=.8),
  plot.margin=margin(12,14,12,14))
cluster_cols <- setNames(grDevices::hcl.colors(nlevels(emb$cluster),'Dark 3'),levels(emb$cluster))
pA <- ggplot(emb,aes(x,y,colour=cluster))+
  geom_point(size=.18,alpha=.75)+
  geom_label(data=centres,aes(x,y,label=cluster),inherit.aes=FALSE,
             fontface='bold',size=3.7,fill='white',linewidth=.18)+
  scale_colour_manual(values=cluster_cols)+
  labs(title='A  Patient-corrected clusters',subtitle='10,446 cells | 10 untreated donors',
       x='UMAP 1',y='UMAP 2',colour='Cluster')+base+
  theme(legend.position='none')
pB <- ggplot(emb,aes(x,y,colour=call_level))+
  geom_point(size=.18,alpha=.78)+
  geom_label(data=centres,aes(x,y,label=cluster),inherit.aes=FALSE,
             fontface='bold',size=3.7,fill='white',linewidth=.18)+
  scale_colour_manual(values=c('Fine subtype'=green,'Broad family'=blue,'Unresolved'=grey))+
  labs(title='B  PCOS-CAM resolution',subtitle='Clusters 3 and 5: Theca & stroma 2 (Wu)',
       x='UMAP 1',y='UMAP 2',colour='Call level')+base

donors$condition <- factor(donors$condition,levels=c('control','pcos'),
                           labels=c('Control','PCOS'))
donors$patient <- factor(donors$patient, levels=donors$patient[order(donors$condition,donors$patient)])
pC <- ggplot(donors,aes(patient,n_cells,fill=condition))+
  geom_col(width=.72,colour=ink,linewidth=.25)+
  scale_fill_manual(values=c('Control'=blue,'PCOS'=orange))+
  labs(title='C  Cluster 5 donor coverage',
       subtitle='Pseudobulk: 5 PCOS + 4 control donors; >=20 cells each',
       x='Donor',y='Cells',fill=NULL)+base+
  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='top')

de$state <- ifelse(!is.na(de$FDR)&de$FDR<.05,
                   ifelse(de$logFC>0,'Higher in PCOS','Lower in PCOS'),'FDR >= 0.05')
de$state <- factor(de$state,levels=c('FDR >= 0.05','Higher in PCOS','Lower in PCOS'))
de$mlog10 <- -log10(pmax(de$FDR,.Machine$double.xmin))
highlight <- de[de$gene %in% c('SOD1','PRDX1','PRDX2','TXN','IQCA1'),]
pD <- ggplot(de,aes(logFC,mlog10,colour=state))+
  geom_point(size=.65,alpha=.72)+
  geom_hline(yintercept=-log10(.05),linetype='dashed',colour=ink,linewidth=.5)+
  geom_text_repel(data=highlight,aes(label=gene),size=3.6,fontface='bold',
                  colour=ink,max.overlaps=Inf,min.segment.length=0,
                  box.padding=.3,seed=7942968,show.legend=FALSE)+
  scale_colour_manual(values=c('FDR >= 0.05'=grey,
                               'Higher in PCOS'=orange,'Lower in PCOS'=blue))+
  labs(title='D  Donor-level differential expression',
       subtitle='Cluster 5 | 618 higher, 407 lower (FDR < 0.05)',
       x='log2 fold change: PCOS / control',y='-log10(FDR)',colour=NULL)+base+
  theme(legend.position='top',legend.text=element_text(size=10))

wanted <- c('REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT',
            'REACTOME_MITOCHONDRIAL_TRANSLATION',
            'REACTOME_CELLULAR_RESPONSES_TO_STRESS',
            'REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS',
            'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES')
path <- ora[ora$ID %in% wanted,]
path$label <- gsub('REACTOME_','',path$ID,fixed=TRUE)
path$label <- gsub('_',' ',path$label,fixed=TRUE)
path$label <- tolower(path$label)
path$label <- paste0(toupper(substr(path$label,1,1)),substring(path$label,2))
path$label <- factor(path$label,levels=rev(path$label[order(path$p.adjust)]))
pE <- ggplot(path,aes(-log10(p.adjust),label))+
  geom_col(fill=purple,width=.65,colour=ink,linewidth=.25)+
  geom_text(aes(label=Count),hjust=-.25,fontface='bold',size=4)+
  scale_x_continuous(expand=expansion(mult=c(0,.18)))+
  labs(title='E  Reactome DEG over-representation',
       subtitle='Upregulated genes; bar labels = overlapping genes',
       x='-log10(ORA FDR)',y=NULL)+base+
  theme(axis.text.y=element_text(size=10,face='bold'))

tf7 <- tf[!is.na(tf$adj.P.Val)&tf$adj.P.Val<.05,]
tf7$TF <- factor(tf7$TF,levels=tf7$TF[order(tf7$logFC)])
tf7$direction <- ifelse(tf7$logFC>0,'Higher in PCOS','Lower in PCOS')
pF <- ggplot(tf7,aes(logFC,TF,fill=direction))+
  geom_col(width=.69,colour=ink,linewidth=.25)+
  geom_vline(xintercept=0,colour=ink,linewidth=.45)+
  scale_fill_manual(values=c('Higher in PCOS'=orange,'Lower in PCOS'=blue))+
  labs(title='F  TRRUST-inferred TF activity',
       subtitle='Donor-level regulon scores; FDR < 0.05',
       x='Activity difference: PCOS - control',y=NULL,fill=NULL)+base+
  theme(legend.position='top')

fig <- ((pA | pB) / (pC | pD) / (pE | pF))+
  plot_annotation(title='Human PCOS cultured-theca application | Zenodo 7942968',
                  subtitle='PCOS-CAM resolves one donor-covered theca/stroma subtype for disease analysis; cultured cells are not intact ovarian tissue',
                  theme=theme(plot.title=element_text(face='bold',size=25,colour=ink),
                              plot.subtitle=element_text(size=15,colour=ink),
                              plot.background=element_rect(fill='white',colour=NA)))
ggsave(file.path(out,'Figure_7_theca_application.png'),fig,width=18,height=23,dpi=260,bg='white')
ggsave(file.path(out,'Figure_7_theca_application.pdf'),fig,width=18,height=23,bg='white')
write.csv(data.frame(metric=c('all_cells','clusters','fine_label_clusters',
                               'cluster5_cells','PCOS_donors','control_donors',
                               'DE_higher_FDR05','DE_lower_FDR05','TRRUST_TF_FDR05',
                               'Reactome_ORA_up_FDR05','Reactome_CAMERA_FDR05'),
                     value=c(nrow(emb),nlevels(emb$cluster),
                             sum(pred$level=='fine'),sum(emb$cluster=='5'),
                             sum(donors$condition=='PCOS'),sum(donors$condition=='Control'),
                             sum(de$state=='Higher in PCOS'),sum(de$state=='Lower in PCOS'),
                             nrow(tf7),sum(ora$p.adjust<.05,na.rm=TRUE),0)),
          file.path(out,'Figure_7_numbers.csv'),row.names=FALSE)
cat('Saved Figure 7; pathway bars are ORA, not CAMERA.\n')
