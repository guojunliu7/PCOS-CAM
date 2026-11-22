# Run from PCOS-CAM root. External cohort; reference subtypes were locally inferred.
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(SingleR);library(celldex)})
source('Global.R');source('hierarchy.R')
out <- 'validation/external_GSE240688'
obj <- readRDS('../../PCOS_JQF/GSE138518/GSE240688_scRNA_refined_annotation/GSE240688_Seurat_refined_annotation.rds')
meta <- obj[[]]
set.seed(17021)
blocks <- split(rownames(meta),interaction(meta$sample_id,meta$seurat_clusters,drop=TRUE))
cells <- unlist(lapply(blocks,function(ids) if(length(ids)>40) sample(ids,40) else ids),use.names=FALSE)
cl <- as.character(meta[cells,'seurat_clusters'])
m <- LayerData(obj,assay='RNA',layer='data')[,cells,drop=FALSE]
stopifnot(identical(colnames(m),cells))
ids <- as.character(sort(as.integer(unique(cl))))
write.csv(data.frame(cell=cells,cluster=cl,sample_id=meta[cells,'sample_id'],condition=meta[cells,'condition']),
          file.path(out,'sampled_cells.csv'),row.names=FALSE)
saveRDS(list(m=m,cl=cl),file.path(out,'sampled_input.rds'))

# Full transcriptome input is essential: PCOS-CAM's denominator includes all genes.
lin <- m;lin@x <- expm1(lin@x)
ave <- do.call(rbind,lapply(ids,function(id) data.frame(cluster=id,gene=rownames(lin),
  expression=Matrix::rowMeans(lin[,cl==id,drop=FALSE]))))
write.csv(ave,file.path(out,'cluster_mean_expression.csv'),row.names=FALSE)
d <- read_pcos_dictionary('data/File_1.xlsx')
h <- read_pcos_hierarchy('data/File_1.xlsx',d,'data/label_hierarchy.csv')
pc <- pcos_cam(ave,d,h,'Human')
write.csv(pc$annotations,file.path(out,'PCOS_CAM_predictions.csv'),row.names=FALSE)

ref <- HumanPrimaryCellAtlasData()
sr <- SingleR(test=m,ref=ref,labels=ref$label.fine,clusters=cl)
write.csv(data.frame(cluster=rownames(sr),SingleR=ifelse(is.na(sr$pruned.labels),'Unknown',sr$pruned.labels),
                     unpruned=sr$labels),file.path(out,'SingleR_predictions.csv'),row.names=FALSE)

local <- read.csv('../../PCOS_JQF/GSE138518/GSE240688_scRNA_refined_annotation/refined_cluster_annotation.csv')
local$seurat_clusters <- as.character(local$seurat_clusters)
audit <- merge(local,pc$annotations,by.x='seurat_clusters',by.y='cluster',all.x=TRUE)
audit <- merge(audit,data.frame(seurat_clusters=rownames(sr),SingleR=ifelse(is.na(sr$pruned.labels),'Unknown',sr$pruned.labels)),
               by='seurat_clusters',all.x=TRUE)
write.csv(audit[order(as.integer(audit$seurat_clusters)),],file.path(out,'cluster_audit.csv'),row.names=FALSE)
cat('cells',length(cells),'clusters',length(ids),'PCOS labels',length(unique(pc$annotations$annotation)),
    'SingleR labels',length(unique(sr$pruned.labels)),'\n')
