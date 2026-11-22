# Batch-aware clustering of untreated theca cultures; disease labels are not
# supplied to Harmony or to cluster annotation.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(Seurat);library(harmony);library(Matrix)})
source('Global.R');source('hierarchy.R')
out <- 'validation/biological_case_Theca_7942968'
query <- readRDS(file.path(out,'untreated_theca_reclustered.rds'))
set.seed(7942968)
query <- RunHarmony(query,group.by.vars='patient',reduction.use='pca',
                    dims.use=1:30,reduction.save='harmony',project.dim=FALSE,
                    verbose=FALSE)
query <- FindNeighbors(query,reduction='harmony',dims=1:30,verbose=FALSE)
query <- FindClusters(query,resolution=.5,verbose=FALSE)
query <- RunUMAP(query,reduction='harmony',dims=1:30,seed.use=7942968,verbose=FALSE)
saveRDS(query,file.path(out,'untreated_theca_harmony.rds'))
cl <- as.character(query$seurat_clusters)
write.csv(data.frame(cell=colnames(query),patient=query$patient,
                     condition=query$condition,cluster=cl),
          file.path(out,'harmony_cells_clusters.csv'),row.names=FALSE)
tab <- as.data.frame(table(patient=query$patient,condition=query$condition,cluster=cl))
write.csv(tab,file.path(out,'harmony_patient_cluster_counts.csv'),row.names=FALSE)
expr <- LayerData(query,assay='RNA',layer='data')
linear <- expr;linear@x <- expm1(linear@x)
ids <- as.character(sort(as.integer(unique(cl))))
means <- do.call(rbind,lapply(ids,function(id) {
  data.frame(cluster=id,gene=rownames(linear),
             expression=Matrix::rowMeans(linear[,cl==id,drop=FALSE]))
}))
d <- read_pcos_dictionary('data/File_1.xlsx')
h <- read_pcos_hierarchy('data/File_1.xlsx',d,'data/label_hierarchy.csv')
calls <- pcos_cam(means,d,h,'Human')
write.csv(calls$annotations,file.path(out,'harmony_PCOS_CAM_predictions.csv'),row.names=FALSE)
write.csv(calls$diagnostics,file.path(out,'harmony_PCOS_CAM_diagnostics.csv'),row.names=FALSE)
cat('Total cells:',ncol(query),'clusters:',length(ids),'\n')
print(calls$annotations[,c('cluster','annotation','level')],row.names=FALSE)
