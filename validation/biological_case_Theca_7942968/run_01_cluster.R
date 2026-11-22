# Independent human theca-cell pilot, using untreated (-C) cultures only.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
source('Global.R');source('hierarchy.R')
out <- 'validation/biological_case_Theca_7942968'
raw <- '../RawData/Zenodo7942968_Theca'
dirs <- sort(list.dirs(raw,full.names=TRUE,recursive=FALSE))
dirs <- dirs[grepl('Mc[0-9]+-C_',basename(dirs))]
stopifnot(length(dirs)==10)
pcos_ids <- c('Mc03','Mc10','Mc16','Mc26','Mc27')
control_ids <- c('Mc02','Mc06','Mc31','Mc40','Mc50')
objects <- lapply(dirs,function(directory) {
  sample <- sub('.*_(Mc[0-9]+)-C_.*','\\1',basename(directory))
  stopifnot(sample %in% c(pcos_ids,control_ids))
  count <- Read10X(file.path(directory,'filtered_feature_bc_matrix'))
  if (is.list(count)) count <- count[['Gene Expression']]
  obj <- CreateSeuratObject(counts=count,project=sample,min.cells=3,min.features=200)
  obj$patient <- sample
  obj$condition <- if (sample %in% pcos_ids) 'pcos' else 'control'
  obj$percent.mt <- PercentageFeatureSet(obj,pattern='^MT-')
  obj <- subset(obj,subset=nFeature_RNA>=200 & percent.mt<5)
  cat(sample,ncol(obj),'cells after QC\n')
  obj
})
names(objects) <- vapply(objects,function(x) unique(x$patient),character(1))
query <- merge(objects[[1]],y=objects[-1],add.cell.ids=names(objects),project='Theca_7942968')
query <- JoinLayers(query,assay='RNA')
rm(objects);gc()
DefaultAssay(query) <- 'RNA'
set.seed(7942968)
query <- NormalizeData(query,verbose=FALSE)
query <- FindVariableFeatures(query,nfeatures=2500,verbose=FALSE)
query <- ScaleData(query,features=VariableFeatures(query),verbose=FALSE)
query <- RunPCA(query,features=VariableFeatures(query),npcs=40,verbose=FALSE)
query <- FindNeighbors(query,dims=1:30,verbose=FALSE)
query <- FindClusters(query,resolution=.5,verbose=FALSE)
query <- RunUMAP(query,dims=1:30,seed.use=7942968,verbose=FALSE)
saveRDS(query,file.path(out,'untreated_theca_reclustered.rds'))
cl <- as.character(query$seurat_clusters)
write.csv(data.frame(cell=colnames(query),patient=query$patient,
                     condition=query$condition,cluster=cl,
                     nFeature_RNA=query$nFeature_RNA,percent_mt=query$percent.mt),
          file.path(out,'cells_clusters_qc.csv'),row.names=FALSE)
write.csv(as.data.frame(table(patient=query$patient,condition=query$condition,
                              cluster=cl)),file.path(out,'patient_cluster_counts.csv'),row.names=FALSE)
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
write.csv(calls$annotations,file.path(out,'PCOS_CAM_predictions.csv'),row.names=FALSE)
write.csv(calls$diagnostics,file.path(out,'PCOS_CAM_diagnostics.csv'),row.names=FALSE)
cat('Total cells:',ncol(query),'clusters:',length(ids),'\n')
print(calls$annotations[,c('cluster','annotation','level')],row.names=FALSE)
