# GSE313660 follicular-fluid application. Run from the PCOS-CAM project root.
# Author labels are retained only for an independent post-hoc concordance audit.
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
source('Global.R')
source('hierarchy.R')

out <- 'validation/biological_case_GSE313660'
dir.create(out, showWarnings=FALSE, recursive=TRUE)
src <- '../RawData/GSE313660/GSE313660_Wegrzynowicz_PCOS_GEO.Rds'
original <- readRDS(src)
meta <- original[[]]
stopifnot(identical(colnames(original), rownames(meta)))
keep <- rownames(meta)[meta$tissue == 'folcel' & !is.na(meta$pcos_status) & !is.na(meta$patient)]
stopifnot(length(keep)>1000, length(unique(meta[keep,'patient']))==8)

counts <- LayerData(original, assay='RNA', layer='counts')[,keep,drop=FALSE]
stopifnot(identical(colnames(counts),keep), all(counts@x == floor(counts@x)))
info <- meta[keep,c('patient','pcos_status','tissue','identities','broad'),drop=FALSE]
names(info) <- c('patient','condition','tissue','author_fine','author_broad')
query <- CreateSeuratObject(counts=counts, meta.data=info, project='GSE313660_FF')
adt <- LayerData(original, assay='ADT', layer='counts')[,keep,drop=FALSE]
stopifnot(identical(colnames(adt),colnames(query)))
query[['ADT']] <- CreateAssay5Object(counts=adt)
query <- NormalizeData(query, assay='ADT', normalization.method='CLR', margin=2, verbose=FALSE)
rm(original, counts); gc()
DefaultAssay(query) <- 'RNA'
set.seed(313660)
query <- NormalizeData(query, verbose=FALSE)
query <- FindVariableFeatures(query, selection.method='vst', nfeatures=2500, verbose=FALSE)
query <- ScaleData(query, features=VariableFeatures(query), verbose=FALSE)
query <- RunPCA(query, features=VariableFeatures(query), npcs=40, verbose=FALSE)
query <- FindNeighbors(query, dims=1:30, verbose=FALSE)
query <- FindClusters(query, resolution=0.5, verbose=FALSE)
query <- RunUMAP(query, dims=1:30, seed.use=313660, verbose=FALSE)
cl <- as.character(query$seurat_clusters)
cat('FF cells:', ncol(query), 'clusters:', length(unique(cl)), '\n')

cell_table <- data.frame(cell=colnames(query), patient=query$patient,
                         condition=query$condition, cluster=cl,
                         author_fine=query$author_fine, author_broad=query$author_broad)
write.csv(cell_table,file.path(out,'cells_and_clusters.csv'),row.names=FALSE)
write.csv(as.data.frame(table(patient=query$patient,condition=query$condition,
                              cluster=cl)),file.path(out,'patient_cluster_counts.csv'),row.names=FALSE)
saveRDS(query,file.path(out,'follicular_reclustered.rds'))

# Use every RNA gene for the PCOS-CAM expression denominator, as in prior benchmarks.
expr <- LayerData(query, assay='RNA', layer='data')
linear <- expr; linear@x <- expm1(linear@x)
ids <- as.character(sort(as.integer(unique(cl))))
means <- do.call(rbind,lapply(ids,function(id) {
  data.frame(cluster=id,gene=rownames(linear),
             expression=Matrix::rowMeans(linear[,cl==id,drop=FALSE]))
}))
saveRDS(means,file.path(out,'cluster_mean_expression.rds'))
d <- read_pcos_dictionary('data/File_1.xlsx')
h <- read_pcos_hierarchy('data/File_1.xlsx',d,'data/label_hierarchy.csv')
calls <- pcos_cam(means,d,h,'Human')
write.csv(calls$annotations,file.path(out,'PCOS_CAM_predictions.csv'),row.names=FALSE)
write.csv(calls$diagnostics,file.path(out,'PCOS_CAM_diagnostics.csv'),row.names=FALSE)
saveRDS(calls,file.path(out,'PCOS_CAM_full_result.rds'))
cat('Predicted labels:\n')
print(calls$annotations[,c('cluster','annotation','level','reason')],row.names=FALSE)
