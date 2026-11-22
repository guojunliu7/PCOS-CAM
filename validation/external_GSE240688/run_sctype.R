suppressPackageStartupMessages({library(Seurat);library(dplyr);library(openxlsx);library(HGNChelper)})
out <- 'validation/external_GSE240688'
repo <- 'validation/sctype_dictionary/upstream'
obj <- readRDS('../../PCOS_JQF/GSE138518/GSE240688_scRNA_refined_annotation/GSE240688_Seurat_refined_annotation.rds')
cells <- read.csv(file.path(out,'sampled_cells.csv'),colClasses='character')$cell
small <- subset(obj,cells=cells);rm(obj);gc()
source(file.path(repo,'R/gene_sets_prepare.R'))
source(file.path(repo,'R/sctype_score_.R'))
code <- readLines(file.path(repo,'R/auto_detect_tissue_type.R'))
code <- sub('seuratObject[[assay]]$data_type','LayerData(seuratObject, assay=assay, layer=data_type)',code,fixed=TRUE)
eval(parse(text=code))
auto_original <- auto_detect_tissue_type
auto_detect_tissue_type <- function(...) {
  z <- auto_original(...)
  write.csv(z,file.path(out,'ScType_tissue_scores.csv'),row.names=FALSE)
  z
}
code <- readLines(file.path(repo,'R/sctype_wrapper.R'))
code <- sub('seurat_object[[assay]]$scale.data','LayerData(seurat_object, assay=assay, layer="scale.data")',code,fixed=TRUE)
eval(parse(text=code))
sctype_source <- function() file.path(repo,'ScTypeDB_full.xlsx')
pdf(tempfile(fileext='.pdf'))
res <- run_sctype(small)
dev.off()
p <- unique(data.frame(cluster=as.character(res$seurat_clusters),ScType=res$sctype_classification))
write.csv(p,file.path(out,'ScType_predictions.csv'),row.names=FALSE)
cat('ScType labels',nrow(p),'tissue',read.csv(file.path(out,'ScType_tissue_scores.csv'))$tissue[1],'\n')
