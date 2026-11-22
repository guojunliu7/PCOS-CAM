# Run from PCOS-CAM root, passing upstream ScType checkout and Seurat RDS.
args <- commandArgs(TRUE)
library(Seurat);library(dplyr);library(HGNChelper)
out <- "validation/native_comparison"
cells <- read.csv("validation/method_comparison/sampled_cells.csv", colClasses="character")
obj <- readRDS(args[2])
cl <- cells$cluster
m <- LayerData(obj,assay="RNA",layer="data")[,cells$cell,drop=FALSE]
stopifnot(identical(colnames(m),cells$cell), identical(as.character(obj[[]][cells$cell,"seurat_clusters"]),cl))
# Independent input for SingleR cluster mode, including genes outside our dictionary.
saveRDS(list(m=m,cl=cl),file.path(tempdir(),"native-input.rds"))
file.copy(file.path(tempdir(),"native-input.rds"),"/tmp/pcos-native-input.rds",overwrite=TRUE)
small <- subset(obj,cells=cells$cell);rm(obj);gc()
repo <- normalizePath(args[1]);db <- file.path(repo,"ScTypeDB_full.xlsx")
source(file.path(repo,"R/gene_sets_prepare.R"))
source(file.path(repo,"R/sctype_score_.R"))
source(file.path(repo,"R/sctype_wrapper.R"))
# Fix literal $data_type accessor in the upstream automatic detector for Assay5.
code <- readLines(file.path(repo,"R/auto_detect_tissue_type.R"))
code <- sub("seuratObject[[assay]]$data_type","LayerData(seuratObject, assay=assay, layer=data_type)",code,fixed=TRUE)
eval(parse(text=code))
native_detect <- auto_detect_tissue_type
auto_detect_tissue_type <- function(...) {
  z <- native_detect(...)
  write.csv(z,file.path(out,"ScType_tissue_scores.csv"),row.names=FALSE);z
}
# Load pinned local upstream assets in lieu of sourcing floating master URLs.
sctype_source <- function() db
pdf(tempfile(fileext=".pdf"))
res <- run_sctype(small)
dev.off()
p <- unique(data.frame(cluster=as.character(res$seurat_clusters),ScType=res$sctype_classification))
stopifnot(nrow(p)==26)
write.csv(p,file.path(out,"ScType_native.csv"),row.names=FALSE)
