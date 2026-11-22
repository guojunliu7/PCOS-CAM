# Run from PCOS-CAM root. Arguments: annotated Seurat RDS, reference CSV,
# official ScType score R script, output folder. Requires UCell, ggplot2, SeuratObject.
if (nzchar(Sys.getenv("PCOS_BENCH_R_LIB"))) .libPaths(c(Sys.getenv("PCOS_BENCH_R_LIB"),.libPaths()))
args <- commandArgs(TRUE)
stopifnot(length(args)==4)
source("Global.R"); source("hierarchy.R"); source(args[3])
suppressPackageStartupMessages({library(Matrix);library(SeuratObject);library(ggplot2)})
out <- args[4];dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read_pcos_dictionary("data/File_1.xlsx")
h <- read_pcos_hierarchy("data/File_1.xlsx",d,"data/label_hierarchy.csv")
obj <- readRDS(args[1])
meta <- obj[[]]
set.seed(17021)
blocks <- split(rownames(meta),interaction(meta$sample_id,meta$seurat_clusters,drop=TRUE))
cells <- unlist(lapply(blocks,function(ids) if(length(ids)>40) sample(ids,40) else ids),use.names=FALSE)
cl <- as.character(meta[cells,"seurat_clusters"])
mat <- LayerData(obj,assay="RNA",layer="data")[,cells,drop=FALSE]
stopifnot(identical(colnames(mat),cells))
rm(obj);gc()
ids <- as.character(sort(as.integer(unique(cl))))
counts <- table(factor(cl,levels=ids))
mean_by_cluster <- function(m) {
  # cells by features -> clusters by features
  sweep(rowsum(as.matrix(m),factor(cl,levels=ids),reorder=FALSE)[ids,,drop=FALSE],1,as.numeric(counts),"/")
}
# All methods use identical cells. A/B consume linear-scale mean expression.
linear <- mat;linear@x <- expm1(linear@x)
x <- do.call(rbind,lapply(ids,function(i) data.frame(cluster=i,gene=rownames(mat),
                                     expression=Matrix::rowMeans(linear[,cl==i,drop=FALSE]))))
b <- pcos_annotate_b(x,d,h,"Mouse")$annotations
as_parent <- function(v) { ans <- unname(h$parent[v]);ifelse(is.na(ans),v,ans) }
pred <- data.frame(cluster=b$cluster,A=b$baseline_A,B=b$annotation)
gs <- split(d$markers$marker_gene,d$markers$cell_cluster)
inverse <- setNames(names(d$ortholog),d$ortholog)
gs <- lapply(gs,function(g) intersect(unname(inverse[g]),rownames(mat)))
write.csv(data.frame(label=names(gs),detected_mouse_markers=lengths(gs)),file.path(out,"signature_availability.csv"),row.names=FALSE)
gs <- gs[lengths(gs)>0]
stopifnot(identical(cells,read.csv("validation/method_comparison/sampled_cells.csv")$cell))
# ScType cluster-level sum with the author's score < nCells/4 rejection.
panel <- as.matrix(mat[unique(unlist(gs)),,drop=FALSE])
panel <- panel[apply(panel,1,sd)>0,,drop=FALSE]
sc <- sctype_score(panel,scaled=FALSE,gs=gs,gs2=NULL,gene_names_to_uppercase=FALSE)
scmeans <- mean_by_cluster(t(sc))
pick <- function(v,cutoff) {
  top <- max(v); k <- which(abs(v-top)<=32*.Machine$double.eps*abs(top))
  if(top<=0 || top<cutoff || length(k)!=1L) return("Unassigned")
  names(v)[k]
}
pred$ScType <- apply(scmeans,1,pick,cutoff=.25)[pred$cluster]
# UCell ranks genes against the whole transcriptome of each cell.
# Use neutral signature IDs so punctuation in a label cannot alter its name.
signature_ids <- setNames(names(gs),sprintf("sig%03d",seq_along(gs)))
names(gs) <- names(signature_ids)
uc <- UCell::ScoreSignatures_UCell(mat,features=gs,name=NULL,maxRank=1500,
                                   chunk.size=100,BPPARAM=BiocParallel::SerialParam())
stopifnot(setequal(rownames(uc),cells))
umeans <- mean_by_cluster(uc[cells,,drop=FALSE])
stopifnot(setequal(colnames(umeans),names(gs)))
colnames(umeans) <- unname(signature_ids[colnames(umeans)])
pred$UCell <- apply(umeans,1,pick,cutoff=0)[pred$cluster]

write.csv(pred,file.path(out,"full_predictions.csv"),row.names=FALSE)
