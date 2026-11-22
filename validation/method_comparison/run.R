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
pred <- data.frame(cluster=b$cluster,A=as_parent(b$baseline_A),B=ifelse(b$level=="fine",as_parent(b$annotation),b$annotation))
gs <- split(h$markers$marker_gene,h$markers$cell_cluster)
inverse <- setNames(names(d$ortholog),d$ortholog)
gs <- lapply(gs,function(g) intersect(unname(inverse[g]),rownames(mat)))
stopifnot(all(lengths(gs)>0))
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
uc <- UCell::ScoreSignatures_UCell(mat,features=gs,name=NULL,maxRank=1500,
                                   chunk.size=100,BPPARAM=BiocParallel::SerialParam())
stopifnot(setequal(rownames(uc),cells))
umeans <- mean_by_cluster(uc[cells,,drop=FALSE])
stopifnot(setequal(colnames(umeans),names(gs)))
pred$UCell <- apply(umeans,1,pick,cutoff=0)[pred$cluster]
ref <- read.csv(args[2],colClasses="character")
label_map <- c(Granulosa="Granulosa","Stromal/fibroblast"="Stromal","Theca/interstitial"="Theca",
 "Proliferating granulosa"="Granulosa",Endothelial="Endothelial",Immune="Immune","Proliferating stromal"="Stromal")
j <- match(pred$cluster,ref$seurat_clusters);stopifnot(!anyNA(j))
pred$n_cells <- unname(as.integer(counts[pred$cluster]))
pred$reference_original <- ref$refined_major_cell_type[j]
pred$reference <- unname(label_map[pred$reference_original])
pred$included <- ref$annotation_confidence[j]=="high" & !is.na(pred$reference)
write.csv(pred,file.path(out,"cluster_predictions.csv"),row.names=FALSE,na="")
write.csv(data.frame(cell=cells,sample=meta[cells,"sample_id"],cluster=cl),file.path(out,"sampled_cells.csv"),row.names=FALSE)
test <- pred[pred$included,];truth <- test$reference
methods <- c("A","B","ScType","UCell")
display <- c("PCOS-CAM A","PCOS-CAM B","ScType","UCell + top label")
detail <- list()
metrics <- do.call(rbind,lapply(seq_along(methods),function(i){
  p <- test[[methods[i]]]; accepted <- !p %in% c("Unassigned","Ambiguous")
  cls <- sort(unique(truth));tp <- vapply(cls,function(k) sum(p==k & truth==k),integer(1))
  n <- vapply(cls,function(k)sum(truth==k),integer(1));calls <- vapply(cls,function(k)sum(p==k),integer(1))
  f <- 2*tp/(n+calls)
  detail[[i]] <<- data.frame(method=display[i],type=cls,n=n,recall=tp/n,F1=f)
  data.frame(method=display[i],n=length(truth),matched=sum(p==truth),called=sum(accepted),
             agreement=mean(p==truth),macro_F1=mean(f),coverage=mean(accepted),
             agreement_called=if(any(accepted))mean(p[accepted]==truth[accepted]) else NA_real_)
}))
write.csv(metrics,file.path(out,"metrics.csv"),row.names=FALSE)
write.csv(do.call(rbind,detail),file.path(out,"class_metrics.csv"),row.names=FALSE)
long <- reshape2::melt(metrics,id.vars="method",measure.vars=c("agreement","macro_F1","coverage"))
long$variable <- factor(long$variable,levels=c("agreement","macro_F1","coverage"),
  labels=c("Reference agreement","Macro-F1","Coverage"))
long$method <- factor(long$method,levels=rev(display))
p <- ggplot(long,aes(value,method,color=method))+geom_segment(aes(x=0,xend=value,yend=method),linewidth=1)+
  geom_point(size=3)+geom_text(aes(label=sprintf("%.2f",value)),nudge_x=.055,size=3.4,show.legend=FALSE)+
  facet_grid(.~variable)+scale_x_continuous(limits=c(0,1.16),breaks=c(0,.5,1))+
  scale_color_manual(values=c("#687584","#c17a28","#057267","#374a69"))+
  labs(title="Mouse ovarian cell-cluster annotation",
       subtitle=sprintf("GSE268919  |  %d reference clusters, %d represented classes  |  %s sampled cells across all 26 clusters",nrow(test),length(unique(truth)),format(length(cells),big.mark=",")),
       x=NULL,y=NULL,caption="Local marker-based reference. Abstentions count as mismatches. UCell uses the highest mean signature score.")+
  theme_minimal(base_size=11)+theme(legend.position="none",panel.grid.major.y=element_blank(),
  panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),plot.title=element_text(face="bold",size=17))
for(ext in c("pdf","png")) ggsave(file.path(out,paste0("comparison.",ext)),p,width=10,height=3.8,dpi=220)
print(metrics)
