# Exploratory agreement with local marker-based annotations, not independent truth.
# Usage (from project root): Rscript validation/compare_GSE268919.R INPUT_DIR OUTPUT.csv
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
source("Global.R")
source("hierarchy.R")
path <- "data/File_1.xlsx"
d <- read_pcos_dictionary(path)
h <- read_pcos_hierarchy(path,d,"data/label_hierarchy.csv")
x <- read_pcos_expression(file.path(args[1],"average_marker_expression_by_cluster.csv"))
b <- pcos_annotate_b(x,d,h,"Mouse")$annotations
ref <- read.csv(file.path(args[1],"refined_cluster_annotation.csv"))
label_map <- c(Granulosa="Granulosa","Stromal/fibroblast"="Stromal","Theca/interstitial"="Theca",
  "Proliferating granulosa"="Granulosa",Endothelial="Endothelial",Immune="Immune","Proliferating stromal"="Stromal")
ix <- match(b$cluster,paste0("g",ref$seurat_clusters))
stopifnot(!anyNA(ix),!anyDuplicated(ix),length(ix)==nrow(ref))
b$reference <- ref$refined_major_cell_type[ix]
b$reference_confidence <- ref$annotation_confidence[ix]
b$reference_parent <- unname(label_map[b$reference])
b$included <- b$reference_confidence=="high" & !is.na(b$reference_parent)
b$A_parent <- unname(h$parent[b$baseline_A])
b$B_parent <- ifelse(b$level=="fine",h$parent[b$annotation],b$annotation)
write.csv(b,args[2],row.names=FALSE,na="")
eval <- b[b$included,]
for (field in c("A_parent","B_parent")) {
  pred <- eval[[field]];truth <- eval$reference_parent
  f1 <- vapply(unique(truth),function(k) 2*sum(pred==k & truth==k,na.rm=TRUE)/
    (sum(pred==k,na.rm=TRUE)+sum(truth==k)),numeric(1))
  print(c(matches=sum(pred==truth,na.rm=TRUE),total=length(truth),macro_F1=mean(f1)))
}
