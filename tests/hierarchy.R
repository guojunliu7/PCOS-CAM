# Run from the project root with Rscript tests/hierarchy.R.
source("Global.R")
source("hierarchy.R")
path <- "data/File_1.xlsx"
d <- read_pcos_dictionary(path)
h <- read_pcos_hierarchy(path,d,"data/label_hierarchy.csv")
stopifnot(length(h$parent)==153,all(abs(tapply(h$markers$percent,h$markers$cell_cluster,sum)-1)<1e-12))
e <- as.data.frame(readxl::read_excel(path,"Evidence",col_types="text"))
for (i in seq_len(nrow(h$markers))) {
  m <- h$markers[i,]
  ok <- !is.na(e$integrated_gene) & e$integrated_gene==m$marker_gene & h$parent[e$cell_cluster]==m$cell_cluster
  stopifnot(length(unique(e$vote_unit[ok]))==m$Number_of_votes_received)
}
# Identical subtype signatures force fallback. Removing a leading gene must also
# remove it from the competitors, with ties among leading genes tested separately.
toy <- function(labels,genes,weights) data.frame(cell_cluster=labels,marker_gene=genes,percent=weights)
fine <- list(markers=toy(rep(c("a1","a2","b1"),each=2),c("X","Y","X","Y","Q","R"),rep(.5,6)),
             ortholog=c(x="X",y="Y"))
hier <- list(parent=c(a1="A",a2="A",b1="B"),markers=toy(rep(c("A","B"),each=2),c("X","Y","Q","R"),rep(.5,4)))
x <- data.frame(cluster=c("c","c","single","none"),gene=c("x","y","x","unknown"),expression=1)
r <- pcos_annotate_b(x,fine,hier,"Mouse")
stopifnot(identical(r$annotations$annotation,c("A","Unassigned","Unassigned")),
          r$diagnostics$stable[1],r$diagnostics$matched[1]==2,
          r$diagnostics$leading_gene[1]=="X | Y",r$diagnostics$weight_coverage[1]==1,
          r$annotations$reason[1]=="within_parent_tie")
hits <- data.frame(cell_cluster=c("a","a","b","b"), marker_gene=c("X","Y","X","Z"),
                   percent=.5,contribution=c(.5,.1,.2,.3))
res <- pcos_evidence_summary(c(a=.6,b=.5),hits)
stopifnot(!res$stable,res$after_removal=="b")
