# Run from the PCOS-CAM root with Rscript validation/sctype_dictionary/compare.R <ScType checkout> <Seurat RDS>
library(SeuratObject)
library(ggplot2)
library(HGNChelper)
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
source(file.path(args[1],"R/sctype_score_.R"))
source(file.path(args[1],"R/gene_sets_prepare.R"))
source("Global.R")
out <- "validation/sctype_dictionary"
native <- "validation/native_comparison"
cells <- read.csv("validation/method_comparison/sampled_cells.csv",colClasses="character")
obj <- readRDS(args[2])
stopifnot(!anyDuplicated(cells$cell),identical(as.character(obj[[]][cells$cell,"seurat_clusters"]),cells$cluster))
x <- LayerData(obj,assay="RNA",layer="scale.data")[,cells$cell,drop=FALSE]
stopifnot(identical(colnames(x),cells$cell))
rm(obj);gc()
cl <- factor(cells$cluster,levels=sort(unique(as.integer(cells$cluster))))
d <- read_pcos_dictionary("data/File_1.xlsx")
inv <- setNames(names(d$ortholog),d$ortholog)
gs <- lapply(split(d$markers$marker_gene,d$markers$cell_cluster),function(g) unique(toupper(na.omit(unname(inv[g])))))
tissue <- read.csv(file.path(native,"ScType_tissue_scores.csv"))$tissue[1]
official <- gene_sets_prepare(file.path(args[1],"ScTypeDB_full.xlsx"),tissue)
score <- function(pos,neg){
  keep <- toupper(rownames(x)) %in% c(unlist(pos),unlist(neg))
  es <- sctype_score(as.matrix(x[keep,,drop=FALSE]),scaled=TRUE,gs=pos,gs2=neg)
  sapply(levels(cl),function(id){
    sums <- sort(rowSums(es[,cl==id,drop=FALSE]),decreasing=TRUE)
    if(sums[1]<sum(cl==id)/4) "Unknown" else names(sums)[1]
  })
}
default <- score(official$gs_positive,official$gs_negative)
old <- read.csv(file.path(native,"ScType_native.csv"))
stopifnot(identical(unname(default),old$ScType[match(names(default),old$cluster)]))
custom <- score(gs,NULL)
map <- read.csv("data/label_hierarchy.csv")
parent <- setNames(map$parent,map$cell_cluster)
audit <- subset(read.csv(file.path(native,"all_cluster_audit.csv")),method=="ScType")
a <- audit[match(names(custom),audit$cluster),]
stopifnot(!anyNA(a$cluster),sum(a$included)==10)
family <- unname(parent[custom])
family[custom=="Unknown"] <- "Unknown"
specific <- !custom %in% c("Unknown","Mesenchymal","Theca / stroma") & (family!=custom | custom=="Immune-like (CD45-)")
family[family=="Mesenchymal"] <- "Stromal"
family[family=="Immune-like (CD45-)"] <- "Immune"
family[family %in% c("Smooth muscle","Perivascular")] <- "Vascular mural"
stopifnot(!anyNA(family))
z <- data.frame(cluster=a$cluster,reviewed=a$included,reference_major=a$reference_major,reference_detail=a$reference_detail,
 custom_label=unname(custom),custom_family=family,custom_specific=unname(specific),
 default_label=unname(default),default_family=a$broad,default_specific=a$specific)
# Use vascular mural as the common parent for smooth muscle and pericytes in both arms.
z$custom_specific[z$custom_label %in% c("Smooth muscle","Pericyte")] <- TRUE
z$custom_match <- z$reviewed & z$custom_family==z$reference_major
z$default_match <- z$reviewed & z$default_family==z$reference_major
write.csv(z,file.path(out,"cluster_predictions.csv"),row.names=FALSE)
methods <- c("ScType + PCOS−CAM","ScType")
stats <- data.frame(method=methods,reviewed=sum(z$reviewed),major_match=c(sum(z$custom_match),sum(z$default_match)),
 distinct_fine=c(length(unique(z$custom_label[z$custom_match & z$custom_specific])),
 length(unique(z$default_label[z$default_match & z$default_specific]))))
write.csv(stats,file.path(out,"metrics.csv"),row.names=FALSE)
panels <- c("Broad-family matches\n(10 reviewed clusters)","Distinct fine labels\nin family-matched clusters")
plotdata <- data.frame(method=factor(rep(methods,2),levels=methods),panel=factor(rep(panels,each=2),levels=panels),
 n=c(stats$major_match,stats$distinct_fine))
p <- ggplot(plotdata,aes(method,n,fill=method))+geom_col(width=.56)+
 geom_text(aes(label=n),vjust=-.45,size=5.6,fontface="bold")+
 facet_wrap(~panel)+scale_fill_manual(values=c("#009E00","#9966B8"))+
 scale_y_continuous(limits=c(0,11),breaks=seq(0,10,2),expand=expansion(mult=c(0,.01)))+
 labs(x=NULL,y="Count")+theme_bw(base_size=14)+
 theme(legend.position="none",panel.grid=element_blank(),panel.border=element_rect(colour="black",linewidth=1),
 strip.background=element_blank(),strip.text=element_text(face="bold"),axis.text=element_text(face="bold",colour="black"),
 axis.title=element_text(face="bold"))
for(ext in c("png","pdf")) ggsave(file.path("../Figures_Tables",paste0("Figure_3.",ext)),p,width=9,height=4.8,dpi=320,bg="white",
 device=if(ext=="pdf") grDevices::cairo_pdf else "png")
print(stats)
