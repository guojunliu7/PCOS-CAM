# Run from PCOS-CAM root. Pass the full GSE268919 Seurat RDS or the saved native input.
library(SingleR)
library(celldex)
library(ggplot2)
source("Global.R")
out <- "validation/singler_dictionary"
cells <- read.csv("validation/method_comparison/sampled_cells.csv",colClasses="character")
obj <- readRDS(commandArgs(TRUE)[1])
if(inherits(obj,"Seurat")){
  x <- SeuratObject::LayerData(obj,assay="RNA",layer="data")[,cells$cell,drop=FALSE]
  cl <- as.character(obj[[]][cells$cell,"seurat_clusters"])
} else {
  x <- obj$m
  cl <- as.character(obj$cl)
}
rm(obj)
stopifnot(identical(colnames(x),cells$cell),identical(cl,cells$cluster),!anyDuplicated(cells$cell))
ref <- MouseRNAseqData()
d <- read_pcos_dictionary("data/File_1.xlsx")
genes <- intersect(intersect(names(d$ortholog)[d$ortholog %in% d$markers$marker_gene],rownames(x)),rownames(ref))
stopifnot(length(genes)>1)
write.csv(data.frame(mouse_gene=genes,human_gene=unname(d$ortholog[genes])),file.path(out,"marker_genes.csv"),row.names=FALSE)
get_calls <- function(features=NULL){
  p <- SingleR(test=x,ref=ref,labels=ref$label.fine,clusters=cl,restrict=features)
  data.frame(cluster=rownames(p),label=ifelse(is.na(p$pruned.labels),"Unknown",p$pruned.labels),unpruned=p$labels)
}
native <- get_calls()
original <- read.csv("validation/native_comparison/SingleR_native.csv")
stopifnot(identical(native$label,original$SingleR[match(native$cluster,original$cluster)]))
custom <- get_calls(genes)
refmap <- unique(data.frame(label=ref$label.fine,main=ref$label.main))
stopifnot(!anyDuplicated(refmap$label))
audit <- subset(read.csv("validation/native_comparison/all_cluster_audit.csv"),method=="SingleR")
z <- audit[order(audit$cluster),c("cluster","included","reference_major","reference_detail")]
for (arm in c("native","custom")){
  pred <- get(arm)
  labels <- pred$label[match(z$cluster,pred$cluster)]
  z[[paste0(arm,"_label")]] <- labels
  z[[paste0(arm,"_main")]] <- refmap$main[match(labels,refmap$label)]
  z[[paste0(arm,"_unpruned")]] <- pred$unpruned[match(z$cluster,pred$cluster)]
}

family <- c("Fibroblasts"="Stromal","Macrophages"="Immune","NK cells"="Immune","B cells"="Immune",
 "Endothelial cells"="Endothelial","Neurons"="Neural","Adipocytes"="Adipocyte")
methods <- c("SingleR + PCOS−CAM","SingleR")
counts <- list()
for(arm in c("custom","native")){
  label <- z[[paste0(arm,"_label")]]
  broad <- unname(family[z[[paste0(arm,"_main")]]])
  broad[label=="Unknown"] <- "Unknown"
  stopifnot(!anyNA(broad))
  match_major <- z$included & broad==z$reference_major
  counted <- match_major & !label %in% c("Endothelial cells","Unknown")
  z[[paste0(arm,"_family")]] <- broad
  z[[paste0(arm,"_match")]] <- match_major
  z[[paste0(arm,"_eligible_fine")]] <- counted
  counts[[arm]] <- c(reviewed=sum(z$included),major_match=sum(match_major),distinct_fine=length(unique(label[counted])))
}
stats <- data.frame(method=methods,do.call(rbind,counts),row.names=NULL)
write.csv(stats,file.path(out,"metrics.csv"),row.names=FALSE)
write.csv(z,file.path(out,"cluster_predictions.csv"),row.names=FALSE)
panels <- c("Broad-family matches\n(10 reviewed clusters)","Distinct fine labels\nin family-matched clusters")
df <- data.frame(method=factor(rep(methods,2),levels=methods),panel=factor(rep(panels,each=2),levels=panels),
 n=c(stats$major_match,stats$distinct_fine))
p <- ggplot(df,aes(method,n,fill=method))+geom_col(width=.56)+geom_text(aes(label=n),vjust=-.45,size=5.6,fontface="bold")+
 facet_wrap(~panel)+scale_fill_manual(values=c("#009E00","#FF8000"))+
 scale_y_continuous(limits=c(0,11),breaks=seq(0,10,2),expand=expansion(mult=c(0,.01)))+
 labs(x=NULL,y="Count",caption=sprintf("Both arms use MouseRNAseqData. Green: marker selection restricted to %d PCOS-CAM genes.",length(genes)))+
 theme_bw(base_size=14)+theme(legend.position="none",panel.grid=element_blank(),panel.border=element_rect(colour="black",linewidth=1),
 strip.background=element_blank(),strip.text=element_text(face="bold"),axis.text=element_text(face="bold",colour="black"),
 axis.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
for(ext in c("png","pdf")) ggsave(file.path("../Figures_Tables",paste0("Figure_4.",ext)),p,width=9,height=4.8,dpi=320,bg="white",
 device=if(ext=="pdf") grDevices::cairo_pdf else "png")
print(stats)
