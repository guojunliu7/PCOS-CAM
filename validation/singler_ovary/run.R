# Run in PCOS-CAM root. Argument: full GSE268919 Seurat RDS or native input cache.
library(rhdf5)
library(Matrix)
library(SingleR)
library(celldex)
library(ggplot2)
source("Global.R")
out <- "validation/singler_ovary"
f <- file.path(out,"GSE267729.h5ad")
get_category <- function(path) {
  h5read(f,paste0(path,"/categories"))[as.integer(h5read(f,paste0(path,"/codes")))+1L]
}
mapping <- read.csv(file.path(out,"label_mapping.csv"),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(mapping[c("source_label","dictionary_label")]))
ref_source <- get_category(paste0("obs/",unique(mapping$source_field)))
labels <- mapping$reference_label[match(ref_source,mapping$source_label)]
stopifnot(!anyNA(labels),setequal(labels,mapping$reference_label))
donor <- get_category("obs/donor_id")
# Balanced, deterministic reference subsample to limit dense training intermediates.
set.seed(17021)
blocks <- split(which(!is.na(labels)),interaction(labels[!is.na(labels)],donor[!is.na(labels)],drop=TRUE))
idx <- sort(unlist(lapply(blocks,function(i) if(length(i)>50) sample(i,50) else i),use.names=FALSE))
read_csr <- function(prefix,nfeatures){
  ptr <- as.integer(h5read(f,paste0(prefix,"/indptr")))
  stopifnot(length(ptr)==length(labels)+1L)
  mat <- new("dgCMatrix",i=as.integer(h5read(f,paste0(prefix,"/indices"))),
      p=ptr,x=as.numeric(h5read(f,paste0(prefix,"/data"))),Dim=as.integer(c(nfeatures,length(labels))))
  mat[,idx,drop=FALSE]
}
cells <- read.csv("validation/method_comparison/sampled_cells.csv",colClasses="character")
input <- readRDS(commandArgs(TRUE)[1])
if(inherits(input,"Seurat")){
 test <- SeuratObject::LayerData(input,assay="RNA",layer="data")[,cells$cell,drop=FALSE]
 cl <- as.character(input[[]][cells$cell,"seurat_clusters"])
} else {test <- input$m; cl <- as.character(input$cl)}
rm(input)
stopifnot(identical(colnames(test),cells$cell),identical(cl,cells$cluster))
symbols <- h5read(f,"var/feature_name")
if(is.list(symbols)) symbols <- get_category("var/feature_name")
ref <- read_csr("X",length(symbols))
symbols <- as.character(symbols)
rownames(ref) <- symbols
barcodes <- as.character(h5read(f,"obs/_index"))[idx]
colnames(ref) <- barcodes
ref <- ref[which(!duplicated(symbols) & symbols %in% rownames(test)),,drop=FALSE]
stopifnot(!any(colnames(ref) %in% colnames(test)),length(unique(labels[idx]))==length(unique(mapping$reference_label)))
if(all(ref@x==round(ref@x))){
 lib <- Matrix::colSums(ref); stopifnot(all(lib>0))
 ref <- ref %*% Diagonal(x=1e4/lib)
 ref@x <- log1p(ref@x)
 colnames(ref) <- barcodes
}
d <- read_pcos_dictionary("data/File_1.xlsx")
inverse <- setNames(names(d$ortholog),d$ortholog)
marker_sets <- lapply(split(mapping$dictionary_label,mapping$reference_label),function(v){
 hits <- d$markers$marker_gene[d$markers$cell_cluster %in% v]
 intersect(unname(inverse[hits]),rownames(ref))
})
stopifnot(all(lengths(marker_sets)>0))
write.csv(data.frame(reference_label=rep(names(marker_sets),lengths(marker_sets)),mouse_gene=unlist(marker_sets,use.names=FALSE)),
          file.path(out,"label_markers.csv"),row.names=FALSE)
write.csv(data.frame(cell=barcodes,donor=donor[idx],source_label=ref_source[idx],reference_label=labels[idx]),file.path(out,"reference_cells.csv"),row.names=FALSE)
predict_calls <- function(r,y,genes="de"){
 p <- SingleR(test=test,ref=r,labels=y,genes=genes,clusters=cl)
 setNames(ifelse(is.na(p$pruned.labels),"Unknown",p$pruned.labels),rownames(p))
}
general <- MouseRNAseqData()
p0 <- predict_calls(general,general$label.fine)
p1 <- predict_calls(ref,labels[idx])
p2 <- predict_calls(ref,labels[idx],marker_sets)
baseline <- read.csv("validation/singler_dictionary/cluster_predictions.csv")
stopifnot(identical(unname(p0[as.character(baseline$cluster)]),baseline$native_label))
key <- unique(mapping[c("reference_label","family","specific")])
stopifnot(!anyDuplicated(key$reference_label))
make_rows <- function(pred,method,default=FALSE){
 z <- baseline[c("cluster","included","reference_major","reference_detail")]
 z$method <- method
 z$predicted <- unname(pred[as.character(z$cluster)])
 i <- match(z$predicted,key$reference_label)
 z$family <- if(default) baseline$native_family else ifelse(z$predicted=="Unknown","Unknown",key$family[i])
 z$specific <- if(default) !z$predicted %in% c("Unknown","Endothelial cells") else !is.na(i) & key$specific[i]
 stopifnot(!anyNA(z$family))
 z$major_match <- z$included & z$family==z$reference_major
 z$eligible_fine <- z$major_match & z$specific
 z
}
methods <- c("SingleR + PCOS−CAM","SingleR (ovary)","SingleR")
rows <- rbind(make_rows(p2,methods[1]),make_rows(p1,methods[2]),make_rows(p0,methods[3],TRUE))
write.csv(rows,file.path(out,"cluster_predictions.csv"),row.names=FALSE)
stats <- do.call(rbind,lapply(methods,function(m){
 z <- rows[rows$method==m,]
 data.frame(method=m,reviewed=sum(z$included),major_match=sum(z$major_match),
 distinct_fine=length(unique(z$predicted[z$eligible_fine])))
}))
write.csv(stats,file.path(out,"metrics.csv"),row.names=FALSE)
# Figure 5 displays the two requested end-to-end configurations.
methods <- c("SingleR + PCOS−CAM","SingleR")
stats <- stats[match(methods,stats$method),]
write.csv(stats,file.path(out,"Figure_5_plot_data.csv"),row.names=FALSE)
panels <- c("Broad-family matches\n(10 reviewed clusters)","Distinct fine labels\nin family-matched clusters")
df <- data.frame(method=factor(rep(methods,2),levels=methods),panel=factor(rep(panels,each=length(methods)),levels=panels),
 n=c(stats$major_match,stats$distinct_fine))
p <- ggplot(df,aes(method,n,fill=method))+geom_col(width=.60)+geom_text(aes(label=n),vjust=-.45,size=5.6,fontface="bold")+
 facet_wrap(~panel)+scale_fill_manual(values=c("#009E00","#FF8000"))+
 scale_y_continuous(limits=c(0,11),breaks=seq(0,10,2),expand=expansion(mult=c(0,.01)))+
 labs(x=NULL,y="Count",caption=paste("Test dataset (both groups): GSE268919.",
 "Reference: SingleR + PCOS−CAM, GSE267729; SingleR, celldex MouseRNAseqData.",
 "Groups differ in both reference data and marker selection.",sep="\n"))+
 theme_bw(base_size=14)+theme(plot.margin=margin(8,12,8,26),legend.position="none",panel.grid=element_blank(),panel.border=element_rect(colour="black",linewidth=1),
 strip.background=element_blank(),strip.text=element_text(face="bold"),axis.text=element_text(face="bold",colour="black"),
 axis.text.x=element_text(angle=20,hjust=1),axis.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
for(ext in c("png","pdf"))ggsave(file.path("../Figures_Tables",paste0("Figure_5.",ext)),p,width=11,height=5.2,dpi=320,bg="white",device=if(ext=="pdf") grDevices::cairo_pdf else "png")
print(stats)
