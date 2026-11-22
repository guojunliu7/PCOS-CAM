# Run from PCOS-CAM root. Pass the full GSE268919 Seurat RDS or the saved native input.
library(SingleR)
library(celldex)
source("Global.R")
out <- "validation/singler_fixed_reference_20261122"
set.seed(17021)
writeLines(capture.output(sessionInfo()), file.path(out,"sessionInfo.txt"))
inputs <- c(query=commandArgs(TRUE)[1], dictionary="../Figures_Tables/File_1.xlsx", cell_manifest="validation/method_comparison/sampled_cells.csv", evaluation="validation/native_comparison/all_cluster_audit.csv", native_predictions="validation/native_comparison/SingleR_native.csv")
write.csv(data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs))),file.path(out,"input_manifest.csv"),row.names=FALSE)
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
d <- read_pcos_dictionary("../Figures_Tables/File_1.xlsx")
frozen <- read_pcos_dictionary("data/File_1.xlsx")
stopifnot(identical(d$markers, frozen$markers),identical(d$ortholog, frozen$ortholog))
write.csv(data.frame(reference_sample=colnames(ref),label_fine=ref$label.fine,label_main=ref$label.main),file.path(out,"reference_sample_labels.csv"),row.names=FALSE)
writeLines(c(paste("Query cells",ncol(x)),paste("Query genes",nrow(x)),paste("Query clusters",length(unique(cl))),paste("Reference",nrow(ref),"genes",ncol(ref),"samples"),paste("Reference fine labels",length(unique(ref$label.fine)))),file.path(out,"run_configuration.txt"))
genes <- intersect(intersect(names(d$ortholog)[d$ortholog %in% d$markers$marker_gene],rownames(x)),rownames(ref))
stopifnot(length(genes)>1)
write.csv(data.frame(mouse_gene=genes,human_gene=unname(d$ortholog[genes])),file.path(out,"marker_genes.csv"),row.names=FALSE)
get_calls <- function(features=NULL){
  p <- SingleR(test=x,ref=ref,labels=ref$label.fine,clusters=cl,restrict=features)
  saveRDS(p,file.path(out,if(is.null(features)) "SingleR_default_scores.rds" else "SingleR_dictionary_scores.rds"))
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

review <- z[z$included,]
review$label_changed <- review$custom_label != review$native_label
review$family_changed <- review$custom_family != review$native_family
write.csv(review,file.path(out,"reviewed_cluster_comparison.csv"),row.names=FALSE)
all_changes <- z[z$custom_label!=z$native_label,]
write.csv(all_changes,file.path(out,"changed_clusters.csv"),row.names=FALSE)
counts_fine <- do.call(rbind,lapply(c("custom","native"),function(arm){
 lab <- z[[paste0(arm,"_label")]][z[[paste0(arm,"_eligible_fine")]]]
 a<-as.data.frame(table(lab),stringsAsFactors=FALSE)
 names(a)<-c("label","n_clusters")
 a$method <- if(arm=="custom") "SingleR + PCOS-CAM" else "SingleR"
 a[,c("method","label","n_clusters")]
}))
write.csv(counts_fine,file.path(out,"counted_fine_labels.csv"),row.names=FALSE)
writeLines(c(paste("PCOS-CAM eligible genes",length(genes)),paste("All-cluster label changes",sum(z$custom_label!=z$native_label)),paste("Reviewed-cluster label changes",sum(review$label_changed)),paste("Reviewed-cluster family changes",sum(review$family_changed)),"Default predictions reproduce the 26 saved native calls.","No figures or manuscript files were modified."),file.path(out,"verification.txt"))
print(stats)
print(review[,c("cluster","reference_major","native_label","custom_label","native_match","custom_match")])
print(counts_fine)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
