source("Global.R")
d <- read_pcos_dictionary("data/File_1.xlsx")
m <- merge(d$markers,read.csv("data/label_hierarchy.csv"))
labels <- c("Endothelial cells","Epithelial cells","Luteal cells","M/B cells","NK/T cells","Oocytes","Stromal cells","Theca cells")
dict <- list("Endothelial","Epithelial","Luteal",c("Macrophage","Monocyte","B cell"),"T / NK","Oocyte","Stromal","Theca")
families <- c("Endothelial","Epithelial","Luteal","Immune","Immune","Oocyte","Stromal","Theca")
rows <- lapply(seq_along(labels),function(i)data.frame(source_label=labels[i],dictionary_label=dict[[i]],family=families[i],specific=families[i]=="Immune"))
for(g in c("Amh","Inhbb","Nppc")){
 labs <- unique(m$cell_cluster[m$parent=="Granulosa" & m$marker_gene==toupper(g)])
 stopifnot(length(labs)>0)
 rows[[length(rows)+1L]] <- data.frame(source_label=paste0("Granulosa cells (",g,"-high)"),dictionary_label=labs,family="Granulosa",specific=TRUE)
}
mapping <- do.call(rbind,rows)
mapping$source_field <- "author_cell_type"
mapping$reference_label <- mapping$source_label
write.csv(mapping,"validation/singler_ovary/label_mapping.csv",row.names=FALSE)
