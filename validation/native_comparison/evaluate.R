library(ggplot2)
out <- "validation/native_comparison"
old <- read.csv("validation/hierarchical_comparison/cluster_audit.csv")
audit <- subset(old,method=="PCOS-CAM B",
                select=c(method,cluster,predicted,broad,reference_major,reference_detail,included,specific))
ref <- subset(audit,method=="PCOS-CAM B")
native <- function(method,table,label,families,coarse_labels) {
  z <- ref
  z$method <- method; z$predicted <- table[[label]][match(z$cluster,table$cluster)]
  z$broad <- unname(families[z$predicted])
  stopifnot(!anyNA(z$broad))
  z$specific <- !z$predicted %in% coarse_labels
  z
}
sc <- read.csv(file.path(out,"ScType_native.csv"));sr <- read.csv(file.path(out,"SingleR_native.csv"))
families <- setNames(rep("Neural / other",length(unique(sc$ScType))),unique(sc$ScType))
families[c("Fibroblast-like","Activated Microglia","Endothelial","Smooth Muscle","Pericytes","Unknown")] <-
  c("Stromal","Immune","Endothelial","Vascular mural","Vascular mural","Unknown")
audit <- rbind(audit,native("ScType",sc,"ScType",families,c("Unknown","Endothelial")),
  native("SingleR",sr,"SingleR",c("Fibroblasts"="Stromal","Fibroblasts activated"="Stromal",
    "Macrophages"="Immune","NK cells"="Immune","B cells"="Immune","Endothelial cells"="Endothelial",
    "Unknown"="Unknown"),c("Unknown","Endothelial cells")))
audit$major_match <- with(audit,included & reference_major==broad)
audit$mixed <- with(audit,included & broad=="Theca / stroma" & reference_major %in% c("Theca","Stromal"))
audit$eligible_fine <- with(audit,major_match & specific)
# Keyed to cluster and concept, identically for every caller. A qualifier added
# beyond the reference receives partial support, not specific-identity support.
audit$support <- ifelse(audit$included,"No supported finer identity","Outside review")
supported <- with(audit,(cluster==2 & predicted %in% c("Fibroblasts","Fibroblast-like","Fibroblast-like stroma")) |
  (cluster==5 & predicted=="Steroidogenic theca") |
  (cluster==11 & predicted=="Macrophages") | (cluster==20 & predicted %in% c("B cells","B cell")))
partial <- with(audit,(cluster==16 & predicted %in% c("NK cells","NK cell")) |
  (cluster==1 & predicted=="GC_Mural-(Preantral) (Mantri)"))
audit$support[audit$eligible_fine & supported] <- "Supported identity"
audit$support[audit$eligible_fine & partial] <- "Partially supported"
write.csv(audit,file.path(out,"all_cluster_audit.csv"),row.names=FALSE,na="")
eligible <- subset(audit,eligible_fine)
write.csv(eligible,file.path(out,"fine_calls_in_correct_family.csv"),row.names=FALSE)
methods <- unique(audit$method)
stats <- do.call(rbind,lapply(methods,function(m) {
  x <- subset(audit,method==m & included)
  count <- function(k) length(unique(x$predicted[k]))
  data.frame(method=m,reviewed=nrow(x),major_match=sum(x$major_match),mixed=sum(x$mixed),
    fine_clusters=sum(x$eligible_fine),distinct_fine=count(x$eligible_fine),
    supported=count(x$support=="Supported identity"),partial_only=count(x$support=="Partially supported"))
}))
print(stats);write.csv(stats,file.path(out,"metrics.csv"),row.names=FALSE)
source(file.path(out,"plot_ternary.R"))
plot_native_ternary(stats,out)
