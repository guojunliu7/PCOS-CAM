# From PCOS-CAM root, after run_fine.R. One output per method and cluster.
library(ggplot2)
out <- "validation/hierarchical_comparison"
original <- read.csv("validation/method_comparison/cluster_predictions.csv",colClasses="character")
pred <- read.csv(file.path(out,"full_predictions.csv"),colClasses="character")
map <- read.csv("data/label_hierarchy.csv")
parent <- setNames(map$parent,map$cell_cluster)
methods <- c("A","B","ScType","UCell")
method_names <- c("PCOS-CAM A","PCOS-CAM B","ScType","UCell + top label")
coarse <- function(x) {
  p <- if (x %in% names(parent)) unname(parent[x]) else x
  if (p=="Immune-like (CD45-)") return("Immune")
  if (p=="Mesenchymal") return("Stromal")
  p
}
specific <- function(x) {
  if (x %in% c("Unassigned","Ambiguous","Mesenchymal","Theca / stroma")) return(FALSE)
  # A finer named phenotype can itself be a B parent fallback.
  if (x=="Immune-like (CD45-)") return(TRUE)
  x %in% names(parent) && unname(parent[x])!=x
}
ref_path <- "../../PCOS_JQF/GSE138518/GSE268919_scRNA_refined_annotation/refined_cluster_annotation.csv"
refs <- read.csv(ref_path,colClasses="character")
# Explicit concordance of at least one reference subtype feature (not confirmation of the whole fine label).
patterns <- c("1"="preantral","2"="fibroblast","5"="androgenic|steroidogenic theca",
              "7"="(^|[^[:alpha:]])antral","8"="mitotic|cycling|proliferat|TOP2A|UBE2C","9"="blood vascular",
              "11"="macrophage","13"="mitotic|cycling|proliferat",
              "16"="^T / NK$|^T cell$|^NK cell$|^CD4\\+ T","20"="^B cell$")
# Cluster 9 is endothelial but its reference does not explicitly establish blood-vs-lymphatic subtype.
patterns["9"] <- "(?!)"
# Neither CD4 nor CD8 status is specified by the T/NK reference.
patterns["16"] <- "^T / NK$|^T cell$|^NK cell$"
stopifnot(setequal(pred$cluster,original$cluster),!anyDuplicated(pred$cluster),
          setequal(names(patterns),original$cluster[original$included=="TRUE"]))
rows <- do.call(rbind,lapply(seq_along(methods),function(i) {
  lab <- pred[[methods[i]]]; broad <- vapply(lab,coarse,character(1),USE.NAMES=FALSE)
  idx <- match(pred$cluster,original$cluster); incl <- original$included[idx]=="TRUE"
  ref <- original$reference[idx]; fine <- vapply(lab,specific,logical(1))
  agree <- incl & broad==ref
  mixed <- incl & broad=="Theca / stroma" & ref %in% c("Stromal","Theca")
  feature <- vapply(seq_along(lab),function(j) {
    if(!agree[j] || !fine[j]) return(FALSE)
    grepl(patterns[[pred$cluster[j]]],lab[j],ignore.case=TRUE,perl=TRUE)
  },logical(1))
  data.frame(method=method_names[i],cluster=pred$cluster,predicted=lab,
             broad=broad,reference_major=ref,
             reference_detail=refs$refined_sub_cell_type[match(pred$cluster,refs$seurat_clusters)],
             included=incl,specific=fine,major_agree=ifelse(incl,agree,NA),
             mixed_compatible=ifelse(incl,mixed,NA),detail_feature_supported=ifelse(incl,feature,NA))
}))
write.csv(rows,file.path(out,"cluster_audit.csv"),row.names=FALSE,na="")
# Same first-stage benchmark as the previous comparison, rescored for broad compatibility only.
stage1 <- do.call(rbind,lapply(seq_along(methods),function(i) {
 x <- original[original$included=="TRUE",];v <- vapply(x[[methods[i]]],coarse,character(1))
 data.frame(method=method_names[i],n=nrow(x),major_match=sum(v==x$reference),
 mixed=sum(v=="Theca / stroma" & x$reference %in% c("Theca","Stromal")))
}))
write.csv(stage1,file.path(out,"previous_parent_rescored.csv"),row.names=FALSE)
summary <- do.call(rbind,lapply(method_names,function(m) {
 x <- subset(rows,method==m);z <- subset(x,included)
 data.frame(method=m,major=sum(z$major_agree),mixed=sum(z$mixed_compatible),compatible=sum(z$major_agree | z$mixed_compatible),
  finer=sum(z$major_agree & z$specific),feature=sum(z$detail_feature_supported),
  unique_correct_major=length(unique(z$predicted[z$major_agree & z$specific])),
  specific_clusters=sum(x$specific),unique_subtypes=length(unique(x$predicted[x$specific])))
}))
write.csv(summary,file.path(out,"metrics.csv"),row.names=FALSE)
panels <- c("Major match\n/ 10 reviewed","Including mixed match\n/ 10 reviewed","More specific, major match\n/ 10 reviewed",
            "Reference feature supported\n/ 10 reviewed","Distinct subtype labels\n/ all 26 clusters")
long <- do.call(rbind,lapply(seq_along(panels),function(i)
 data.frame(method=summary$method,panel=panels[i],n=summary[[c("major","compatible","finer","feature","unique_subtypes")[i]]])))
long$panel <- factor(long$panel,levels=panels);long$method <- factor(long$method,levels=rev(method_names))
p <- ggplot(long,aes(n,method,color=method))+geom_segment(aes(x=0,xend=n,yend=method),linewidth=1)+
 geom_point(size=3)+geom_text(aes(label=n),nudge_x=.55,size=3.7,show.legend=FALSE)+
 facet_grid(.~panel,scales="free_x",space="fixed")+expand_limits(x=0)+
 scale_x_continuous(expand=expansion(mult=c(.04,.22)),breaks=scales::breaks_pretty(4))+
 scale_color_manual(values=c("#687584","#bc7326","#08746b","#374a69"))+
 labs(title="Ovarian cluster annotation at two levels",subtitle="GSE268919 mouse ovary   |   4,913 cells   |   Shared fine-label dictionary",
 x=NULL,y=NULL,caption="Mixed theca/stroma calls are tabulated separately. Subtype diversity counts distinct reported labels, irrespective of correctness.")+
 theme_minimal(base_size=11)+theme(legend.position="none",panel.grid.major.y=element_blank(),
 panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),plot.title=element_text(face="bold",size=17))
for (ext in c("png","pdf")) ggsave(file.path(out,paste0("hierarchical_comparison.",ext)),p,width=14,height=4.0,dpi=220)
print(stage1);print(summary)
audit <- subset(rows,included)
audit$status <- with(audit,ifelse(detail_feature_supported,"Reference feature",
 ifelse(major_agree & specific,"More specific",ifelse(major_agree,"Major match",
 ifelse(mixed_compatible,"Mixed compatible",ifelse(predicted %in% c("Ambiguous","Unassigned"),"Unresolved","Different major"))))))
audit$method <- factor(audit$method,levels=method_names)
audit$key <- paste0("g",audit$cluster,"  ",audit$reference_major)
audit$key <- factor(audit$key,levels=rev(unique(audit$key[order(as.integer(audit$cluster))])))
audit$text <- vapply(audit$predicted,function(v)paste(strwrap(v,width=25),collapse="\n"),character(1))
palette <- c("Reference feature"="#88c6a7","More specific"="#bcded7","Major match"="#dee9e9",
             "Mixed compatible"="#f5dda5","Unresolved"="#eceef0","Different major"="#efc5c1")
p2 <- ggplot(audit,aes(method,key))+geom_tile(aes(fill=status),color="white",linewidth=1)+
 geom_text(aes(label=text),size=3.2)+scale_fill_manual(values=palette,breaks=names(palette))+
 scale_x_discrete(position="top")+labs(title="Which finer labels were actually assigned?",subtitle="GSE268919   |   10 reviewed clusters   |   Row labels show the local reference family",
 x=NULL,y=NULL,fill=NULL)+theme_minimal(base_size=11)+theme(panel.grid=element_blank(),legend.position="bottom",
 plot.title=element_text(face="bold",size=17),axis.text.x=element_text(face="bold"))
for (ext in c("pdf","png")) ggsave(file.path(out,paste0("cluster_label_audit.",ext)),p2,width=12,height=7,dpi=180)
