# Run from validation/native_comparison. Use native comparator predictions.
library(ggplot2)
source("../../hierarchy.R")
read <- function(path) read.csv(path, stringsAsFactors=FALSE)
pred <- read("../hierarchical_comparison/full_predictions.csv")
hier <- read("../../data/label_hierarchy.csv")
audit <- read("all_cluster_audit.csv")
six <- subset(pred, B=="Granulosa", select=c(cluster,B,A))
six$A_parent <- hier$parent[match(six$A,hier$cell_cluster)]
stopifnot(nrow(six)==6, !anyNA(six$A_parent))
six$candidate <- pcos_granulosa_candidate(six$B,six$A,setNames(hier$parent,hier$cell_cluster))
audit$method[audit$method=="PCOS-CAM B"] <- "PCOS-CAM"
methods <- c("PCOS-CAM","SingleR","ScType")
for (m in methods) {
  z <- subset(audit,method==m)
  stopifnot(!anyDuplicated(z$cluster),nrow(z)==26)
  idx <- match(six$cluster,z$cluster)
  if (m=="PCOS-CAM") {
    stopifnot(all(six$B==z$predicted[idx]))
    six$reviewed <- z$included[idx]
    six$reference <- z$reference_detail[idx]
  } else six[[m]] <- z$predicted[idx]
}
names(six)[match(c("B","A","A_parent"),names(six))] <- c("annotation","global_top_label","global_top_parent")
write.csv(six,"granulosa_six_cluster_review.csv",row.names=FALSE,na="")
counts <- do.call(rbind,lapply(methods,function(m) {
  z <- subset(audit,method==m)
  proposed <- rep(NA_character_,nrow(z))
  if(m=="PCOS-CAM") proposed <- six$candidate[match(z$cluster,six$cluster)]
  combined <- ifelse(z$specific,z$predicted,proposed)
  data.frame(method=m,reviewed=sum(z$included),major_match=sum(z$major_match),
    distinct_reported=length(unique(na.omit(z$predicted[z$eligible_fine]))),
    distinct_including_candidates=length(unique(na.omit(combined[z$major_match]))),
    clusters_reported=sum(z$specific),clusters_including_candidates=sum(!is.na(combined)))
}))
write.csv(counts,"granulosa_candidate_comparison.csv",row.names=FALSE)
panels <- c("Broad-family matches\n(10 reviewed clusters)","Distinct fine labels\nin family-matched clusters")
plot_data <- do.call(rbind,lapply(seq_along(panels),function(i) {
  total <- counts[[c("major_match","distinct_including_candidates")[i]]]
  data.frame(method=factor(counts$method,levels=methods),panel=factor(panels[i],levels=panels),
    total=total)
}))
g <- ggplot(plot_data,aes(method,fill=method))+
  geom_col(aes(y=total),width=.65)+
  geom_text(aes(y=total,label=total),vjust=-.5,fontface="bold",size=5)+
  facet_wrap(~panel,nrow=1,scales="free_y")+scale_fill_manual(values=setNames(c("#009E00","#FF8000","#9966B8"),methods))+
  scale_y_continuous(expand=expansion(mult=c(0,.22)))+
  labs(x=NULL,y="Count",caption="The PCOS-CAM total of 5 includes one candidate granulosa subtype.")+
  theme_bw(base_size=16)+theme(legend.position="none",panel.grid=element_blank(),strip.background=element_blank(),
  text=element_text(face="bold"),axis.text=element_text(colour="black",face="bold"),axis.text.x=element_text(angle=25,hjust=1),
  panel.border=element_rect(colour="black",fill=NA,linewidth=1),plot.caption=element_text(hjust=0,size=11))
for(ext in c("png","pdf")) ggsave(paste0("granulosa_candidate_comparison.",ext),g,width=9,height=5,dpi=300,bg="white")
