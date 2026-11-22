library(ggplot2)
x <- read.csv("metrics.csv")
methods <- c("PCOS-CAM B", "SingleR", "ScType")
stopifnot(setequal(methods,x$method),all(x$reviewed==10))
x$method <- factor(x$method,levels=methods)
audit <- read.csv("all_cluster_audit.csv")
stopifnot(all(table(audit$method)==26), !anyDuplicated(audit[c("method","cluster")]))
totals <- tapply(audit$specific, audit$method, sum)
x$all_fine_calls <- unname(totals[as.character(x$method)])
write.csv(x[order(x$method),c("method","reviewed","major_match","distinct_fine","all_fine_calls")],"three_panel_counts.csv",row.names=FALSE)
z <- rbind(data.frame(method=x$method,panel="Broad-family matches",count=x$major_match,label=paste0(x$major_match,"/",x$reviewed)),
           data.frame(method=x$method,panel="Distinct subtypes in\nfamily-matched clusters",count=x$distinct_fine,label=as.character(x$distinct_fine)),
           data.frame(method=x$method,panel="Clusters with a fine label",count=x$all_fine_calls,label=paste0(x$all_fine_calls,"/26")))
z$panel <- factor(z$panel,levels=unique(z$panel))
p <- ggplot(z,aes(method,count,fill=method))+
  geom_col(width=.66)+geom_text(aes(label=label),vjust=-.45,fontface="bold",size=6)+
  facet_wrap(~panel,nrow=1,scales="free_y")+scale_fill_manual(values=c("PCOS-CAM B"="#009E00","SingleR"="#FF8000","ScType"="#9966B8"))+
  scale_y_continuous(breaks=function(lim) if (lim[2]>15) seq(0,26,5) else seq(0,10,if(lim[2]>7) 2 else 1),expand=expansion(mult=c(0,.20)))+
  labs(x=NULL,y="Count",subtitle="GSE268919   |   Panels 1-2: 10 reviewed clusters   |   Panel 3: all 26 clusters")+
  theme_classic(base_size=17)+theme(legend.position="none",axis.line=element_blank(),
  panel.border=element_rect(fill=NA,colour="black",linewidth=1.1),strip.background=element_blank(),
  strip.text=element_text(face="bold",size=17,margin=margin(b=9)),
  axis.text=element_text(colour="black",face="bold"),axis.text.x=element_text(angle=35,hjust=1),
  axis.title.y=element_text(face="bold",margin=margin(r=10)),axis.ticks=element_line(linewidth=.75),
  panel.spacing=grid::unit(1.4,"lines"),plot.subtitle=element_text(size=13,margin=margin(b=12)),
  plot.margin=margin(12,16,10,12))
for(ext in c("png","pdf"))ggsave(paste0("native_subtype_bars.",ext),p,width=13.2,height=5.2,dpi=300,bg="white")
