# Native SingleR, with the general mouse reference supplied by celldex.
library(SingleR);library(celldex)
out <- "validation/native_comparison"
ref <- MouseRNAseqData()
z <- readRDS(commandArgs(TRUE)[1])
p <- SingleR(test=z$m,ref=ref,labels=ref$label.fine,clusters=z$cl)
mapping <- unique(data.frame(fine=ref$label.fine,main=ref$label.main))
stopifnot(!anyDuplicated(mapping$fine))
tab <- data.frame(cluster=rownames(p),SingleR=ifelse(is.na(p$pruned.labels),"Unknown",p$pruned.labels),unpruned=p$labels)
tab$reference_main <- mapping$main[match(tab$SingleR,mapping$fine)]
write.csv(tab,file.path(out,"SingleR_native.csv"),row.names=FALSE)
