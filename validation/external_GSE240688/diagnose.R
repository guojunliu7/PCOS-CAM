# Run from PCOS-CAM root; unchanged settings, diagnostics only.
source('Global.R');source('hierarchy.R')
out <- 'validation/external_GSE240688'
d <- read_pcos_dictionary('data/File_1.xlsx')
h <- read_pcos_hierarchy('data/File_1.xlsx',d,'data/label_hierarchy.csv')
x <- read.csv(file.path(out,'cluster_mean_expression.csv'))
z <- pcos_cam(x,d,h,'Human')
e <- z$parent_evidence
rows <- lapply(seq_len(nrow(z$annotations)),function(j) {
 a <- z$annotations[j,]
 id <- as.character(a$cluster)
 fam <- a$parent_candidates
 hits <- e[e$cluster==id & e$cell_cluster==fam,,drop=FALSE]
 hits <- hits[order(-hits$contribution),,drop=FALSE]
 if(!nrow(hits)) return(NULL)
 data.frame(cluster=id,predicted_family=fam,parent_score=a$parent_score,
   top_gene=hits$marker_gene[1],top_contribution=hits$contribution[1],
   top_fraction=hits$contribution[1]/a$parent_score,
   second_gene=if(nrow(hits)>1) hits$marker_gene[2] else NA_character_,
   second_contribution=if(nrow(hits)>1) hits$contribution[2] else NA_real_,
   parent_stable=a$parent_stable)
})
write.csv(do.call(rbind,rows),file.path(out,'PCOS_CAM_parent_contribution_audit.csv'),row.names=FALSE)
