# Reactome overrepresentation among patient-level FDR-significant genes.
# This answers a different question from CAMERA's competitive gene-set test.
# Run from the PCOS-CAM project root.
suppressPackageStartupMessages({library(msigdbr);library(clusterProfiler)})
out <- 'validation/biological_case_Theca_7942968'
sets <- msigdbr(species='Homo sapiens',collection='C2',subcollection='CP:REACTOME')
term <- unique(sets[,c('gs_name','gene_symbol')])
for(prefix in c('theca','fine5')) {
  path <- if(prefix=='theca') 'theca_PCOS_vs_control_DE.csv' else 'fine5_PCOS_vs_control_DE.csv'
  de <- read.csv(file.path(out,path))
  sig <- de$gene[de$FDR<.05]
  result <- enricher(gene=sig,universe=de$gene,
                     TERM2GENE=term,pAdjustMethod='BH',
                     pvalueCutoff=1,qvalueCutoff=1,
                     minGSSize=10,maxGSSize=500)
  tab <- as.data.frame(result)
  write.csv(tab,file.path(out,paste0(prefix,'_Reactome_DEG_overrepresentation.csv')),
            row.names=FALSE)
  cat(prefix,'DEGs:',length(sig),'ORA pathways adjusted p<.05:',
      sum(tab$p.adjust<.05),'\n')
  if(nrow(tab)) print(head(tab[,c('ID','Count','p.adjust')],8),row.names=FALSE)
  for(direction in c('up','down')) {
    selected <- de$gene[de$FDR<.05 & if(direction=='up') de$logFC>0 else de$logFC<0]
    directional <- enricher(gene=selected,universe=de$gene,
                            TERM2GENE=term,pAdjustMethod='BH',
                            pvalueCutoff=1,qvalueCutoff=1,
                            minGSSize=10,maxGSSize=500)
    direction_tab <- as.data.frame(directional)
    write.csv(direction_tab,
              file.path(out,paste0(prefix,'_Reactome_DEG_ORA_',direction,'.csv')),
              row.names=FALSE)
    cat(prefix,direction,'DEGs:',length(selected),'ORA FDR<.05:',
        sum(direction_tab$p.adjust<.05),'\n')
  }
}
