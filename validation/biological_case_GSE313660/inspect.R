suppressPackageStartupMessages(library(Seurat))
base <- '/Users/LGJ-15/Documents/PCOS/01_PCOS生信/RawData/GSE313660'
obj <- readRDS(file.path(base, 'GSE313660_Wegrzynowicz_PCOS_GEO.Rds'))
cat('class:', class(obj), '\n')
cat('dimensions:', paste(dim(obj), collapse=' x '), '\n')
cat('assays:', paste(Assays(obj), collapse=', '), '\n')
for (a in Assays(obj)) {
  cat('assay', a, 'features', nrow(obj[[a]]), 'cells', ncol(obj[[a]]), '\n')
  cat('first features:', paste(head(rownames(obj[[a]])), collapse=', '), '\n')
}
m <- obj[[]]
cat('metadata columns:', paste(colnames(m), collapse=', '), '\n')
for (k in intersect(c('pcos_status','tissue','patient','identities','broad','seurat_clusters','hashing'), colnames(m))) {
  cat('\n', k, '\n')
  print(table(m[[k]], useNA='ifany'))
}
cat('\npatient by condition and tissue\n')
print(with(m, table(patient, pcos_status, tissue)))
cat('cell names match metadata:', identical(colnames(obj), rownames(m)), '\n')
cat('RNA count layer:', paste(Layers(obj[['RNA']]), collapse=', '), '\n')
