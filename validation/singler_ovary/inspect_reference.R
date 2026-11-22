library(rhdf5)
f <- "GSE267729.h5ad"
items <- h5ls(f,recursive=2)
print(subset(items,group %in% c("/","/obs","/var","/raw","/X"))[,c("group","name","dim")],row.names=FALSE)
cat("\nCategories in reference annotations\n")
for (key in items$name[items$group=="/obs" & items$otype=="H5I_GROUP"]) {
 v <- h5read(f,paste0("obs/",key,"/categories"))
 if(length(v)<150) print(setNames(list(v),key))
}
