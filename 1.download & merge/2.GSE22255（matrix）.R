rm(list = ls()) 
geo.data.series.number = "GSE22255"
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge")


Sys.setenv(LIB_XML = "$(MINGW_PREFIX)") 

library(Biobase)
library(GEOquery)
library(limma) 
library(gplots)
library(stringr)
library(impute)
library(dplyr)
library(limma)
library(tibble)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)
gset <- getGEO(geo.data.series.number, GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
samples <- colnames(exprs(gset))
geo.platform.number = gset@annotation   
probes_expr = exprs(gset);dim(probes_expr)
phenoDat = pData(gset) 
write.csv(phenoDat, file=paste(geo.data.series.number,"_Type.csv",sep=""), row.names=F)
ex = exprs(gset)
imputed_gene_exp = impute.knn(ex,k=10,rowmax = 0.5,
                              colmax=0.8,maxp =3000, rng.seed=362436069)
ex = imputed_gene_exp$data
ex = as.data.frame(ex)
ex$ID = rownames(ex)
ids = as.data.frame(gset@featureData@data[["ID"]])
ids$Gene.symbol = gset@featureData@data[["Gene.symbol"]]  
colnames(ids) = c("ID","Gene.Symbol")
ids = ids[-grep('///',ids$Gene.Symbol),]         
exprSet = inner_join(ids,ex,by = 'ID')         
exprSet = na.omit(exprSet)
exprSet = exprSet[-which(exprSet$Gene.Symbol==''),]
exprSet= avereps(exprSet[,-c(1,2)],              
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)
exprSet = normalizeBetweenArrays(exprSet)
exprSet = as.data.frame(exprSet)
exprSet <- rownames_to_column(exprSet, var="geneNames")
exprSet <- exprSet[order(exprSet[,1]),]
write.table(exprSet, file=paste(geo.data.series.number,".csv",sep=""), row.names=F, col.names=T, quote=F, sep=",")



