rm(list = ls()) 
geo.data.series.number = "GSE16561"
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge")

########################
Sys.setenv(LIB_XML = "$(MINGW_PREFIX)") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Biobase", force = TRUE)
#BiocManager::install("pkgname")
#BiocManager::install("GEOquery")
#BiocManager::install("limma")
#install.packages("gplots")
#install.packages("stringr")
#install.packages("Rcpp", dependencies = TRUE) 
#install.packages(c("devtools","usethis"))
#install.packages("backports")
#install.packages("dplyr")


#devtools::install_github("vqv/ggbiplot", force = TRUE)
######################################
# library package
library(Biobase)
library(GEOquery)
library(limma) 
library(stringr)


########################################
## For easier use, define arguments first
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)
# Differential expression analysis with limma
# load series and platform data from GEO
gset <- getGEO(geo.data.series.number, GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))
# group names for all samples
samples <- colnames(exprs(gset))
geo.platform.number = gset@annotation    #gset@annotation
probes_expr = exprs(gset);dim(probes_expr) 
# probes_expr = log2(probe_expr+1)  
# boxplot(probes_expr2, las = 2)


####################phenoType.csv########################
phenoDat = pData(gset) 
phenoDat = phenoDat[c(40:63,1:39),]

write.csv(phenoDat, file=paste(geo.data.series.number,"_Type.csv",sep=""), row.names=F)

############################## 


############################## log2 transform
# log2 transform
ex = exprs(gset)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("impute")
library(impute)
#KNN,to impute missing expression data, using nearest neighbor averaging
imputed_gene_exp = impute.knn(ex,k=10,rowmax = 0.5,
                              colmax=0.8,maxp =3000, rng.seed=362436069)
ex = imputed_gene_exp$data

ex = as.data.frame(ex)
ex$ID = rownames(ex)

ids = as.data.frame(gset@featureData@data[["ID"]])
ids$Gene.symbol = gset@featureData@data[["Gene.symbol"]] 
colnames(ids) = c("ID","Gene.Symbol")

library(dplyr)
exprSet = inner_join(ids,ex,by = 'ID')         
exprSet = na.omit(exprSet)
exprSet = exprSet[-which(exprSet$Gene.Symbol==''),]

library(limma)
exprSet= avereps(exprSet[,-c(1,2)],              
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)

exprSet = normalizeBetweenArrays(exprSet)

exprSet = as.data.frame(exprSet)

exprSet = exprSet[,c(40:63,1:39)] 
library(tibble)

exprSet <- rownames_to_column(exprSet, var="geneNames")
exprSet <- exprSet[order(exprSet[,1]),]
write.table(exprSet, file=paste(geo.data.series.number,".csv",sep=""), row.names=F, col.names=T, quote=F, sep=",")



