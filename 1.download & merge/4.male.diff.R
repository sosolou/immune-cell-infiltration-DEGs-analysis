rm(list = ls()) 

setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge")
library(tidyverse)
library(dplyr)
library(limma)
library("clusterProfiler")
library("org.Hs.eg.db")


normalize_GSE_male <- read.table("normalize_GSE_male.txt", header=T, sep="\t", check.names=F) %>%
  column_to_rownames("Gene")
  
normalize_GSE_type_male <- read.table("normalize_GSE_type_male.txt", header=T, sep="\t", check.names=F)
table(normalize_GSE_type_male$data,normalize_GSE_type_male$Group)

conNum <- table(normalize_GSE_type_male$Group)[1]
treatNum <- table(normalize_GSE_type_male$Group)[2]


normalize_GSE <- normalize_GSE_male

exp <- as.matrix(normalize_GSE) 
dimnames <- list(rownames(exp),colnames(exp))
rt <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt <- avereps(rt)

modType <- c(rep("normal",conNum),rep("patient",treatNum)) 
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("normal","patient")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(patient-normal,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff <- topTable(fit2,adjust='fdr',number=200000)
rm(design,cont.matrix,fit,fit2)

logFoldChange=0.58
adjustP=0.05

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
diffSig <- diffSig[order(-diffSig$logFC),]



pvalueFilter=0.05      
qvalueFilter=0.05      

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


rt <- diffSig %>% rownames_to_column("Gene")
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        

kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="male.GO.txt",sep="\t",quote=F,row.names = F)
