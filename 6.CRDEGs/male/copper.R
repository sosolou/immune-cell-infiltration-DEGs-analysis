rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/6.CRDEGs/male")


library(dplyr)
library(tidyverse)
coppergene=read.table("coppergene.txt", header=F, sep="\t", check.names=F) %>%
  distinct()

norm=read.table("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/normalize_GSE_male.txt", header=T, sep="\t", check.names=F)
copper.norm=norm[which(norm$Gene %in% coppergene$V1),]  %>% as_tibble() %>%
  column_to_rownames("Gene")


library(limma)
library(pheatmap)
exp <- as.matrix(copper.norm) 
dimnames <- list(rownames(exp),colnames(exp))
rt <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt <- avereps(rt)
group=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
conNum <- table(group)[1]
treatNum <- table(group)[2]

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
diffSigOut <- rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="male.diff.txt", sep="\t", quote=F, col.names=F)



