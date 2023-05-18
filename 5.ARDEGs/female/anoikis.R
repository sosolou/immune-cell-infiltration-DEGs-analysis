rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/5.ARDEGs/female")


library(dplyr)
library(tidyverse)
anoikisgene=read.table("anoikisgene.txt", header=F, sep="\t", check.names=F) %>%
  distinct()

norm=read.table("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/normalize_GSE_female.txt", header=T, sep="\t", check.names=F)
anoikis.norm=norm[which(norm$Gene %in% anoikisgene$V1),]  %>% as_tibble() %>%
  column_to_rownames("Gene")

library(limma)
library(pheatmap)
exp <- as.matrix(anoikis.norm) 
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
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]

diffGeneExp=rt[rownames(diffSig),]



Sig=ifelse((allDiff$adj.P.Val<adjustP) & (abs(allDiff$logFC)>logFoldChange), ifelse(allDiff$logFC>logFoldChange,"Up","Down"), "Stable")
rt=cbind(allDiff, Sig=Sig) %>%
  rownames_to_column("id")

library(ggplot2)
library(ggrepel)
p =  ggplot(
  rt, aes(x = logFC, y = -log10(adj.P.Val), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue","black","red"))+ 
  geom_vline(xintercept=c(-logFoldChange,logFoldChange),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(adjustP),lty=3,col="black",lwd=0.5) +
  labs(x="log2(Fold Change)",
       y="-log10 (P-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) +
  geom_label_repel(data=filter(rt, ((rt$adj.P.Val<adjustP) & (abs(rt$logFC)>logFoldChange))),
                      box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                      size=1.8, aes(label=id),show.legend = FALSE) 
  

pdf(file="Fig5F.vol.pdf", width=7, height=6.1)
print(p)
dev.off()


anoikis.exp = diffGeneExp
anoikis.t <- t(anoikis.exp) 
library(limma)
library(pheatmap)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(anoikis.exp))
names(Type)=colnames(anoikis.exp)
Type=as.data.frame(Type)
pdf(file="Fig5G.heatmap.pdf", width=8, height=6)
pheatmap(anoikis.exp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()


library("RCircos") 
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)


pdf(file="Fig5H.RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
RCircos.Gene.Label.Data=RCircos.Gene.Label.Data[which(RCircos.Gene.Label.Data$Gene %in% colnames(anoikis.t)),]
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


GO=read.table("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/female.GO.txt", header=T, sep="\t", check.names=F)
GO.anoikis <-GO %>% as_tibble() %>% 
  separate_rows(geneID, sep = "/") 
GO.anoikis1 = GO.anoikis[which(GO.anoikis$geneID %in% rownames(anoikis.exp)),]
GO.anoikis1 = GO[which(GO$ID %in% GO.anoikis1$ID),] %>%
  arrange(ONTOLOGY,desc(Count))


GO.anoikis <- GO.anoikis1 %>%
  group_by(ONTOLOGY) %>% do(head(., n = 8)) %>%
  arrange(ONTOLOGY,Count)

library(ggplot2)
GO_term_order=factor(as.integer(rownames(GO.anoikis)),labels=GO.anoikis$Description)
pdf(file="Fig5I.anoikis-related TOP8 GO term.pdf", width=7, height=5)
ggplot(data=GO.anoikis, aes(x=GO_term_order,y=Count, fill=ONTOLOGY)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GO term") + ylab("Num of Genes") + theme_bw()
dev.off()


library(limma)
library(reshape2)
library(ggplot2)

immFile="C:/Users/Administrator/Desktop/ischemic stroke & cell death/2.immune cell infiltration/female.CIBERSORT-Results.txt"     


group=gsub("(.*)\\_(.*)", "\\2", colnames(anoikis.exp))
data=anoikis.exp[,group=="patient",drop=F]
data=t(data)


immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]


outTab=data.frame()
for(cell in colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene in colnames(data)){
    x=as.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}


outTab$cor=as.numeric(outTab$cor)
pdf(file="Fig5J.immun-cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   
        axis.text.y = element_text(size = 8, face = "bold")) +       
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
  scale_x_discrete(position = "bottom")      
dev.off()

