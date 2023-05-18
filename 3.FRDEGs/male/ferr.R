rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/3.FRDEGs/male")

library(dplyr)
library(tidyverse)
ferrgene=read.table("ferrgene.txt", header=F, sep="\t", check.names=F) %>%
  distinct()

norm=read.table("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/normalize_GSE_male.txt", header=T, sep="\t", check.names=F)
ferr.norm=norm[which(norm$Gene %in% ferrgene$V1),]  %>% as_tibble() %>%
  column_to_rownames("Gene")

library(limma)
library(pheatmap)
exp <- as.matrix(ferr.norm) 
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
allDiffOut <- rbind(id=colnames(allDiff),allDiff)
rm(design,cont.matrix,fit,fit2)


logFoldChange=0.58
adjustP=0.05

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
diffSig <- diffSig[order(-diffSig$logFC),]
diffSigOut <- rbind(id=colnames(diffSig),diffSig)

diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]

diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]

rt1 <- diffSig %>% rownames_to_column("id")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

pvalueFilter=0.05       
qvalueFilter=0.05       


colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="p.adjust"
}


colnames(rt1)[1]="Gene"
genes=as.vector(rt1[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        


kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$p.adjust<pvalueFilter & GO$qvalue<qvalueFilter),]



showNum=8
if(nrow(GO)<30){
  showNum=nrow(GO)
}



ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]


BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5


pdf("Fig3E.GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                    
                    })
circos.clear()

middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)

main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)

logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

diffGeneExp=rt[rownames(diffSig),]
diffGeneExpOut=rbind(id=colnames(diffGeneExp),diffGeneExp)


Sig=ifelse((allDiff$adj.P.Val<adjustP) & (abs(allDiff$logFC)>logFoldChange), ifelse(allDiff$logFC>logFoldChange,"Up","Down"), "Stable")
rt=cbind(allDiff, Sig=Sig) %>%
  rownames_to_column("id")


library(ggplot2)
library(ggrepel)
p =  ggplot(
  rt, aes(x = logFC, y = -log10(adj.P.Val), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("black","red"))+ 
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
  

pdf(file="Fig3A.vol.pdf", width=7, height=6.1)
print(p)
dev.off()


ferr.exp = diffGeneExp

ferr.t <- t(ferr.exp) 

library(reshape2)
exp=as.data.frame(ferr.t)
Type=gsub("(.*)\\_(.*)", "\\2", rownames(exp))
Type=gsub("normal","healthy control group",Type)
Type=gsub("patient","patient group",Type)
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")


library(ggpubr)
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="*** p<0.001, ** p<0.01, * p<0.05",
            ylab="Gene expression",
            legend.title="Type",
            palette = "jama",
            add="point",
            width=0.8)+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")


pdf(file="Fig3B.boxplot.pdf", width=7, height=5)
print(p1)
dev.off()

cutoff=0.1                     
rt=ferr.exp[,-1:-conNum]
data=t(rt)
cordata=cor(data)

mydata = cordata
upper = upper.tri(mydata)
mydata[upper] = NA

df = data.frame(gene=rownames(mydata),mydata)
dfmeltdata = melt(df,id="gene")
dfmeltdata = dfmeltdata[!is.na(dfmeltdata$value),]
dfmeltdata = dfmeltdata[dfmeltdata$gene!=dfmeltdata$variable,]
dfmeltdata = dfmeltdata[abs(dfmeltdata$value)>cutoff,]

library(igraph)
corweight = dfmeltdata$value
weight = corweight+abs(min(corweight))+5
d = data.frame(p1=dfmeltdata$gene,p2=dfmeltdata$variable,weight=dfmeltdata$value)
g = graph.data.frame(dfmeltdata,directed = FALSE)


E(g)$weight = weight
E(g)$color = ifelse(corweight>0,rgb(254/255,67/255,101/255,abs(corweight)),rgb(0/255,0/255,255/255,abs(corweight)))
V(g)$size = 8
V(g)$shape = "circle"
V(g)$lable.cex = 1.2
V(g)$color = "white"


pdf("Fig3D.co-expression network among differential ferr-genes.pdf", width=7, height=6)
layout(matrix(c(1,1,1,0,2,0),byrow=T,nc=3),height=c(6,1),width=c(3,4,3))
par(mar=c(1.5,2,2,2))
vertex.frame.color = NA
plot(g,layout=layout_nicely,vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,vertex.label.color="black",vertex.frame.color=vertex.frame.color,edge.color=E(g)$color,vertex.label.cex=V(g)$lable.cex,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.4)


color_legend = c(rgb(254/255,67/255,101/255,seq(1,0,by=-0.01)),rgb(0/255,0/255,255/255,seq(0,1,by=0.01)))
par(mar=c(2,2,1,2),xpd = T,cex.axis=1.6,las=1)
barplot(rep(1,length(color_legend)),border = NA, space = 0,ylab="",xlab="",xlim=c(1,length(color_legend)),horiz=FALSE,
        axes = F, col=color_legend,main="")
axis(3,at=seq(1,length(color_legend),length=5),c(1,0.5,0,-0.5,-1),tick=FALSE)
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

pdf(file="Fig3C.RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()


RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
RCircos.Gene.Label.Data=RCircos.Gene.Label.Data[which(RCircos.Gene.Label.Data$Gene %in% colnames(ferr.t)),]
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


library(limma)
library(reshape2)
library(tidyverse)
library(dplyr)
library(ggplot2)

immFile="C:/Users/Administrator/Desktop/ischemic stroke & cell death/2.immune cell infiltration/male.CIBERSORT-Results.txt"     


group=gsub("(.*)\\_(.*)", "\\2", colnames(ferr.exp))
data=ferr.exp[,group=="patient",drop=F]
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
pdf(file="Fig3F.immun-cor.pdf", width=7, height=5)
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

