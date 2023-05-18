###############Female immune immune cell infiltration####################
rm(list = ls()) 
inputFile="C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/normalize_GSE_female.txt"     
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/2.immune cell infiltration")     
source("geoCRG.CIBERSORT.R")       

outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)

outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab1=rbind(id=colnames(outTab),outTab)
write.table(outTab1, file="female.CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

con=grepl("_normal", rownames(outTab), ignore.case=T)
treat=grepl("_patient", rownames(outTab), ignore.case=T)
conData=outTab[con,]
treatData=outTab[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))
rt=rbind(conData,treatData)

pdf(file="Fig2B.female.barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"normal",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"patient",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

library(reshape)
library(ggpubr)
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data1=cbind(as.data.frame(t(data)), Type)
data1=melt(data1, id.vars=c("Type"))
colnames(data1)=c("Type", "Immune", "Expression")

library(vioplot)      
outTab2=data.frame()
pdf(file="Fig2D.female.vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab2=rbind(outTab2,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("normal", "patient"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

library(dplyr)
mean = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),mean)
sd = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),sd)
immuneDiff = left_join(mean,sd,by=c("Group.1","Group.2"))

immuneDiff <- immuneDiff[which(immuneDiff$Group.2 %in% outTab2$Cell),] %>%
  group_by(Group.2) %>%
  dplyr::select(Cell=Group.2, group=Group.1,mean=x.x,sd=x.y) %>%
  left_join(outTab2,by="Cell")

write.table(immuneDiff,file="Table3.female_immuneDiff.xls",sep="\t",row.names=F,quote=F)

###############Male immune immune cell infiltration####################
rm(list = ls()) 
inputFile="C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/normalize_GSE_male.txt"     
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/2.immune cell infiltration")     
source("geoCRG.CIBERSORT.R")       

outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)

outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab1=rbind(id=colnames(outTab),outTab)
write.table(outTab1, file="male.CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

con=grepl("_normal", rownames(outTab), ignore.case=T)
treat=grepl("_patient", rownames(outTab), ignore.case=T)
conData=outTab[con,]
treatData=outTab[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))
rt=rbind(conData,treatData)

pdf(file="Fig2A.male.barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"normal",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"patient",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

library(reshape)
library(ggpubr)
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data1=cbind(as.data.frame(t(data)), Type)
data1=melt(data1, id.vars=c("Type"))
colnames(data1)=c("Type", "Immune", "Expression")

library(vioplot)      
outTab2=data.frame()
pdf(file="Fig2C.male.vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab2=rbind(outTab2,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("normal", "patient"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

library(dplyr)
mean = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),mean)
sd = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),sd)
immuneDiff = left_join(mean,sd,by=c("Group.1","Group.2"))

immuneDiff <- immuneDiff[which(immuneDiff$Group.2 %in% outTab2$Cell),] %>%
  group_by(Group.2) %>%
  dplyr::select(Cell=Group.2, group=Group.1,mean=x.x,sd=x.y) %>%
  left_join(outTab2,by="Cell")


write.table(immuneDiff,file="Table2.male_immuneDiff.xls",sep="\t",row.names=F,quote=F)