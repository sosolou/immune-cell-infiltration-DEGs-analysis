rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/7.Machine learning")


library(dplyr)
library(tidyverse)
norm=read.table("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge/normalize_GSE_female.txt", header=T, sep="\t", check.names=F)
ml=c("IL6","PDK4","LCN2","FAR1","SLC40A1","TLR4","OLFM4",
     "CEACAM6","MMP9","TNF","C5AR1","CAMP","NLRC4","ELANE","CTSG",
     "CD163","CD96","GZMK","IL32","CD27","CD3E","CD3D","CD2","MT1X"
)
ml.norm=norm[which(norm$Gene %in% ml),]  %>% as_tibble() %>%
  column_to_rownames("Gene")


library(limma)
library(pheatmap)
exp <- as.matrix(ml.norm) 
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

diffGeneExp=rt[rownames(diffSig),]
ml.exp = diffGeneExp

library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

data=t(ml.exp)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group


inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]


control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)


mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)


mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)


mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)



p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="normal", 0, 1)


explainer_rf=explain(mod_rf, label = "RF",
                     data = test, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_rf=model_performance(explainer_rf)

explainer_svm=explain(mod_svm, label = "SVM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_svm=model_performance(explainer_svm)

explainer_xgb=explain(mod_xgb, label = "XGB",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)

explainer_glm=explain(mod_glm, label = "GLM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_glm=model_performance(explainer_glm)


pdf(file="Fig6C.boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()

pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
ciVec1=as.numeric(ci.auc(roc1, method="delong"))
ciVec2=as.numeric(ci.auc(roc2, method="delong"))
ciVec3=as.numeric(ci.auc(roc3, method="delong"))
ciVec4=as.numeric(ci.auc(roc4, method="delong"))
pdf(file="Fig6D.ROC.pdf", width=8, height=8)
plot(smooth(roc1,method="density"), print.auc=F, legacy.axes=T, main="", col="#8ECFC9")
plot(smooth(roc2,method="density"), print.auc=F, legacy.axes=T, main="", col="#FFBE7A", add=T)
plot(smooth(roc3,method="density"), print.auc=F, legacy.axes=T, main="", col="#BEB8DC", add=T)
plot(smooth(roc4,method="density"), print.auc=F, legacy.axes=T, main="", col="#FF8884", add=T)
legend('bottom',
       c(paste0('   RF: ',sprintf("%.03f",roc1$auc),' (',paste0("95%CI: ",sprintf("%.03f",ciVec1[1]),"-",sprintf("%.03f",ciVec1[3])),')'),
         paste0('SVM: ',sprintf("%.03f",roc2$auc),' (',paste0("95%CI: ",sprintf("%.03f",ciVec2[1]),"-",sprintf("%.03f",ciVec2[3])),')'),
         paste0('XGB: ',sprintf("%.03f",roc3$auc),' (',paste0("95%CI: ",sprintf("%.03f",ciVec3[1]),"-",sprintf("%.03f",ciVec3[3])),')'),
         paste0('GLM: ',sprintf("%.03f",roc4$auc),' (',paste0("95%CI: ",sprintf("%.03f",ciVec4[1]),"-",sprintf("%.03f",ciVec4[3])),')')),
       col=c("#8ECFC9","#FFBE7A","#BEB8DC","#FF8884"), lwd=2, bty = 'n')
dev.off()


importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)

pdf(file="Fig6F.importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-3):(ncol(data)+1)),],
     importance_svm[c(1,(ncol(data)-3):(ncol(data)+1)),],
     importance_xgb[c(1,(ncol(data)-3):(ncol(data)+1)),],
     importance_glm[c(1,(ncol(data)-3):(ncol(data)+1)),])
dev.off()

geneNum=5     
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="female.importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="female.importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="female.importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="female.importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)

