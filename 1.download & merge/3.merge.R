rm(list = ls()) 

setwd("C:/Users/Administrator/Desktop/ischemic stroke & cell death/1.download & merge")
library(tidyverse)
library(dplyr)
library(sva)
library(limma)

GSE1 <- read.table(file = "GSE22255.csv", header = TRUE, sep="," )
GSE2 <- read.table(file = "GSE16561.csv", header = TRUE, sep="," )
GSE <- inner_join(GSE1, GSE2, by = "geneNames")

GSE1_type <- read.table(file = "GSE22255_type.csv", header = TRUE, sep="," ) %>%
  dplyr::select(sample = geo_accession, gender = characteristics_ch1.1, group = "affected.status..disease.state..ch1") %>%
  mutate(data="GSE22255")

GSE2_type <- read.table(file = "GSE16561_type.csv", header = TRUE, sep="," ) %>%
  dplyr::select(sample = geo_accession, gender = characteristics_ch1, group = description) %>%
  mutate(data="GSE16561")



GSE_type <- rbind(GSE1_type,GSE2_type) %>%
  distinct()  %>%
  mutate(
    gender = case_when(
      gender == "gender: female" | gender == "gender: Female" ~ "female",
      gender == "gender: male" | gender == "gender: Male" ~ "male")
  ) %>%
  mutate(
    Group = case_when(
      group == "IS patient" | group == "Stroke" ~ "patient",
      group == "Control" | group == "control" ~ "normal")
  ) 


table(GSE_type$data,GSE_type$Group,GSE_type$gender)



rt <- as.matrix(GSE)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType <- c(rep(1,40),rep(2,63))
modType <- c(rep("Control",20),rep("IS",20),rep("Control",24),rep("IS",39))
mod <- model.matrix(~as.factor(modType))
normalize_GSE <- ComBat(data, batchType, mod, par.prior=TRUE)#对数据进行批次校正
normalize_GSE <- as.data.frame(normalize_GSE)

rm(mod,rt,batchType,data,dimnames,exp)

normalize_GSE_type <- GSE_type %>%
  arrange(Group, data) %>%
  mutate(id=paste0(sample,"_",Group))

normalize_GSE <- normalize_GSE[,c(normalize_GSE_type$sample)]
colnames(normalize_GSE) <- paste0(colnames(normalize_GSE),"_",normalize_GSE_type$Group) 


normalize_GSE_type_male <- normalize_GSE_type %>%
  filter(gender == "male")

normalize_GSE_type_female <- normalize_GSE_type %>%
  filter(gender== "female")

normalize_GSE_male <- normalize_GSE[,c((colnames(normalize_GSE) %in% normalize_GSE_type_male$id))]

normalize_GSE_female <- normalize_GSE[,(colnames(normalize_GSE) %in% normalize_GSE_type_female$id)]


normalize_GSE <- normalize_GSE %>% rownames_to_column("Gene") 
write.table(normalize_GSE, file="normalize_GSE.txt", sep="\t", quote=F, col.names=T,row.names =F)
write.table(normalize_GSE_type, file="normalize_GSE_type.txt", sep="\t", quote=F, col.names=T,row.names =F)



normalize_GSE_male <- normalize_GSE_male %>% rownames_to_column("Gene") 
write.table(normalize_GSE_male, file="normalize_GSE_male.txt", sep="\t", quote=F, col.names=T,row.names =F)
write.table(normalize_GSE_type_male, file="normalize_GSE_type_male.txt", sep="\t", quote=F, col.names=T,row.names =F)


normalize_GSE_female <- normalize_GSE_female %>% rownames_to_column("Gene") 
write.table(normalize_GSE_female, file="normalize_GSE_female.txt", sep="\t", quote=F, col.names=T,row.names =F)
write.table(normalize_GSE_type_female, file="normalize_GSE_type_female.txt", sep="\t", quote=F, col.names=T,row.names =F)
