
###################
library(Hmisc)

setwd("D:\\research\\tf_activities\\append\\BLCA_plot")

rm(list = ls())
load("D:\\research\\Chromosome instability\\BLCA_PLOT\\BLCA.rda")
load("D:/research/tf_activities/append/BLCA_plot/TCGA_group.rda")

sur.cat$score<-capitalize(sur.cat$score)

pdata<-pData(BLCA)
rownames(sur.cat)<-gsub("\\.","-",rownames(sur.cat))
colnames(pdata)[1]<-"sample"
sur.cat$sample<-rownames(sur.cat)
TCGA<-merge(sur.cat,pdata,by="sample")


BLCA_PLOT<-TCGA[,c(1,4,5,221,222,223,224,225,226,227,268,260,261,270,271,275)]
Mutation<-TCGA[,c(1,2,3,4,5,215,216)]
Mutation$Mutation.load<-Mutation[,6]+Mutation[,7]
BLCA_PLOT$Mutation.load<-Mutation$Mutation.load

colnames(BLCA_PLOT)[16]<-"Subclonal"

save(BLCA_PLOT,file = "BLCA_PLOT.rda")
