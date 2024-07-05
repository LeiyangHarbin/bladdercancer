
#################################################################
rm(list=ls())

library(viper)
library(bcellViper)
library(mixtools)
library(corto)
library(HGNChelper)
library(diggit)

setwd("D:\\research\\tf_activities\\master")

load("D:/research/data_new/Bladder_7.rda")
load("TCGA_group.rda")
BLCA<-Bladder_7[["TCGA-BLCA_FPKM"]]
exp<-t(exprs(BLCA))
high<-sur.cat[sur.cat$score == "high",]
low<-sur.cat[sur.cat$score == "low",]
exp<-data.frame(exp)
high_exp<-exp[rownames(high),]
low_exp<-exp[rownames(low),]
high_exp<-as.matrix(t(high_exp))
low_exp<-as.matrix(t(low_exp))

signature<-rowTtest(low_exp,high_exp)
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) *sign(signature$statistic))[, 1]
nullmodel <- ttestNull(low_exp, high_exp, per = 1000,repos = TRUE, verbose = FALSE)

save(signature,file = "signature.rda")
save(nullmodel,file = "nullmodel.rda")


############################
rm(list=ls())

library(viper)
library(dplyr)

setwd("D:\\research\\tf_activities\\master")

load("signature.rda")
load("nullmodel.rda")


data(dorothea_hs,package="dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A","B","C","D","E"))

TF<-rownames(as.matrix(table(regulons$tf)))

regulon_1<-list()
for (i in 1:length(TF)) {
  A1<-which(regulons$tf==TF[i])
  A2<-unlist(regulons[A1,4])
  names(A2)<-unlist(regulons[A1,3])
  regulon_1[i]<-list(A2)
  names(regulon_1)[i]<-TF[i]
}

mrs <- msviper(signature, regulon_1, nullmodel, verbose = FALSE)
summary(mrs)


mrs1<-mrs
save(mrs1,file = "mrs1.rda")



