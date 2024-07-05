
rm(list=ls())

library(survminer)
library(survival)
library(survcomp)


setwd("D:\\research\\tf_activities\\cox")

load("D:/research/tf_activities/tf_activities_7.rda")
load("D:/research/tf_activities/cox/risk_TCGA.rda")
tf<-read.table("D:\\research\\tf_activities\\cox\\ÏÔÖøÄ£¿é»ùÒò.txt")

rownames(risk_TCGA)<-gsub("\\.","-",rownames(risk_TCGA))
risk_TCGA<-risk_TCGA[,c(2,3)]


exp<-t(result[["TCGA-BLCA_FPKM"]])
exp<-exp[,tf[,1]]


cox_data<-merge(risk_TCGA,exp,by = "row.names")

cox_result<-matrix(nrow=ncol(cox_data)-3,ncol=6)
for (i in 4:ncol(cox_data)) {
  Bcox<-coxph(Surv(times, status)~cox_data[,i],data=cox_data)
  summcph<-summary(Bcox)
  cox_result[i-3,1]<-summcph$conf.int[1]
  cox_result[i-3,2]<-summcph$conf.int[3]
  cox_result[i-3,3]<-summcph$conf.int[4]
  cox_result[i-3,4]<-as.matrix(summcph$logtest)[3]
  cox_result[i-3,5]<-as.matrix(summcph$sctest)[3]
  cox_result[i-3,6]<-summcph$coefficients[5]
  print(i)
}
rownames(cox_result)=colnames(cox_data)[4:ncol(cox_data)]
colnames(cox_result)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
cox_sig<-cox_result[cox_result[,6] < 0.01,]

save(cox_result,file = "cox_result.rda")
save(cox_sig,file = "cox_sig.rda")

