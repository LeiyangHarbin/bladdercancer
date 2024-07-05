
rm(list=ls())

library(factoextra)

setwd("D:\\research\\tf_activities\\cox")

load("D:/research/tf_activities/tf_activities_7.rda")
load("cox_sig.rda")

exp<-t(result[["TCGA-BLCA_FPKM"]])
data<-exp[,rownames(cox_sig)]

p=fviz_nbclust(data, kmeans, method = "silhouette",linecolor = "darkred")
p1<-p+geom_vline(xintercept = 2, linetype = 2,col="blue")
p1=p1+theme(axis.title.x = element_text(size = 15,family = "serif",face = "bold"),
            axis.title.y = element_text(size = 15,family = "serif",face = "bold"),
            axis.text.x = element_text(size = 12,family = "serif",face = "bold"),
            axis.text.y = element_text(size =12,family = "serif",face = "bold")
            ,title=element_text(size = 18,family = "serif",face = "bold"))
p1=p1+theme(plot.title = element_text(hjust = 0.5))
p1
ggsave(p1,filename = "K-mean.pdf")

km<- kmeans(data, 2, nstart = 2)
label=km$cluster
label=data.frame(label)
save(label,file = "k-mean.rda")


load("D:/research/tf_activities/cox/risk_TCGA.rda")

rownames(risk_TCGA)<-gsub("\\.","-",rownames(risk_TCGA))
surdata<-merge(risk_TCGA,label,by = "row.names")

fit1 <- survfit(Surv(times, status)~label, data = surdata)  


ggsurvplot(fit1,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


