
rm(list = ls())

library(Hmisc)
library(reshape2)
library(ggplot2)
library(ggpubr)

setwd("D:\\research\\tf_activities\\append\\BLCA_plot")

load("BLCA_plot.rda")
data<-BLCA_PLOT[,-c(1,3)]
colnames(data)<-gsub("\\."," ",colnames(data))
colnames(data)[1]<-"Subtype"
data$Subtype<-capitalize(data$Subtype)
plotA=melt(data,id="Subtype",variable.name = "Immune.Checkpoint",value.name = "Expression")
plotAA<-na.omit(plotA)

score<-colnames(data)[2:15]
plot_list<-list()

for (i in 1:length(score)) {
  a<-plotAA[which(plotAA$Immune.Checkpoint == score[i]),]
  
  OutVals = boxplot(a$Expression)$out
  a<-a[which(! a$Expression %in% OutVals),]#É¾³ýÀëÉ¢µã
  
  p<-ggplot(data = a,aes(x = Subtype,y =Expression,fill = Subtype))+geom_violin()+
    geom_boxplot(width=0.2)+
    scale_x_discrete(name =NULL) +
    scale_y_continuous(name = paste0(score[i]))+
    labs(fill="Rick score")+
    scale_fill_manual(values = c("#FC4E07","#00AFBB"))+
    theme_bw()+
    theme(legend.text = element_text(size=14),
          axis.title = element_text(size=15),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=10),
          text = element_text(size = 12,face =  "bold",family="serif"),
          panel.background = element_rect(fill="lightgoldenrodyellow",colour="lightgoldenrodyellow",size=0.5,linetype="solid"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          legend.position = "none"
    )+ 
    stat_compare_means(aes(x = Subtype,y = Expression),
                       method="wilcox.test",label = "p.format",
                       label.x = 1.4)
  plot_list[[i]]<-p
  names(plot_list)[i]<-paste0(score[i])
  p
  ggsave(p,filename = paste0(score[i],".pdf"),width = 5,height = 4)
}