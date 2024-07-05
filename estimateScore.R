
rm(list=ls())

library(Hmisc)

setwd("D:\\research\\tf_activities\\append\\estimate")

load("D:/research/tf_activities/append/estimate/BLCA_estimateScore.rda")
load("D:/research/tf_activities/append/estimate/TCGA_group.rda")

BLCA<-estimateScore[["TCGA-BLCA_FPKM"]]
BLCA$sample<-rownames(BLCA)
sur.cat$sample<-rownames(sur.cat)

data<-merge(sur.cat,BLCA)
dat<-data[,c(4,6,7,8,9,10)]
plotA=melt(dat,id="score",variable.name = "Immune.Checkpoint",value.name = "Expression")
plotAA<-na.omit(plotA)
plotAA$score<-capitalize(plotAA$score)
colnames(plotAA)[1]<-"Subtype"
plotAA$Expression<-as.numeric(plotAA$Expression)

score<-colnames(data)[6:10]
plot_list<-list()

col<-c("cadetblue1","palevioletred1")
for (i in 1:length(score)) {
  a<-plotAA[which(plotAA$Immune.Checkpoint == score[i]),]
  
  OutVals = boxplot(a$Expression)$out
  a<-a[which(! a$Expression %in% OutVals),]#É¾³ýÀëÉ¢µã
  
  p<-ggplot(data = a,aes(x = Subtype,y =Expression,fill = Subtype))+geom_violin()+
    geom_boxplot(width=0.2)+
    scale_x_discrete(name =NULL) +
    scale_y_continuous(name = paste0(score[i]))+
    scale_fill_manual(values=col,name="Rick score")+
    theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
          axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
          axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
          legend.title=element_text(size=13,family="serif",color="black",face="bold"),
          legend.text=element_text(size=10,family="serif",color="black",face="bold"),
          panel.background = element_rect(fill="lightgoldenrodyellow",colour="lightgoldenrodyellow",size=0.5,linetype="solid"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          legend.position="none")+
    stat_compare_means(aes(x = Subtype,y = Expression),
                       method="wilcox.test",label = "p.format",
                       label.x = 1.3)
  plot_list[[i]]<-p
  names(plot_list)[i]<-paste0(score[i])
}

figure<-plot_grid(
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  plot_list[[5]],
  ncol=2)

pdf("violinplot.pdf",width = 10,height = 8)
figure
dev.off()

pdf("violinplot_CYT.pdf",width = 5,height = 4)
plot_list[[1]]
dev.off()