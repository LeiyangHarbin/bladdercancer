
############################
rm(list=ls())

library(survminer)
library(survival)
library(survcomp)

setwd("D:\\research\\tf_activities\\master")

load("mrs1.rda")
load("D:/research/tf_activities/tf_activities_7.rda")
load("D:/research/tf_activities/master/TCGA_group.rda")

BLCA<-result[["TCGA-BLCA_FPKM"]]
BLCA<-BLCA[,rownames(sur.cat)]

sig<-names(which(mrs1[["es"]][["p.value"]]<0.001))

#############################

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)

x<-sig

eg <- bitr(x, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")

go <- enrichGO(gene = eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='BP',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID',
               readable = T)
head(go)


g1<-barplot(go,drop=TRUE,showCategory=10)
g2<-g1+theme(axis.title=element_text(size=15,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=15,family="serif",color="black"),
             axis.text.y=element_text(size=15,family="serif",color="black"),
             legend.text=element_text(size=15,family="serif",color="black"),
             legend.title=element_text(size=15,family="serif",color="black",face="bold")
)
g3<-g2+ggtitle("TCGA")+
  xlab("Count")+
  ylab("GO")+
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
g3
ggsave(g3,filename = "master_regulon_GO.pdf",height = 5,width = 8)

####################################

kegg <- enrichKEGG(gene = eg$ENTREZID,
                   organism = 'hsa', #KEGG可以用organism = 'hsa'
                   pvalueCutoff = 0.05)
head(kegg,2)



k1<-barplot(kegg,drop=TRUE,showCategory=10)
k2<-k1+theme(axis.title=element_text(size=15,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=15,family="serif",color="black"),
             axis.text.y=element_text(size=15,family="serif",color="black"),
             legend.text=element_text(size=15,family="serif",color="black"),
             legend.title=element_text(size=15,family="serif",color="black",face="bold")
)

k3<-k2+ggtitle("TCGA")+
  xlab("Count")+
  ylab("KEGG pathway")+
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
k3
ggsave(k3,filename = "master_regulon_KEGG.pdf",height = 5,width = 8)


####################################

sig_TF<-data.frame(t(BLCA[sig,]))

data<-merge(sur.cat,sig_TF,by = "row.names")


cox_TF<-matrix(nrow=ncol(data),ncol=6)
for (i in 6:length(data)) {
  Bcox<-coxph(Surv(times, status)~data[,i],data=data)
  summcph<-summary(Bcox)
  cox_TF[i,1]<-summcph$conf.int[1]
  cox_TF[i,2]<-summcph$conf.int[3]
  cox_TF[i,3]<-summcph$conf.int[4]
  cox_TF[i,4]<-as.matrix(summcph$logtest)[3]
  cox_TF[i,5]<-as.matrix(summcph$sctest)[3]
  cox_TF[i,6]<-summcph$coefficients[5]
  print(i)
}
rownames(cox_TF)<-colnames(data)
cox_TF<-na.omit(cox_TF)

colnames(cox_TF)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")

save(cox_TF,file = "cox_TF.rda")

#############################################
load("cox_TF.rda")

cox_TF<-data.frame(cox_TF[which(cox_TF[,6]<0.01),])

hr=c(sprintf("%0.2f",as.numeric(cox_TF$HR)))
hr1=format(hr,digits=3)
low=c(round(as.numeric(cox_TF$Lower.95),2))
low1=format(low,digits = 3)
low2=gsub(" ","", low1)
upp=c(round(as.numeric(cox_TF$Upper.95),2))
upp1=format(upp,digits = 3)
upp2=gsub(" ","", upp1)

pvalue=c(as.numeric(cox_TF$p_value))
value=format(pvalue,scientific = T,digits = 3)
HR=paste0(hr1,"(",low2,"-",upp2,")")
set<-rownames(cox_TF)
dat=cbind(c("Feature",set),c("HR (95% CI)",HR),c("P-value",value))
datA<-dat[-1,]
datA<-data.frame(datA)

HR_cox<-data.frame(cox_TF)
HR_cox$Group<-row.names(HR_cox)
HR_cox$var<-row.names(HR_cox)
HR_cox$CI<-datA$X2
HR_cox$P.value<-datA$X3
HR_cox<-HR_cox[order(HR_cox$HR),]
HR_cox$gene<-factor(rownames(HR_cox),levels = rownames(HR_cox))


####################################
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(cowplot)

HR_cox1<-HR_cox[which(HR_cox$HR>1),]

cbbPalette<-colorRampPalette(brewer.pal(12,"Set3"))(nrow(HR_cox1))

p1<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$gene,family="serif",color=c("black"),
            hjust=0,fontface="bold",inherit.aes = FALSE,size=4) +
  ggtitle("Master regulon")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.065))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()
  )

p1<-p1+scale_fill_manual(values=cbbPalette)
p1


p2<-ggplot(HR_cox1,aes(gene,HR,fill =gene)) +
  theme(panel.background = element_rect(fill='transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black',size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face="bold",colour = "black",hjust = 0.5)) +
  coord_flip() +
  xlab("") +
  ylab("") +
  #labs(title="95% confidence intervals")+
  ggtitle("95% confidence intervals")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))


p2<-p2+geom_errorbar(aes(ymin=Lower.95,ymax=Upper.95,color=gene), 
                     width=0.5,size = 0.8) +
  geom_point(aes(color=gene),shape=22,size=6)+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  geom_hline(aes(yintercept=median(HR_cox1$HR)),linetype='dashed')+
  theme(#y轴刻度内容调整
    axis.text.x=element_text(size=13,family="serif"),
  )

p2 

p3<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$P.value,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("P-value")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p3


p4<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$CI,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("HR (95% CI)")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p4

p<-p1+p2+p4+p3+plot_layout(widths=c(2,5,2.5,1.5))

forest<-p
forest

save(forest,file = "forestplot_cox1.rda")

pdf('forestplot_cox1.pdf',width = 10,height = 8)
forest
dev.off()

################################

HR_cox1<-HR_cox[which(HR_cox$HR<1),]

cbbPalette<-colorRampPalette(brewer.pal(12,"Set3"))(nrow(HR_cox1))

p1<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$gene,family="serif",color=c("black"),
            hjust=0,fontface="bold",inherit.aes = FALSE,size=4) +
  ggtitle("Master regulon")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.065))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()
  )

p1<-p1+scale_fill_manual(values=cbbPalette)
p1


p2<-ggplot(HR_cox1,aes(gene,HR,fill =gene)) +
  theme(panel.background = element_rect(fill='transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black',size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face="bold",colour = "black",hjust = 0.5)) +
  coord_flip() +
  xlab("") +
  ylab("") +
  #labs(title="95% confidence intervals")+
  ggtitle("95% confidence intervals")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))


p2<-p2+geom_errorbar(aes(ymin=Lower.95,ymax=Upper.95,color=gene), 
                     width=0.5,size = 0.8) +
  geom_point(aes(color=gene),shape=22,size=6)+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  geom_hline(aes(yintercept=median(HR_cox1$HR)),linetype='dashed')+
  theme(#y轴刻度内容调整
    axis.text.x=element_text(size=13,family="serif"),
  )

p2 

p3<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$P.value,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("P-value")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p3


p4<-ggplot(HR_cox1,aes(HR_cox1$gene,HR_cox1$HR,fill=gene))+
  geom_text(aes(y=0,x=gene),label=HR_cox1$CI,family="serif",
            hjust=0,fontface = "bold",inherit.aes = FALSE,size=4) +
  ggtitle("HR (95% CI)")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))+
  coord_flip()+
  ylim(c(0,1))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())
p4

p<-p1+p2+p4+p3+plot_layout(widths=c(2,5,2.5,1.5))

forest<-p
forest

save(forest,file = "forestplot_cox2.rda")

pdf('forestplot_cox2.pdf',width = 10,height = 8)
forest
dev.off()