

############################################
rm(list= ls())

library(MEGENA)
library(GEOquery)
library(WGCNA)

setwd("D:\\research\\MEGENA\\BLCA\\enrich_diffmodules")

load("sig_modules.rda")
load("D:/research/MEGENA/BLCA/enrich_diffmodules/megena_res_BLCA.rda")

sig_modules<-substring(sig_modules,3)

modulehub<-megena_res[["full"]][["hub.output"]][["module.degreeStat"]]


hub<-data.frame(matrix(NA,ncol = 4,nrow = 1))
colnames(hub)<-colnames(modulehub[[1]])
for (i in sig_modules) {
  a<-modulehub[[i]]
  a<-a[order(a$degree,decreasing = T),]
  a<-a[1,]
  hub<-rbind(hub,a)
}

hub<-hub[-1,]


########################

library(ggplot2)

hub<-hub[order(hub$degree,decreasing = F),]
hub$gene<-factor(hub$gene,levels = rownames(hub))

p1<-ggplot(hub,aes(x = degree ,y = gene))+
  geom_segment(aes(x = 0 ,xend = degree,y = gene,yend = gene ))+
  geom_point(color = "darkorange",size = 3)+
  labs(size = "Degree",color = "P-value",y = "Gene", x = "Degree of genes in modules")+
  theme(legend.title=element_text(size=12,family="serif",color="black",face="bold"),
        legend.text=element_text(size=8,family="serif",color="black",face="bold"),
        axis.line = element_line(size=1, colour = "black"),
        panel.background = element_rect(fill="white",colour="white",size=0.5,linetype="solid"),
        axis.text.x=element_text(colour="black", size = 12,family="serif",face= "bold"),
        axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
        axis.text.y=element_text(colour="black", size = 12,family="serif",face= "bold"),
        
        axis.line.x = element_line(size=0.5, colour = "black"),
        axis.line.y = element_line(size=0.5, colour = "black"))
p1

ggsave(p1,filename = "degree_segment.pdf",width = 7,height = 8)
