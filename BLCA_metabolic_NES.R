################################
rm(list = ls())
library(GEOquery)
library(GSVA)
#############################################
distCor <- function(x) as.dist(1-cor(x))
zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
  if (scale=="row") z <- t(scale(t(x))) else z<-x
  if (scale=="col") z <- scale(x)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  hcl_row <- hclust(distCor(t(z)), method=method)
  hcl_col <- hclust(distCor(z), method=method)
  return(list(data=z, Rowv=as.dendrogram(hcl_row), Colv=as.dendrogram(hcl_col)))
}

###############################################
setwd("D:\\research\\metabolic")
load("D:/research/data_new/Bladder_7.rda")

library(clusterProfiler)
library(tidyr)

metabolic<-read.table("metabolic_gene-114.txt",sep = "\t",fill = T,header = F)
colnames(metabolic)<-c("gene","term")
term<-unique(as.character(metabolic[,2]))
metabolic_type=list()
for (i in 1:length(term)){
  b=as.character(metabolic[(as.character(metabolic$term))==term[i],1])
  metabolic_type[i]=list(b)
  names(metabolic_type)[i]=term[i]
}


#################对osdata计算每个数据集的ssGSEA
BLCA_metabolic=list()
for (k in 1:length(Bladder_7)) 
{
  MyGSE<-Bladder_7[[k]]
  dat<-exprs(MyGSE)
  dat<-t(scale(t(dat)))
  gsva_es<-gsva(dat,metabolic_type,method="ssgsea",abs.ranking=F)#计算免疫细胞的富集得分
  z<- zClust(gsva_es)
  NES<-z$data
  BLCA_metabolic[k]=list(NES)
  names(BLCA_metabolic)[k]=names(Bladder_7)[k]
  print(k)
}

save(BLCA_metabolic,file="BLCA_metabolic_NES.rda")

