
############################
rm(list=ls())

library(GSVA)
library(Biobase)

setwd("D:\\research\\tf_activities\\master")

load("D:/research/tf_activities/tf_activities_7.rda")
load("cox_TF.rda")
load("D:\\research\\data_new\\Bladder_7.rda")

exp<-Bladder_7[["TCGA-BLCA_FPKM"]]
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

cell_gene<-read.csv("immune.csv",header = T)
cell_gene <- cell_gene[,c(1,6)]
cell<-unique(as.character(cell_gene[,2]))
cellgene_type=list()
for (i in 1:length(cell)){
  b=as.character(cell_gene[(as.character(cell_gene$Category))==cell[i],1])
  cellgene_type[i]=list(b)
  names(cellgene_type)[i]=cell[i]
}

dat<-exprs(exp)
dat<-t(scale(t(dat)))
gsva_es<-gsva(dat,cellgene_type,method="ssgsea",abs.ranking=F)
z<- zClust(gsva_es)
NES<-z$data


##############################################
cox_TF<-data.frame(cox_TF[which(cox_TF[,6]<0.01),])

a<-result[["TCGA-BLCA_FPKM"]]

a<-a[rownames(cox_TF),]
b<-t(NES[,colnames(a)])

cor<-matrix(0,nrow = nrow(a),ncol = ncol(b))
for (i in 1:nrow(a)) {
  for (j in 1:ncol(b)) {
    A1<-cor(a[i,],b[,j],method = "pearson")
    cor[i,j]<-A1
  }
}
rownames(cor)<-rownames(a)
colnames(cor)<-colnames(b)

library(corrplot)
library(ggplot2)

pdf(file = "cor_immport_sig.pdf",width = 12,height = 7)
par(family="serif")
corrplot(t(cor),method="pie",is.corr=FALSE, tl.col = "black",mar = c(0,0,1.5,0),
         cl.pos = 'b',tl.cex=0.7,cl.ratio=0.15,cl.length=6,cl.cex=0.7,
         col= rev(RColorBrewer::brewer.pal(6,"RdBu")))

dev.off()

##########################################

metabolic<-read.table("metabolic_gene-114.txt",sep = "\t",fill = T)
metabolic<-metabolic[,c(2,1)]
me<-unique(as.character(metabolic[,1]))
metabolic_type=list()
for (i in 1:length(me)){
  b=as.character(metabolic[(as.character(metabolic[,1]))==me[i],2])
  metabolic_type[i]=list(b)
  names(metabolic_type)[i]=me[i]
}

gsva_es<-gsva(dat,metabolic_type,method="ssgsea",abs.ranking=F)
z<- zClust(gsva_es)
NES<-z$data

c<-t(NES[,colnames(a)])

cor1<-matrix(0,nrow = nrow(a),ncol = ncol(c))
for (i in 1:nrow(a)) {
  for (j in 1:ncol(c)) {
    A1<-cor(a[i,],c[,j],method = "pearson")
    cor1[i,j]<-A1
  }
}
rownames(cor1)<-rownames(a)
colnames(cor1)<-colnames(c)

library(corrplot)
library(ggplot2)

pdf(file = "cor_metabolic_sig.pdf",width = 12,height = 20)
par(family="serif")
corrplot(t(cor1),method="pie",is.corr=FALSE, tl.col = "black",mar = c(0,0,1.5,0),
         cl.pos = 'b',tl.cex=0.7,cl.ratio=0.15,cl.length=6,cl.cex=0.7,
         col= rev(RColorBrewer::brewer.pal(6,"RdBu")))

dev.off()