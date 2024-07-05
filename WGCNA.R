#1.加载包

rm(list=ls())
library(WGCNA)
library(data.table)
library(stringr)
library(gplots)
library(Biobase)
allowWGCNAThreads()

setwd("D:\\research\\metabolic\\WGCNA")

load("D:/research/data_new/Bladder_7.rda")
load("D:/research/tf_activities/tf_activities_7.rda")
load("D:/research/metabolic/WGCNA/k-mean.rda")


exp<-t(result[["TCGA-BLCA_FPKM"]])
pdata<-pData(Bladder_7[["TCGA-BLCA_FPKM"]])
pdata1<-pdata[,c(4,5,7)]
pdata1<-na.omit(pdata1)
pdata2<-merge(pdata1,label,by = "row.names")
rownames(pdata2)<-pdata2[,1]
pdata2<-pdata2[,-1]
colnames(pdata2)<-c("Age","Gender","Stage","Cluster")
exp<-exp[rownames(pdata2),]
# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 获得各个阈值下的 R方 和平均连接度
sft = pickSoftThreshold(exp, powerVector = powers, verbose = 5)

pdf("1.软阈值筛选.pdf",height = 6,width = 10)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



#4.根据β值获得临近矩阵和拓扑矩阵

# 获得临近矩阵：
softPower <- sft$powerEstimate
adjacency = adjacency(exp, power = softPower)
# 将临近矩阵转为 Tom 矩阵
TOM = TOMsimilarity(adjacency)
# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")

pdf("2.检验选定的β值.pdf",height = 6,width = 10)

ADJ1_cor <- abs(WGCNA::cor( exp,use = "p" ))^softPower
# 基因少（<5000）的时候使用下面的代码：
k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
# # 基因多的时候使用下面的代码：
# k <- softConnectivity(datE=exp,power=softPower) 

par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
#可以看出k与p(k)成负相关(相关性系数0.9),说明选择的β值能够建立基因无尺度网络


#5.一步法构建共表达矩阵
net = blockwiseModules(
  exp,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  #saveTOMs = TRUE,
  #saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)

#6. 模块可视化
pdf("3.基因聚类模块.pdf",height = 6,width = 10)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.
dev.off()

pdf("0.样本聚类.pdf",height = 6,width = 10)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(exp), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
dev.off()

#7.模块和性状的关系
design<-pdata2
a<-model.matrix(~0+pdata2[,2])
b<-model.matrix(~1+pdata2[,3])
design[,2]<-a[,2]
design[,3]<-b[,2]+b[,3]*2+b[,4]*3


moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(exp))

#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


pdf("4.模块和性状的关系.pdf")
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),####
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# 首先计算模块与基因的相关性矩阵
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(exp, MEs, use = "p"))
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exp)))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


# 再计算性状与基因的相关性矩阵


Cluster<- as.data.frame(design[,4])
names(Cluster)<-"Subtype"
geneTraitSignificance = as.data.frame(cor(exp, Cluster, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exp)))
names(geneTraitSignificance) = paste("GS.", names(Cluster), sep="")
names(GSPvalue) = paste("p.GS.", names(Cluster), sep="")

pdf("5.grey_Cluster.pdf",height = 6,width = 6)

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "grey"#更改模块
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),#选择对应性状的相关性列
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Drug sensitivity for Cluster score",#更改该性状名称
                   main = paste("Module membership vs. Drug sensitivity\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

###############
pdf("5.blue_Cluster.pdf",height = 6,width = 6)

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "blue"#更改模块
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),#选择对应性状的相关性列
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Drug sensitivity for Cluster score",#更改该性状名称
                   main = paste("Module membership vs. Drug sensitivity\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##################

pdf("5.green_Cluster.pdf",height = 6,width = 6)

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "green"#更改模块
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),#选择对应性状的相关性列
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Drug sensitivity for Cluster score",#更改该性状名称
                   main = paste("Module membership vs. Drug sensitivity\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
#9.网络的可视化
# 首先针对所有基因画热图

geneTree = net$dendrograms[[1]]
dissTOM = 1-TOMsimilarityFromExpr(exp, power = sft$powerEstimate)
plotTOM = dissTOM^7
diag(plotTOM) = NA

pdf("6.allgene_network.pdf",height = 6,width = 6)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()


#10.提取指定模块的基因名
# Select module
module = "grey";
# Select module probes
probes = colnames(exp) ## 我们例子里面的probe就是基因名
# inModule = (moduleColors==module);
# modProbes = probes[inModule];
modProbes = cbind(probes,moduleColors)
write.table(modProbes,"8.模块基因.txt",quote = F)
# modProbes

grey_gene<-modProbes[modProbes[,2] %in% c("grey","blue","green"),]
write.table(grey_gene,"显著模块基因.txt",quote = F)
