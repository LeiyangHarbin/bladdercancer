

######################################
rm(list=ls())
library(DGCA)
library(MEGENA)
library(doMC)

setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\DGCA-WGCNA\\")
load("exp.rda")
load("pdata.rda")

exp1<-t(exp)
exp1<-exp1[,rownames(pdata2)]
datExpr<-log(exp1+0.001)
counts<-unique(datExpr)


design=model.matrix(~0+pdata2[,4])
rownames(design)<-colnames(counts)
colnames(design)<-c("high","low")
ddcor_res=ddcorAll(inputMat = counts, design = design,sigThresh=0.05,
                   corr_cutoff = 0.99,corSigThresh=0.05,
                   compare = c("high","low"),
                   adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")



save(ddcor_res,file="ddcor_res_BLCA.rda")
#的相关性再进行ddMEGENA分析，这样就找到了2类之间富集的模块
#######################################################
n.cores<-16; # 内核/线程数
doPar <-TRUE; # 是否需要并行

#################################
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}
###############################################################
ddcor_res1<-na.omit(ddcor_res) 

megena_res=ddMEGENA(ddcor_res1, pval_gene_thresh = 0.05,adjusted =TRUE, 
                    nPerm = 100, hubPVal = 0.05,modulePVal = 0.05, 
                    minModSize = 10, maxModSize = 1000,
                    evalCompactness=TRUE,parallelize=TRUE,nCores=n.cores)



save(megena_res,file = "megena_res_BLCA.rda")

#####################################

#####################################

######################################
rm(list=ls())
library(DGCA)
library(MEGENA)
library(doMC)
library(WGCNA)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\DGCA-WGCNA\\")
load("exp.rda")
load("pdata.rda")

load("exp_high.rda")
# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 获得各个阈值下的 R方 和平均连接度
sft = pickSoftThreshold(exp_high, powerVector = powers, verbose = 5)

sft$powerEstimate<-ifelse(is.na(sft$powerEstimate == T),6,sft$powerEstimate)

softPower = sft$powerEstimate
adjacency = adjacency(exp_high, power = softPower)
# 将临近矩阵转为 Tom 矩阵
TOM = TOMsimilarity(adjacency)
# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")

net = blockwiseModules(
  exp_high,
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

high_net<-net

save(high_net,file = "high_net_BLCA.rda")

#######################################
load("exp_low.rda")

# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 获得各个阈值下的 R方 和平均连接度
sft = pickSoftThreshold(exp_low, powerVector = powers, verbose = 5)

sft$powerEstimate<-ifelse(is.na(sft$powerEstimate == T),6,sft$powerEstimate)

softPower = sft$powerEstimate
adjacency = adjacency(exp_low, power = softPower)
# 将临近矩阵转为 Tom 矩阵
TOM = TOMsimilarity(adjacency)
# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")

net = blockwiseModules(
  exp_low,
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
low_net<-net

save(low_net,file = "low_net_BLCA.rda")



#################################
load("exp.rda")
load("pdata.rda")

setLabels = c("high","low")
nSets = 2
# Object that will contain the expression data
multiExpr = list()
high<-rownames(pdata2[which(pdata2$Risk_score == "high"),])
low<-rownames(pdata2[which(pdata2$Risk_score == "low"),])
multiExpr[[1]] = list(data = exp[high,])
multiExpr[[2]] = list(data = exp[low, ])
names(multiExpr) = setLabels


colorList = list(high_net, low_net)
names(colorList) = setLabels


system.time( {
  mp = modulePreservation(multiExpr, colorList,
                          referenceNetworks = c(1:2),
                          loadPermutedStatistics = FALSE,
                          nPermutations =100,randomSeed = 12345,
                          verbose = 3)
} )

save(mp, file = "modulePreservation_BLCA.rda")

#####################################
###########################################
nSets=2

library(impute)

impExpr = list();
for (set in 1:nSets)
{
  impExpr[[set]] = list(data = t(impute.knn(t(multiExpr[[set]]$data))$data));
}
eigengenes = list();
for (set in 1:nSets)
{
  eigengenes[[set]] = multiSetMEs(impExpr, universalColors = colorList[[set]], excludeGrey = TRUE);
  for (ss in 1:nSets)
  {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExpr[[ss]]$data);
  }
}
# Here comes the IGP calculation
library(clusterRepro)
cr=list();
set.seed(20);
for (ref in 1:nSets)
{
  cr[[ref]] = list();
  for (test in 1:nSets)
  {
    printFlush(system.time({
      cr[[ref]][[test]] = clusterRepro(Centroids = as.matrix(eigengenes[[ref]][[test]]$data),
                                       New.data = as.matrix(impExpr[[test]]$data),
                                       Number.of.permutations = 100); }));
    collectGarbage();
  }
}
# Save the results
save(cr,file="clusterRepro_BLCA.rda")

#######################################