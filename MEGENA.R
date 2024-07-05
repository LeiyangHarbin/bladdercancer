
##########################
rm(list = ls())
library(MEGENA)
library(Biobase)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\TF activity\\")
load("tf_activities_7.rda")
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\MEGENA-new\\")

TF<-result[["TCGA-BLCA_FPKM"]]
counts<-unique(TF)


n.cores <- 15; # number of cores to be used for parallelization (need to MATCH those requested in your LSF).
doPar <- TRUE; # TRUE = perform parallel PFN -> MCA. FALSE = no parallelization
method = "pearson" # Currently “pearson” (Pearson’s correlation) and “spearman” (Spearman’s correlation available) for correlation.
FDR <- 0.05; # FDR threshold to identify significant interactions.
n.cor.perm = 100; # Number of permutations in correlation screening.
n.hub.perm = 100; # Number of permutations to calculate hub significance in MHA
mod.pval = 0.05; # module p-value threshold in MCA
hub.pval = 0.05; # hub p-value threshold in MHA
min.size = 10 # minimum module size
max.size = NULL # maximum module size
# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2


##########################Generate Correlation Matrix

ijw <- calculate.correlation(counts,doPerm = n.cor.perm,output.corTable = FALSE,output.permFDR = FALSE)


###########################Register multiple cores if needed: note that set.parallel.backend() is deprecated.
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}


####################################Major Heavy Memory Calculation (will take a lot of time)
# Compute PFN
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores)
# graph dataframe
g <- graph.data.frame(el,directed = FALSE)

## perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
                           mod.pval = mod.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 10,max.size = vcount(g)/2,
                           doPar = doPar,num.cores = n.cores,n.perm = n.hub.perm,
                           save.output = FALSE)

# unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


####################Summary of MEGENA results

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = mod.pval,hub.pvalue = hub.pval,
                                       min.size = min.size,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)

module.output <- module_convert_to_table(MEGENA.output,mod.pval = 0.05,
                                         hub.pval = 0.05,min.size = 10,max.size=vcount(g)/2)

save(summary.output,MEGENA.output,module.output,g,file="MEGENA.Results.rda")

#####################Gene Ontology Enrichment Analysis for Each Module



##########################
rm(list=ls())
library(MEGENA)
library(Biobase)
library(GOstats)
library(HGNChelper)
library(org.Hs.eg.db)
library(DGCA)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\TF activity\\")
load("tf_activities_7.rda")
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\MEGENA-new\\")
load("MEGENA.Results.rda")
TF<-result[["TCGA-BLCA_FPKM"]]
counts<-unique(TF)
# run ontology using vectors of gene symbols and module IDs
moduleGO_res = moduleGO(
  genes = module.output$id, 
  labels = module.output$module,
  universe = rownames(counts), 
  pval_GO_cutoff = 0.05)

ontology = extractModuleGO(moduleGO_res)

ontology <- ontology[,-length(colnames(ontology))]

### Filter top gene ontology for each module

pvalues <- ontology[,c(1:4,grep("pVal",colnames(ontology)))]

## Get p top Ontology terms output from MEGENA


## Reaad in ontology p value file if needed

#ontology <- read.csv("module_ontology_pvalues.csv",header=T)

## make it just whatever process you care about (e.g. biological process)

pvalues <- pvalues[pvalues$Ontology=="BP",]
modules <- colnames(pvalues[5:length(colnames(pvalues))])
top.ontology <- data.frame()
for (module in modules) {
  temps <- pvalues[order(pvalues[module]),]
  temp <- temps[1,1:4]
  x <- temps[module][[1]]
  x <- x[1]
  temp$module <- module
  temp$pval <- x
  top.ontology <- rbind(top.ontology,temp)
}
top.ontology$module <- unique(module.output$module)
write.csv(top.ontology,file="top.ontology.per.module.csv")
save(top.ontology,file="top.ontology.rda")

##########################



#####################Compute modules enriched for DEGs
rm(list=ls()) 
library(GeneOverlap)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\MEGENA-new\\")
load("MEGENA.Results.rda")
load("sig.rda")

### Generate List of Module Lists and determine enrichment of DEGs
### Need module.IDs and two column "modules" dataframe with module and genes
modules<-module.output[,c(1,6)]
module_list = list()
for (module in modules$module){
  k <- modules[modules$module==module,]
  ids <- as.vector(k$id)
  module_list[[module]] <- ids
}
# load in GeneOverlap 

## DEGs is a list of your DEGs
## DEGs should only be genes in your modules
DEGs <- list()
degs <- intersect(rownames(sig),modules$id)
DEGs[["DEGs"]] <- degs
total.genes <- length(unique(modules$id)) ## Make this the total number of genes across all modules
Object <- newGOM(DEGs,module_list,total.genes)
overlap.intersect <- getMatrix(Object, name="intersection")
overlap.pval <- getMatrix(Object, name="pval")
overlap.OR <- getMatrix(Object, name="odds.ratio")
Overlap.sum <- cbind(overlap.intersect,overlap.pval,overlap.OR)

#head(Overlap.sum)
#        overlap.intersect overlap.pval overlap.OR
#c1_2.ok                 0 1.0000000000    0.00000
#c1_3.ok                 1 0.0473062048   24.20475

Overlap.sum <- data.frame(Overlap.sum)
write.csv(Overlap.sum,file="Enrichment.for.DEGs.in.Modules.csv")
save(Overlap.sum,file="Enrichment.for.DEGs.in.Modules.rda")
##############################Visualize significant modules



# Compile module stats
top.ontology<-read.csv("top.ontology.per.module.csv",header = T,row.names=1)

SUM<-cbind(as.character(top.ontology$module),as.numeric(top.ontology$Size),as.character(top.ontology$Term),Overlap.sum$overlap.pval,as.numeric(Overlap.sum$overlap.OR))
SUM<-data.frame(SUM)
colnames(SUM) <- c("module","size","term","pval","OR")
SUM$pval <- -log10(as.numeric(as.character(SUM$pval)))
SUM$OR <- as.numeric(SUM$OR)
SUM$size <- as.numeric(SUM$size)
SUM$term <- as.character(SUM$term)

library(ggplot2)
ggplot(SUM,aes(x=OR,y=pval,size=size,label=term))+
  geom_point(alpha=.25,shape=21,colour="black",fill="blue",stroke=2)+
  geom_text(check_overlap=FALSE,aes(label=ifelse(pval>3,term,""),size=1))+ #,fontface="bold"))+ #hjust=ifelse(pval>3,.1,.3315),vjust= ifelse(OR<60,.6,-.4),fontface="bold"))+
  scale_size_continuous(range = c(3, 20),breaks=c(10,100,1000))+
  ylab("-log10(pvalue)")+
  xlab("Odds Ratio")+
  labs(size="Module Size")+
  theme(
    panel.background = element_rect(colour="black",fill="white"),
    legend.title=element_text(size=20,family="serif"),
    legend.text=element_text(size=18,family="serif"),
    axis.text.x=element_text(size=22,face="bold",family="serif"),
    axis.text.y=element_text(size=22,face="bold",family="serif"),
    axis.title.x=element_text(size=22,face="bold",family="serif"),
    axis.title.y = element_text(size=22,face="bold",family="serif"))






############################Determine Module Eigengene-Trait Relationships
rm(list=ls())
library(WGCNA) ## Eigengene calculation uses WGCNA code
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\TF activity\\")
load("tf_activities_7.rda")
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BLCA-dorothea\\MEGENA-new\\")
load("pdata.rda")
load("MEGENA.Results.rda")
load("TCGA_group.rda")
TF<-result[["TCGA-BLCA_FPKM"]]
counts<-unique(TF)
counts<-t(counts)
counts<-counts[rownames(pdata2),]
modules <- module.output[,c(1,6)]
sur.cat<-sur.cat[rownames(pdata2),]
pdata2[,4]<-sur.cat$risk

design=model.matrix(~0+pdata2[,3])
design<-cbind(design,as.numeric(pdata2[,4]))
colnames(design)=c("StageI","StageII","StageIII","StageIV","Risk_score")


nSamples<-nrow(counts)## Number of subjects
moduleColors <- modules$module ## A vector of module ids corresponding to every gene in counts (note counts need to be genes in columns, subjects in rows)
## Will need counts dataframe to feature column for gene corresponding to each module it is in
expression.table<-counts[,as.vector(modules$id)] ## counts should already be transposed before this
MEs0 = moduleEigengenes(expression.table, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design, use = "p",method="spearman") ## Probably use spearman correlation to control for outliers, both MEs and metadata uses subjects for rows
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("模块和性状的关系.pdf",height = 20)
par(mar = c(6, 8.5, 3, 3),family="serif");
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




#################挑选risk_score相关性排名前十的模块画图

moduleTraitCor1<-moduleTraitCor[order(abs(moduleTraitCor[,5]),decreasing = T),]
moduleTraitCor1<-moduleTraitCor1[c(1:10),]
moduleTraitPvalue1<-moduleTraitPvalue[rownames(moduleTraitCor1),]
MEs1<-MEs[,rownames(moduleTraitCor1)]

textMatrix = paste(signif(moduleTraitCor1, 2), "\n(",
                   signif(moduleTraitPvalue1, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor1)

pdf("模块和性状的关系top10.pdf")
par(mar = c(6, 8.5, 3, 3),family="serif");


labeledHeatmap(Matrix = moduleTraitCor1,
               xLabels = colnames(design),####
               yLabels = names(MEs1),
               ySymbols = names(MEs1),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#############################################
adjusted.p <- data.frame(1:67) ## long winded way of getting FDR-adjusted pvalues
for (i in colnames(moduleTraitPvalue)){
  adjust <- p.adjust(moduleTraitPvalue[,i],method="BH") ## can choose other methods like Bonferroni
  adjusted.p <- cbind(adjusted.p,adjust)
}
adjusted.p <- adjusted.p[,-1]
colnames(adjusted.p)<-colnames(moduleTraitPvalue)


##################################Determine Modules differentially connected using DGCA
counts <- t(counts) 
load("pdata.rda")

design1=model.matrix(~0+pdata2[,4])
colnames(design1)<-c("high","low")
metas = data.matrix(design1)   ## metas is a column for your group with a 0 or 1

#Control <- ""  ## Use name of column
#Treatment <- "" ## Use name of column
compare = c("high","low")

modules <- modules 
  
library(DGCA)
diff.connect<-moduleDC(counts,metas, compare, modules$id, modules$module, 
                         corr_cutoff = 0.99,
                         signType = "none", 
                         corrType = "pearson", 
                         nPerms = 50,
                         oneSidedPVal = FALSE, 
                         gene_avg_signif = 0.05, 
                         number_DC_genes = 3,
                         dCorAvgMethod = "median")

write.csv(diff.connect,file="Differential.module.connectivity.analysis.csv")
save(diff.connect,file="Differential.module.connectivity.analysis.rda")

##################Extracting weights of gene-gene relationships for cytoscape

pnet.obj <- plot_module(output = summary.output,PFN = g,subset.module = "comp1_3",
                        layout = "kamada.kawai",label.hubs.only = FALSE,
                        gene.set = NULL,color.code =  "grey",
                        output.plot = FALSE,out.dir = "modulePlot",col.names = c("magenta","green","cyan"),label.scaleFactor = 2,
                        hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)