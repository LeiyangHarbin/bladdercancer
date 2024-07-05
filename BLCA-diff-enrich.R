
###########################################
rm(list=ls())
library("edgeR")
setwd("D:\\reserach1\\enrich-diff")
load("TCGA-BLCA_RNASeq_Counts.rda")
load("pdata.rda")


BLCA<-expDataCounts
count<-BLCA[apply(BLCA, MARGIN=2,FUN=function(xxx)
{
  (sum(xxx==0)/length(xxx))<=0.5
}),]

high<-pdata2[pdata2$Risk_score == "high",]
low<-pdata2[pdata2$Risk_score == "low",]

group_list = factor(c(high$Risk_score,low$Risk_score))
design <- model.matrix(~0+group_list)
rownames(design) = colnames(data)
colnames(design) <- levels(group_list)

data<-count[,c(rownames(high),rownames(low))]


#差异表达矩阵
DGElist <- DGEList( counts = data, group = group_list)
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 ## 自定义
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1)) 
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
save(nrDEG_edgeR,file = 'TCGA_BLCA_edgeR.rda')
head(nrDEG_edgeR)


#提取差异显著的差异矩阵
padj = 0.05# 自定义
foldChange= 1.0 # 自定义
nrDEG_edgeR_signif  = nrDEG_edgeR[(nrDEG_edgeR$FDR < padj & 
                                     (nrDEG_edgeR$logFC>foldChange | nrDEG_edgeR$logFC<(-foldChange))),]
nrDEG_edgeR_signif = nrDEG_edgeR_signif[order(nrDEG_edgeR_signif$logFC),]

save(nrDEG_edgeR_signif,file = 'TCGA_BLCA_edgeR_signif.rda')

########################################

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)

load("TCGA_BLCA_edgeR_signif.rda")

x<-rownames(nrDEG_edgeR_signif)
eg <- bitr(x, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")


#go
go <- enrichGO(gene = eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID',
               readable = T)
head(go)


library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(3,"Dark2")
colorl <- rep(color,each=10)

pdata<-go@result
pdata<-pdata[,c(1,3,10)]
count1<-pdata[which(pdata$ONTOLOGY=="BP"),]
pdata1<-count1[order(count1$Count,decreasing = T),][1:10,]
count2<-pdata[which(pdata$ONTOLOGY=="CC"),]
pdata2<-count2[order(count2$Count,decreasing = T),][1:10,]
count3<-pdata[which(pdata$ONTOLOGY=="MF"),]
pdata3<-count3[order(count2$Count,decreasing = T),][1:10,]
pdata4<-rbind(pdata1,pdata2,pdata3)
pdata4$Description <- factor(pdata4$Description,levels=pdata4$Description)

p1<-ggplot(pdata4) +
  aes(x = Description, y = Count, fill = ONTOLOGY) +
  geom_bar(stat = "identity",colour="black") +
  xlab("GO term")+
  coord_flip()+
  scale_fill_manual(values =color)+
  theme(
    axis.title=element_text(size=15,face="bold",family = "serif",color="black"),
    axis.text = element_text(size=12,face="plain",family = "serif",color="black"),
    axis.text.y = element_text(colour = colorl,hjust=1,vjust=0.6),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold",family = "serif"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.direction = "vertical",
    legend.background = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "black"),
    plot.background = element_blank()
  )
p1
ggsave(p1,filename = "BLCA_diff_GO.pdf",height = 5,width = 8)
##########################################
kegg <- enrichKEGG(gene = eg$ENTREZID,
                   organism = 'hsa', #KEGG可以用organism = 'hsa'
                   pvalueCutoff = 0.05)
head(kegg,2)


pdf("BLCA_diff_KEGG.pdf",height = 5,width = 8)
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
dev.off()