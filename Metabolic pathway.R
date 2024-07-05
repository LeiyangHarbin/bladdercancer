
#################################################################
rm(list=ls())

library(survminer)
library(survival)
library(survcomp)

setwd("D:\\research\\metabolic")
load("D:/research/metabolic/BLCA_metabolic_NES.rda")
load("risk_TCGA.rda")
TCGA<-BLCA_metabolic[["TCGA-BLCA_FPKM"]]
TCGA<-data.frame(t(TCGA))
risk_TCGA$sample<-gsub("\\.","-",rownames(risk_TCGA))
TCGA$sample<-rownames(TCGA)
risk_TCGA<-risk_TCGA[,c(2,3,5)]
data<-merge(risk_TCGA,TCGA,by = "sample")
rownames(data)<-data[,1]
data1<-data[,-1]

cox<-matrix(nrow=ncol(data1),ncol=6)
for (i in 3:length(data1)) {
  Bcox<-coxph(Surv(times, status)~as.numeric(as.character(data1[,i]))>median(as.numeric(as.character(data1[,i]))),data=data1)
  summcph<-summary(Bcox)
  cox[i,1]<-summcph$conf.int[1]
  cox[i,2]<-summcph$conf.int[3]
  cox[i,3]<-summcph$conf.int[4]
  cox[i,4]<-as.matrix(summcph$logtest)[3]
  cox[i,5]<-as.matrix(summcph$sctest)[3]
  cox[i,6]<-summcph$coefficients[5]
  print(i)
}
rownames(cox)=colnames(data1)
colnames(cox)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
cox<-cox[-c(1,2),]
cox<-data.frame(cox)

#########

library(factoextra)

data<-TCGA[,rownames(cox)]

p=fviz_nbclust(data, kmeans, method = "silhouette",linecolor = "darkred")
p1<-p+geom_vline(xintercept = 2, linetype = 2,col="blue")
p1=p1+theme(axis.title.x = element_text(size = 15,family = "serif",face = "bold"),
            axis.title.y = element_text(size = 15,family = "serif",face = "bold"),
            axis.text.x = element_text(size = 12,family = "serif",face = "bold"),
            axis.text.y = element_text(size =12,family = "serif",face = "bold")
            ,title=element_text(size = 18,family = "serif",face = "bold"))
p1=p1+theme(plot.title = element_text(hjust = 0.5))
p1
ggsave(p1,filename = "K-mean.pdf")

km<- kmeans(data, 2, nstart = 2)
label=km$cluster
label=data.frame(label)


#############################

library(ComplexHeatmap)
library(circlize)

heatmap<-data.frame(BLCA_metabolic[["TCGA-BLCA_FPKM"]])
label$sample<-rownames(label)
label<-label[order(label[,1]),]
label[label == 1,1] <-"Cluster 1"
label[label == 2,1] <-"Cluster 2"
rownames(label)<-gsub("-","\\.",rownames(label))
heatmap<-heatmap[,rownames(label)]

col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("#00AFFF", "white", "#FC0E00")
)



col_an<-columnAnnotation(Subtype = factor(label[,1]), 
                         col = list(Subtype = c("Cluster 1" = "#FC0E00", "Cluster 2" ="#00AFFF")),
                         annotation_legend_param = list(title = "Subtype",
                                                        title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                                                        labels_gp = gpar(fontsize = 10,fontfamily = "serif"),
                                                        grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")),
                         border = F,
                         annotation_name_side = "right",
                         annotation_name_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold")
)

pdf("Metabolic_pathway_heatmap.pdf",,width = 9,height = 12)
Heatmap(heatmap,
        col = col_fun,
        show_column_names = F,
        top_annotation = col_an,
        row_title = "Metabolic pathway",
        column_title = "TCGA",
        column_title_side = "top",
        column_split = label[,1],
        column_title_gp = gpar(fontsize = 17,fontfamily = "serif",fontface = "bold"),
        row_title_gp = gpar(fontsize = 17,fontfamily = "serif",fontface = "bold"),
        cluster_columns = F,
        cluster_rows = F,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8,fontfamily = "serif"),
        heatmap_legend_param = list(title = "NES",
                                    title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                                    labels_gp = gpar(fontsize = 10,fontfamily = "serif"),
                                    legend_height = unit(8, "cm"),
                                    grid_width = unit(0.8, "cm")),
        border = TRUE
)
dev.off()
