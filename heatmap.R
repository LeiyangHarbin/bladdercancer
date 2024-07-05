
############################
rm(list=ls())

library(ComplexHeatmap)
library(circlize)
library(Hmisc)

setwd("D:\\research\\tf_activities\\master")

load("cox_TF.rda")
load("D:/research/tf_activities/tf_activities_7.rda")

cox_TF<-data.frame(cox_TF[which(cox_TF[,6]<0.01),])
BLCA<-result[["TCGA-BLCA_FPKM"]]

dat<-BLCA[rownames(cox_TF),]

load("D:/research/tf_activities/master/TCGA_group.rda")


label<-sur.cat
label<-label[order(label$score),]
library(Hmisc)
label$score<-capitalize(label$score)

cluster_high<-gsub("\\.","-",rownames(label[label$score == "High",]))
cluster_low<-gsub("\\.","-",rownames(label[label$score == "Low",]))
viper_TCGA_high<-dat[,cluster_high]
viper_TCGA_low<-dat[,cluster_low]

########################################

Sigmark<-function(plot1){
  plot1[which(plot1[,1]>0.05|plot1[,1]==1),2]="NS"
  plot1[which(0.01<plot1[,1]&plot1[,1]<0.05),2]="*"
  plot1[which(0.001<plot1[,1]&plot1[,1]<0.01),2]="**"
  plot1[which(plot1[,1]<0.001),2]="***"
  return(plot1[,2])
}

TCGA_wilcox<-data.frame()
for(i in 1:nrow(dat)){
  TCGA_wilcox[i,1]<-wilcox.test(viper_TCGA_high[i,],viper_TCGA_low[i,],
                                exact=FALSE,correct=FALSE)$p.value
  TCGA_wilcox[i,2]<-mean(viper_TCGA_high[i,])
  TCGA_wilcox[i,3]<-mean(viper_TCGA_low[i,])
  if (TCGA_wilcox[i,2]>TCGA_wilcox[i,3]){
    TCGA_wilcox[i,4]<-"High"
  }else if (TCGA_wilcox[i,2]<=TCGA_wilcox[i,3]){
    TCGA_wilcox[i,4]<-"Low"
  }
}

TCGA_wilcox[,5]<-Sigmark(TCGA_wilcox)
row.names(TCGA_wilcox)<-row.names(dat)
colnames(TCGA_wilcox)<-c("P-value","High","Low","label","Signifi")
high<-subset(rownames(TCGA_wilcox),TCGA_wilcox[,4] == "High")
low<-subset(rownames(TCGA_wilcox),TCGA_wilcox[,4] == "Low")

######################################################



heatmap<-cbind(viper_TCGA_high,viper_TCGA_low)
heatmap<-heatmap[c(high,low),c(cluster_high,cluster_low)]

subtype<-c(rep("High",length(cluster_high)),rep("Low",length(cluster_low)))
col_fun <- colorRamp2(
  c(-7, 0, 7),
  c("#00AFFF", "white", "#FC0E00")
)

ha11=rowAnnotation(Pvalue=row_anno_text(as.matrix(TCGA_wilcox[rownames(heatmap),5]), rot =0,offset=unit(0, "mm")))

col_an<-columnAnnotation(Subtype = factor(label$score), 
                         col = list(Subtype = c("High" = "#FC0E00", "Low" ="#00AFFF")),
                         annotation_legend_param = list(title = "Subtype",
                                                        title_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold"),
                                                        labels_gp = gpar(fontsize = 10,fontfamily = "serif",fontface = "bold"),
                                                        grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")),
                         border = F,
                         annotation_name_side = "right",
                         annotation_name_gp = gpar(fontsize = 13,fontfamily = "serif",fontface = "bold")
)

ha_cn1=Heatmap(heatmap,row_names_side="left",
               cluster_rows = F,cluster_columns =F,
               top_annotation=col_an,
               show_column_names=F,show_row_names =T,col=col_fun,
               row_split = rep(c("high","low"),c(length(high),length(low))),
               column_split = label$score,
               row_title="Master regulon",column_title="TCGA",
               row_title_gp=gpar(fontsize=15,fontface = "bold",fontfamily="serif"),
               row_names_gp = gpar(fontsize=8,fontface = "bold",fontfamily="serif"),
               column_title_gp=gpar(fontsize=15,fontface = "bold",fontfamily="serif"),
               heatmap_legend_param=list(title="Viper",
                                         title_gp = gpar(fontsize = 11,fontface = "bold",fontfamily="serif"),
                                         labels_gp = gpar(fontsize = 8,fontface = "bold",fontfamily="serif"),
                                         legend_height = unit(6, "cm"),
                                         grid_width = unit(0.6, "cm"))
)
ha_cn1A<-ha_cn1+ha11
draw(ha_cn1A)
pdf("regulon-sigmark.pdf",width=12,height=8)
draw(ha_cn1A)
dev.off()
