
rm(list=ls())

setwd("D:\\research\\tf_activities")

library(dorothea)

data(dorothea_hs, package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A","B","C","D","E"))

regulon <- split(regulons, regulons$tf)


tf<-unique(regulons[['tf']])

regulon<-list()
for (i in 1:length(tf)) {
  data1<-regulons[which(regulons[,1] == tf[i]),]
  mor<-data1[['mor']]
  names(mor)<-data1[['target']]
  regulon[[i]]<-list(mor)
  names(regulon[[i]])<-"tfmode"
  regulon[[i]]$confidence<-data1[['confidence']]
  names(regulon)[i]<-tf[i]
  print(i)
}

save(regulon,file = "dorothea_viper.rda")

      
      