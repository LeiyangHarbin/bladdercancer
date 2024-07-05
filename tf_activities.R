rm(list=ls())

setwd("D:\\research\\tf_activities")

library(MethReg)
library(GSVA)
library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
library(GEOquery)

#####################################

data(dorothea_hs, package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A","B","C","D","E"))

save(regulons,file = "regulons.rda")

load("D:/research/data_new/Bladder_7.rda")


result<-list()
for (i in 1:length(Bladder_7)) {
  dset<-exprs(Bladder_7[[i]])
  tf_activities<-run_viper(dset,regulons,
                            options = list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))
  result[[i]]<-tf_activities
  names(result)[i]<-names(Bladder_7)[i]
  print(i)
}

save(result,file = "tf_activities_7.rda")
