#R4.0.4
rm(list=ls())
library(Seurat)
library(tidyverse)
library(cellcallEXT)
load("DADA2Data.RData")
countData<-countData[names(sort(rowSums(countData), decreasing = TRUE))[1:10000],]
#Let us normalize the count data
source("../cellCallEXT/counts2rpm.R")
normalizedCounts<-counts2normalized_10X(countData, "CPM", scale.factor = 10^6)
metaData<-metaData[metaData$CD48MonoLarge %in% c("CD48", "Mono"),]
normalizedCounts<-normalizedCounts[, rownames(metaData)]
rm(countData)
status<-rep("CASE", dim(metaData)[1]);status[metaData$status=="HD"]<-"CONTROL"

colnames(normalizedCounts)<-paste0(metaData$CD48MonoLarge, "_celltype_", colnames(normalizedCounts))
mt <- CreateNichConObject(data=normalizedCounts, min.feature = 3,names.field = 1,names.delim = "_",
                          source = "UMI", scale.factor = 10^6, Org = "Homo sapiens",project = "Microenvironment",status = status,
                          expDirection = "BOTH")
mt <- TransCommuProfile(object = mt,pValueCor = 0.05,CorValue = 0.05,topTargetCor=1,p.adjust = 0.1,
                        use.type="mean",method="weighted",IS_core = TRUE,Org = 'Homo sapiens')

n <- mt@data$expr_l_r_log2_scale

pathway.hyper.list <- lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
p <- plotBubble(myPub.df)
save(n, p, pathway.hyper.list, mt, file = "Debug_Result_BOTH_10000cells.RData")

