#R4.0.4
rm(list=ls())
library(Seurat)
library(tidyverse)
library(cellcallEXT)
load("../ourDataSeuratAfterPCAtSNE.RData")
set.seed(888); indexCells<-sample(1:184971, 10000)
countData<-seuratObj@assays$RNA@counts[, indexCells]
metaData<-seuratObj@meta.data[indexCells,]
save(countData, metaData, file = "DADA2Data.RData")
