library(Seurat)
library(dplyr)
library(tidyverse)
library(harmony)
library(patchwork)
rm(list = ls())
load('scRNA_new_celltype.Rdata')
set.seed(12345)
dir.create('Subset_CD8+T3')

scRNAsub <- subset(scRNA, idents = 'CD8+_T3')
DefaultAssay(scRNAsub) <- 'RNA'
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[, c('orig', 'nCount_RNA','nFeature_RNA','celltype',
                                             'S.Score','G2M.Score','Phase')]
scRNAsub$main.celltype <- scRNAsub$celltype
scRNAsub$celltype <- NULL


var.sub <- FindVariableFeatures(scRNAsub, assay = 'RNA') %>%
  VariableFeatures(assay = 'RNA')

length(var.sub)

scRNAsub <- NormalizeData(scRNAsub) %>% FindVariableFeatures() %>%
  ScaleData(features = rownames(scRNAsub))
scRNAsub <- RunPCA(scRNAsub, npcs = 50, verbose = FALSE)
DimPlot(scRNAsub, reduction = 'pca', group.by = 'orig')
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = 'orig', assay.use = 'RNA',
                       max.iter.harmony = 50)
