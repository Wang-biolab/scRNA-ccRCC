# load packages
library(Seurat)
library(tidyverse)
rm(list=ls())


dir <- dir("RawData/")
dir <- paste0("RawData/", dir)

# Name the samples in file order（do not start the name with a number, and have no spaces in between）
group = c('p1','p2','p3','t1','t2','t3','rb1','rb2','rb3','rn1','rn2','rn3','rt1','rt2','rt3')

# Create seurat objects in batches using loop command
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNA <- CreateSeuratObject(counts, project= group[i])
  save(scRNA, file= paste('RawData/scRNA_',group[i],'.Rdata',sep=''))}

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i],gene.column = 1)
  counts <- as.data.frame(counts)
  mat <- counts[rowSums(counts)!=0,]
  rm(counts)

  scRNAlist[[i]] <- CreateSeuratObject(mat,
                                       project=group[i],
                                       min.cells=3,
                                       min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = group[i])
}


# Name the list and save the data
names(scRNAlist) <- group
scRNA <- merge(scRNAlist[[1]],
               y= scRNAlist[2:length(scRNAlist)])