# Decoding the Immune Microenvironment of Clear Cell Renal Cell Carcinoma by Single-cell Profiling to Aid Immunotherapy

## Introduction
Data and downstream analysis for 2021 submission "Decoding the Immune Microenvironment of Clear Cell Renal Cell Carcinoma by Single-cell Profiling to Aid Immunotherapy" Yang Wang et. al
There are 6 main steps involved in data processing which organized in different folders.

## Download and read the datasets
```
cd "1-ReadData"
Rscript readDATA.R
cd ../
```
All samples were obtained from the National Center for Biotechnology Information GEO dataset.Raw data was converted into a Seurat object by R package Seurat (v 3.1.2)

## Normalization of upstream analysis results
```
cd "2-Normalization"
Rscript normalization.R 
cd ../
```
After quality control, the Seurat object was normalized by the function SCTransform of Seurat package.

## Unsupervised clustering of normalized results
```
cd "3-Clustering"
Rscript clustering.R 
cd ../
```
We performed principal component analysis (PCA) on an integrated data matrix.And then they were visualized with 2D UMAP plots.


## Pathway and Functional annotation analysis
```
cd "4-PathwayFunctionalAnnotation"
Rscript pathway.R 
cd ../
```
A gene functional enrichment analysis was performed based on the marker genes in each cell cluster using KEGG pathway enrichment analysis. These differential expressed genes were loaded into clusterProfiler for the GO and KEGG pathway enrichment analysis

## Single-cell trajectory analysis
```
cd 5-TrajectoryAnalysis
Rscript monocle.R 
cd ../
```
We used the 'monocle 2' for cell trajectory analysis.Based on the ‘DDRTree’ method,the data were reduced to two dimensional, and then, the cells were ordered along the trajectory.

## Quantify tumor and normal tissue interactions
```
cd 6-CellChat
Rscript chellchat.R 
cd ../
```
We used CellChat R package to characterize cell-cell communication across all cell types in ccRCC.

## Contact
yangwang@hubu.edu.cn
