# Decoding the Immune Microenvironment of Clear Cell Renal Cell Carcinoma by Single-cell Profiling to Aid Immunotherapy

## Introduction
Data and downstream analysis for 2021 submission "Decoding the Immune Microenvironment of Clear Cell Renal Cell Carcinoma by Single-cell Profiling to Aid Immunotherapy" Yang Wang et. al
There are () main steps involved in data processing which organized in different folders.

## Download and read the datasets
```
cd "1-ReadData"
bash download_dataset.sh
Rscript readData_head_neck.R
Rscript readData_melanoma.R
cd ../
```
All samples were obtained from the National Center for Biotechnology Information GEO dataset.Raw data was converted into a Seurat object by R package Seurat (v 3.1.2)

## QC
```
cd "2-QC"
Rscript QC.R 
cd ../
```

## Normalization of upstream analysis results
```
cd "3-Normalization"
Rscript normalization.R 
cd ../
```
After quality control, the Seurat object was normalized by the function SCTransform of Seurat package.

## Unsupervised clustering of normalized results
```
cd "4-Clustering"
Rscript metabolic_landscape.R 
cd ../
```
We performed principal component analysis (PCA) on an integrated data matrix.And then they were visualized with 2D UMAP plots.

## Differential expressed genes
```
cd 5-DE
Rscript TCGA_pathway_activity.R
cd ../
```
The differentially expressed genes (DEGs) between pathological tissues (obtained by TCGA website) and normal tissues were initially screened by R software Limma package and EdgeR package.

The bulk RNA-seq data was downloaded from TCGA website.

## Metabolic pathway heterogeneity
cd 6-PathwayHeterogeneity
Rscript intra_malignant_heterogeneity.R melanoma
Rscript intra_malignant_heterogeneity.R head_neck
Rscript intra_non-malignant_heterogeneity.R melanoma
Rscript intra_non-malignant_heterogeneity.R head_neck
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R melanoma
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R head_neck
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot-CCLE.R
Rscirpt GeneSignature-of-Low_OXPHOS_Glycolysis_Hypoxia.R melanoma
Rscript GeneSignature-of-Low_OXPHOS_Glycolysis_Hypoxia.R head_neck
cd ..
In this step, the PCA and GSEA analysis will be performed to investigate the metabolic pathway heterogeneity across single cells in malignant and non-malignant cell populations. The scatter plots will be performed to compare activities of OXPHOS, glycolysis and response to hypoxia in single malignant cells and cultured cell lines from CCLE database. The gene signatures in single cells with low OXPHOS/glycolysis/hypoxia activity will be identified and stored as the text files, which can be used as the input of GO analysis on the website: http://metascape.org

Metabolic features of nonmalignant cell subtypes
cd 7-MetabolicFeatures
Rscript non-malignant_subtype.R melanoma
Rscript non-malignant_subtype.R head_neck
cd ..
The metabolic features in different T-cell subtypes and fibroblast subtypes will be identified in this step.

Contact
