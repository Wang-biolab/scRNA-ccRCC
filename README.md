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

## Landscape of the metabolic gene expression profile
cd "4-Clustering"
Rscript metabolic_landscape.R melanoma
Rscript metabolic_landscape.R head_neck
Rscript inter_tumor_distance.R melanoma
Rscript inter_tumor_distance.R head_neck
cd ../
The t-SNE algorithm will be performed in this step for visualizing metabolic gene expression in millions of cells (The result may be slightly different with the figure in manuscript due to the random initialization). The spearman correlation matrix will aslo be generated to show the inter-tumor heterogeneity using the metabolic genes.

Metabolic pathway activities in different cell types
cd 5-PathwayActivity
Rscript scRNA_pathway_activity.R melanoma
Rscript scRNA_pathway_activity.R head_neck
Rscript TCGA_pathway_activity.R
cd ..
This step will calculate the metabolic pathway activities for different single cell populations or bulk tumor/normal samples. The scatter plot will show the discrepancy of pathway activities between single malignant cells and bulk tumors. The violin plot will show the distribution of metabolic pathway activities in single cell populations or bulk tumor/normal samples.

The bulk RNA-seq data was downloaded from TCGA website, please see the instruction of data downloading and preprocessing in Data/TCGA/README.md

Metabolic pathway heterogeneity
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
