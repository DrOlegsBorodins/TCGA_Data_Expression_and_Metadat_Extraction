# TCGA_Data_Expression_and_Metadat_Extraction
TCGA_Data_Expression_and_Metadat_Extraction


# TCGA-HNSC Gene Expression and Clinical Data Integration

This R script automates the process of querying, downloading, and processing gene expression and clinical metadata from the TCGA Head and Neck Squamous Cell Carcinoma (HNSC) dataset using the `TCGAbiolinks` package. It also integrates survival data and prepares a clean, merged dataset suitable for downstream analysis, including survival analysis and machine learning.

## 📦 Requirements

Before running the script, ensure the following R packages are installed:

```r
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
install.packages(c("SummarizedExperiment", "survival", "survminer", "GEOquery", "ggplot2", "svglite", "dplyr"))
