# Install TCGAbiolinks package from GitHub using Bioconductor's BiocManager and devtools
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

library(BiocManager)
library(devtools)

# Install packages
BiocManager::install("TCGAbiolinks")
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")

# Load necessary libraries for data analysis
library(TCGAbiolinks)
library(SummarizedExperiment)
library(methylGSA)
library(sesame)
library(sesameData)
library(tidyverse)
library(maftools)
library(pheatmap)

# Check versions of the packages to ensure compatibility
packageVersion("methylGSA")
packageVersion("sesame")
packageVersion("sesameData")
packageVersion("TCGAbiolinks")
packageVersion("tidyverse")
packageVersion("maftools")
packageVersion("pheatmap")
packageVersion("SummarizedExperiment")

# Retrieve a list of projects from GDC
gdcprojects <- getGDCprojects()

# View project metadata (e.g., in RStudio's Viewer)
view(gdcprojects)

# Get a summary of the TCGA-PAAD project
getProjectSummary('TCGA-PAAD')

# Build a query to retrieve gene expression data for TCGA-PAAD
query_TCGA <- GDCquery(
  project = 'TCGA-PAAD',
  data.category = 'Transcriptome Profiling',
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  access = 'open'
)

# Execute query and get result table
output_query_TCGA <- getResults(query_TCGA)

# Save results to a CSV file
output_directory <- "your_output_directory"  # e.g., "C:/xxxx/xxxx/results/TCGA_PAAD"
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
write.csv(output_query_TCGA, file.path(output_directory, "output_query_TCGA_PAAD.csv"), row.names = TRUE, quote = TRUE)

# Debug: Preview first few rows
head(output_query_TCGA)

# Download data files (to a chosen directory)
GDCdownload(query_TCGA, directory = output_directory)

# Prepare data from the downloaded files
tcga_PAAD_data <- GDCprepare(
  query = query_TCGA,
  save = FALSE,
  directory = output_directory,
  summarizedExperiment = TRUE,
  remove.files.prepared = FALSE
)

# Extract TPM expression matrix from the prepared data
PAAD_matrix <- assay(tcga_PAAD_data, 'tpm_unstrand')

# Check matrix dimensions
dim(PAAD_matrix)

# Save expression matrix to CSV
write.csv(PAAD_matrix, file.path(output_directory, "PAAD_expression_matrix.csv"), row.names = TRUE, quote = TRUE)
