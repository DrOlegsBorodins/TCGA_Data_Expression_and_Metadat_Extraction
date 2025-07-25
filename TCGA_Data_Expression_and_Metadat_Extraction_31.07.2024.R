# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# Load libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)  # For Kaplan-Meier plots
library(GEOquery)
library(ggplot2)
library(svglite)
library(dplyr)

# Define working directories (change as needed)
data_directory <- "data_directory"         # e.g., "C:/xxxx/xxxx/Desktop/data/TCGA_HNSC"
output_directory <- "output_directory"     # e.g., "C:/xxxx/xxxx/results/"



# Query TCGA HNSC data
query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download and prepare the data
GDCdownload(query, directory = data_directory)
data <- GDCprepare(query, directory = data_directory)

# Extract expression data
exprSet <- as.data.frame(SummarizedExperiment::assay(data))

# Extract clinical data (metadata)
clinicalData <- colData(data)

# Process survival data
survivalTime <- as.numeric(clinicalData$days_to_death)
survivalEvent <- ifelse(clinicalData$vital_status == "Dead", 1, 0)

valid_indices <- !is.na(survivalTime) & !is.na(survivalEvent)
survivalTime <- survivalTime[valid_indices]
survivalEvent <- survivalEvent[valid_indices]
exprSet <- exprSet[, valid_indices]

survival_df <- data.frame(
  SampleID = colnames(exprSet),
  SurvivalTime = survivalTime,
  SurvivalEvent = survivalEvent
)

# List of sample IDs to exclude --> codeing for TCGA https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
columns_to_remove <- c(
  "TCGA-CV-7091-11A-01R-2016-07", "TCGA-CV-6961-11A-01R-1915-07",
  "TCGA-CV-7242-11A-01R-2016-07", "TCGA-HD-A6HZ-11A-11R-A31N-07",
  "TCGA-CV-7425-11A-01R-2081-07", "TCGA-WA-A7GZ-11A-11R-A34R-07",
  "TCGA-CV-7238-11A-01R-2016-07", "TCGA-CV-6943-11A-01R-1915-07",
  "TCGA-CV-6962-11A-01R-1915-07", "TCGA-CV-7177-11A-01R-2016-07",
  "TCGA-CV-7183-11A-01R-2016-07", "TCGA-CV-7432-11A-01R-2132-07",
  "TCGA-CV-6959-11A-01R-1915-07", "TCGA-CV-7416-11A-01R-2081-07",
  "TCGA-CV-7252-11A-01R-2016-07", "TCGA-CV-7235-11A-01R-2016-07",
  "TCGA-CV-6934-11A-01R-1915-07", "TCGA-CV-7434-11A-01R-2132-07",
  "TCGA-CV-7178-11A-01R-2016-07", "TCGA-CV-6956-11A-01R-2016-07",
  "TCGA-CV-7097-11A-01R-2016-07", "TCGA-CV-7101-11A-01R-2016-07",
  "TCGA-CV-7103-11A-01R-2016-07", "TCGA-CV-6935-11A-01R-1915-07",
  "TCGA-H7-A6C5-11A-11R-A30B-07", "TCGA-CV-7250-11A-01R-2016-07",
  "TCGA-CV-6936-11A-01R-1915-07", "TCGA-CV-7245-11A-01R-2016-07",
  "TCGA-CV-6955-11A-01R-2016-07", "TCGA-CV-7406-11A-01R-2081-07",
  "TCGA-HD-8635-11A-01R-2403-07", "TCGA-CV-7255-11A-01R-2016-07",
  "TCGA-CV-6938-11A-01R-1915-07", "TCGA-HD-A6I0-11A-11R-A31N-07",
  "TCGA-CV-7437-11A-01R-2132-07", "TCGA-CV-7440-11A-01R-2187-07",
  "TCGA-CV-7438-11A-01R-2132-07", "TCGA-CV-6933-11A-01R-1915-07",
  "TCGA-CV-6960-11A-01R-2016-07", "TCGA-CV-6939-11A-01R-1915-07",
  "TCGA-H7-A6C4-11A-21R-A466-07", "TCGA-CV-7423-11A-01R-2081-07",
  "TCGA-CV-7424-11A-01R-2081-07", "TCGA-CV-7261-11A-01R-2016-07"
)

# Remove unwanted samples
columns_indices <- which(colnames(exprSet) %in% columns_to_remove)
exprSet2 <- exprSet[, -columns_indices]
survival_df_2 <- survival_df[!survival_df$SampleID %in% columns_to_remove, ]

# List of genes of interest
genes <- c('ENSG00000135424.18')  # Replace with your gene ID(s)

# Align expression and survival data
aligned_exprSet <- exprSet2[, match(survival_df_2$SampleID, colnames(exprSet2))]
gene_expr <- exprSet2[rownames(exprSet2) %in% genes, ]

if (!all(rownames(gene_expr) %in% genes)) {
  stop("Some gene IDs in 'genes' are not present in 'exprSet2'")
}

gene_of_interest_expr <- as.data.frame(t(gene_expr))

# Merge expression with survival data
combined_df <- merge(survival_df_2, gene_of_interest_expr, by.x = "SampleID", by.y = "row.names")
colnames(combined_df)[colnames(combined_df) == "SampleID"] <- "PatientID"
colnames(combined_df)[ncol(combined_df)] <- "GeneExpression"
combined_df <- combined_df[, c("PatientID", "GeneExpression", "SurvivalTime", "SurvivalEvent")]

# Extract and clean clinical features
selected_clinical_features <- clinicalData %>%
  as.data.frame() %>%
  select(
    PatientID = barcode,
    patient,
    sample,
    shortLetterCode,
    definition,
    sample_submitter_id,
    sample_type_id,
    tumor_descriptor,
    sample_id,
    sample_type,
    composition,
    days_to_collection,
    state,
    initial_weight,
    preservation_method,
    intermediate_dimension,
    pathology_report_uuid,
    submitter_id,
    shortest_dimension,
    oct_embedded,
    specimen_type,
    longest_dimension,
    is_ffpe,
    tissue_type,
    synchronous_malignancy,
    ajcc_pathologic_stage,
    days_to_diagnosis,
    treatments,
    last_known_disease_status,
    tissue_or_organ_of_origin,
    days_to_last_follow_up,
    age_at_diagnosis,
    primary_diagnosis,
    ajcc_clinical_stage,
    prior_malignancy,
    year_of_diagnosis,
    prior_treatment,
    ajcc_staging_system_edition,
    ajcc_pathologic_t,
    morphology,
    ajcc_clinical_m,
    ajcc_pathologic_n,
    ajcc_pathologic_m,
    ajcc_clinical_n,
    ajcc_clinical_t,
    classification_of_tumor,
    diagnosis_id,
    icd_10_code,
    site_of_resection_or_biopsy,
    tumor_grade,
    progression_or_recurrence,
    cigarettes_per_day,
    alcohol_history,
    exposure_id,
    years_smoked,
    pack_years_smoked,
    race,
    gender,
    ethnicity,
    vital_status,
    age_at_index,
    days_to_birth,
    year_of_birth,
    demographic_id,
    year_of_death,
    days_to_death,
    bcr_patient_barcode,
    primary_site,
    project_id,
    disease_type,
    name,
    releasable,
    released
  )

# Print missing values per column
missing_values <- sapply(selected_clinical_features, function(x) sum(is.na(x)))
print("Missing values in clinical data:")
print(missing_values)

# Replace NAs with "Unknown"
selected_clinical_features[is.na(selected_clinical_features)] <- "Unknown"

# Merge full clinical info with gene and survival data
combined_df <- combined_df %>%
  left_join(selected_clinical_features, by = "PatientID")

# Flatten list-columns (if any) and save
dir.create(output_directory, showWarnings = FALSE)

combined_df <- combined_df %>%
  mutate(across(where(~is.list(.)), ~sapply(., function(x) paste(unlist(x), collapse = ";"))))

write.csv(combined_df, file = file.path(output_directory, "patient_gene_clinical_data.csv"), row.names = FALSE)
