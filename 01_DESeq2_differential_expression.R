#===============================================================================
# Script: 01_DESeq2_differential_expression.R
# Description: Differential gene expression analysis of RNA-Seq data from 
#              herbicide-resistant (KCTR/KCTR-G2) vs susceptible (KSS) 
#              Amaranthus palmeri populations using DESeq2.
#
# Input files:
#   - raw_counts_star.csv: Gene count matrix (genes x samples)
#   - separate sample list.csv: Sample metadata for batch-specific analysis
#
# Output: CSV files with differential expression results for each comparison
#
# Comparisons performed:
# Batch 1
#   - Non-treated (NT) R vs S
#   - 2,4-D treated R vs S
#   - Mesotrione treated R vs S
# Batch 2
#   - Non-treated (NT) R vs S
#   - Atrazine treated R vs S
#   - Chlorsulfuron treated R vs S
# Batch 3
#   - Non-treated (NT) R vs S
#   - Lactofen treated R vs S
#
#===============================================================================

#-------------------------------------------------------------------------------
# Setup: Set working directory (modify path for your system)
#-------------------------------------------------------------------------------
setwd("/path/to/your/working/directory")

#-------------------------------------------------------------------------------
# Install and load required packages
#-------------------------------------------------------------------------------
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DESeq2 (only needed once)
# BiocManager::install("DESeq2")

# Load required libraries
library(DESeq2)      # For differential expression analysis
library(tidyverse)
library(readxl)# For data manipulation

#===============================================================================
# STEP 1: Load and prepare count data
#===============================================================================

# Read gene count matrix
# Format: rows = genes, columns = samples
pa_counts <- as.matrix(
  read.csv("raw-counts_star.csv",
           header = TRUE,
           row.names = 1)
)

# Pick up the desired treatments from to analyze for Batch 1, 2 or 3
# This script will follow as an example for batch 1

# Preview count matrix
head(pa_counts)

#-------------------------------------------------------------------------------
# Load sample metadata
#-------------------------------------------------------------------------------
# Option 1: Separate batch analysis
# Pick up the metadata for the desired treatments from to analyze for Batch 1, 2 or 3
df.sample <- read.csv(
  #"metadata for batch 1,2 or 3.csv"
  )
rownames(df.sample) <- df.sample$Sample

#-------------------------------------------------------------------------------
# Verify sample-count matrix correspondence
#-------------------------------------------------------------------------------
# Check if all sample names match between count matrix and metadata
all(colnames(pa_counts) %in% rownames(df.sample))
all(colnames(pa_counts) == rownames(df.sample))

#===============================================================================
# STEP 2: Create DESeq2 dataset objects
#===============================================================================

# Create DESeqDataSet for batch-specific analysis
# Design formula: ~ group (where group = population_treatment combination)
dds_separate <- DESeqDataSetFromMatrix(countData = pa_counts,
                                        colData = df.sample,
                                        design = ~ group)

#===============================================================================
# STEP 3: Run DESeq2 differential expression pipeline
#===============================================================================
# DESeq() performs: estimation of size factors, estimation of dispersion,
# and negative binomial GLM fitting

dds_separate <- DESeq(dds_separate)

#===============================================================================
# STEP 4: Extract pairwise comparisons
#===============================================================================
#-------------------------------------------------------------------------------
# Helper function for extracting and saving DE results
#-------------------------------------------------------------------------------
extract_and_save_results <- function(dds, contrast_vector, output_name) {
  # Extract results for specified contrast
  res <- results(dds, contrast = contrast_vector)
  
  # Sort by adjusted p-value
  res_ordered <- res[order(res$padj), ]
  
  # Export to CSV
  write.csv(res_ordered, paste0(output_name, ".csv"))
  
  # Return ordered results
  return(res_ordered)
}

#-------------------------------------------------------------------------------
# NON-TREATED COMPARISONS (NT)
#-------------------------------------------------------------------------------

# KCTR NT vs KSS NT (combined batches) - Primary comparison
KCTR_NT_vs_KSS_NT <- extract_and_save_results(
  dds_separate, 
  c("group", "KCTR_NT", "KSS_NT"),
  "results/KCTR_NT_vs_KSS_NT_combined"
)

#-------------------------------------------------------------------------------
# 2,4-D TREATMENT COMPARISONS
#-------------------------------------------------------------------------------

# Resistant vs Susceptible after 2,4-D treatment
KCTR_D_vs_KSS_D <- extract_and_save_results(
  dds_separate,
  c("group", "KCTR_24D", "KSS_24D"),
  "results/KCTR_24D_vs_KSS_24D_combined"
)

# Treatment effect within populations
KCTR_D_vs_KCTR_NT <- extract_and_save_results(
  dds_separate,
  c("group", "KCTR_24D", "KCTR_NT"),
  "results/KCTR_24D_vs_KCTR_NT"
)

KSS_D_vs_KSS_NT <- extract_and_save_results(
  dds_separate,
  c("group", "KSS_24D", "KSS_NT"),
  "results/KSS_24D_vs_KSS_NT"
)

#-------------------------------------------------------------------------------
# MESOTRIONE TREATMENT COMPARISONS
#-------------------------------------------------------------------------------

# Resistant vs Susceptible after mesotrione treatment
KCTR_Meso_vs_KSS_Meso <- extract_and_save_results(
  dds_separate,
  c("group", "KCTR_MESO", "KSS_MESO"),
  "results/KCTR_MESO_vs_KSS_MESO"
)

# Treatment effect within populations
KCTR_Meso_vs_KCTR_NT <- extract_and_save_results(
  dds_separate,
  c("group", "KCTR_MESO", "KCTR_NT"),
  "results/KCTR_MESO_vs_KCTR_NT"
)

KSS_Meso_vs_KSS_NT <- extract_and_save_results(
  dds_separate,
  c("group", "KSS_MESO", "KSS_NT"),
  "results/KSS_MESO_vs_KSS_NT"
)

#===============================================================================
# STEP 5: Visualization (MA Plot)
#===============================================================================

# Example: MA plot for chlorsulfuron comparison
res_viz <- results(dds_separate, 
                   contrast = c("group", "KCTR_CHL", "KSS_CHL"), 
                   alpha = 0.05)
plotMA(res_viz, main = "KCTR vs KSS - Chlorsulfuron Treatment")

#===============================================================================
# Session Info
#===============================================================================
sessionInfo()

#-------------------------------------------------------------------------------
# Candidate gene selection criteria (for downstream filtering):
# 1. Adjusted p-value < 0.05
# 2. Absolute log2 fold change >= 2
# 3. Functional annotation associated with herbicide metabolism
# 4. Upregulation in resistant (R) vs susceptible (S) population
#-------------------------------------------------------------------------------
