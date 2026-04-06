#===============================================================================
# Script: 03_WGCNA_coexpression_analysis.R
# Description: Weighted Gene Co-expression Network Analysis (WGCNA) to identify 
#              clusters of co-expressed genes associated with herbicide resistance
#              in Amaranthus palmeri. WGCNA reveals gene network connections that
#              cannot be identified by differential expression analysis alone.
#
# Input files:
#   - raw_counts_star.csv: Raw gene count matrix (genes x samples)
#   - metadata1+2+3.csv: Sample metadata with batch and treatment information
#
# Output:
#   - Module eigengenes and gene assignments
#   - Module-trait correlation heatmaps
#   - Gene lists for modules of interest
#
# WGCNA Parameters Used:
#   - Network type: signed
#   - Soft-thresholding power: 14
#   - Maximum block size: 8,000 genes
#   - Module merge cut height: 0.25
#   - Minimum gene count filter: 10 counts in 75% of samples
#
# Tutorial references: 
#   - https://www.youtube.com/watch?v=PvBf65Y8Cqk
#   - https://www.youtube.com/watch?v=gYE59uEMXT4
#===============================================================================

#-------------------------------------------------------------------------------
# Setup: Set working directory
#-------------------------------------------------------------------------------
setwd("/path/to/your/working/directory")

#-------------------------------------------------------------------------------
# Install required packages (run once)
#-------------------------------------------------------------------------------
# install.packages("BiocManager")
# BiocManager::install("WGCNA", force = TRUE)
# BiocManager::install("DESeq2")
# BiocManager::install("CorLevelPlot")
# install.packages("tidyverse")
# install.packages("patchwork")

#-------------------------------------------------------------------------------
# Load required packages
#-------------------------------------------------------------------------------
library(WGCNA, quietly = TRUE)   # WGCNA analysis
library(tidyverse)                # Data manipulation
library(DESeq2)                   # Normalization
library(CorLevelPlot)             # Module-trait correlation visualization
library(gridExtra)                # Grid plotting
library(dplyr)                    # Data manipulation
library(ggplot2)                  # Plotting
library(readxl)                   # Read Excel files
library(patchwork)                # Combine plots
library(writexl)                  # Write Excel files

# Enable multi-threading for WGCNA
allowWGCNAThreads()

#===============================================================================
# STEP 1: DATA LOADING AND PREPARATION
#===============================================================================

#-------------------------------------------------------------------------------
# Load count matrix
#-------------------------------------------------------------------------------
counts_df <- read_excel("batch_all_counts.xlsx")

# Convert to matrix with gene IDs as rownames
counts_mat <- as.matrix(counts_df[, -1])
rownames(counts_mat) <- counts_df[[1]]
storage.mode(counts_mat) <- "integer"

pa_counts <- counts_mat

# Load metadata
meta <- read.csv("metadata1+2+3.csv", row.names = 1, stringsAsFactors = FALSE)

#===============================================================================
# STEP 2: QUALITY CONTROL - OUTLIER DETECTION
#===============================================================================

#-------------------------------------------------------------------------------
# Detect outlier genes using goodSamplesGenes function
#-------------------------------------------------------------------------------
# Note: WGCNA requires samples as rows and genes as columns (transposed)
gsg <- goodSamplesGenes(t(pa_counts))
summary(gsg)

# Check if all samples pass quality control
gsg$allOK    # Should be TRUE

# View summary of good/bad genes and samples
table(gsg$goodGenes)
table(gsg$goodSamples)

# Remove outlier genes
pa_counts <- pa_counts[gsg$goodGenes == TRUE, ]
cat("Genes remaining after QC:", nrow(pa_counts), "\n")

#-------------------------------------------------------------------------------
# Detect outlier samples using hierarchical clustering
#-------------------------------------------------------------------------------
# Cluster samples to identify outliers
htree <- hclust(dist(t(pa_counts)), method = "average")
plot(htree, main = "Sample Clustering to Detect Outliers")

# Define samples to exclude based on clustering
# (Modify this list based on your clustering results)
out.samples <- c('KCTR8NT1', 'KCTR7NT1')

# Remove outlier samples
pa_counts <- pa_counts[, !(colnames(pa_counts) %in% out.samples)]
cat("Samples remaining after outlier removal:", ncol(pa_counts), "\n")

# Re-plot clustering after outlier removal
htree <- hclust(dist(t(pa_counts)), method = "average")
plot(htree, main = "Sample Clustering After Outlier Removal")

#===============================================================================
# STEP 3: NORMALIZATION
#===============================================================================

#-------------------------------------------------------------------------------
# Update metadata to exclude outlier samples
#-------------------------------------------------------------------------------
meta <- meta %>%
  filter(!row.names(.) %in% out.samples)

# Verify sample correspondence
all(colnames(pa_counts) %in% rownames(meta))
all(rownames(meta) == colnames(pa_counts))

# Convert variables to factors
meta$group <- factor(meta$group)
meta$batch <- factor(meta$batch)

#-------------------------------------------------------------------------------
# Create DESeq2 object for normalization
#-------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = pa_counts,
  colData = meta,
  design = ~ group + batch    # Include batch as covariate
)

cat("Initial gene count:", nrow(dds), "\n")

#-------------------------------------------------------------------------------
# Filter low-count genes
# Remove genes with < 10 counts in at least 75% of samples
#-------------------------------------------------------------------------------
min_samples <- round(ncol(pa_counts) * 0.75)
dds <- dds[rowSums(counts(dds) >= 10) >= min_samples, ]
cat("Genes after filtering:", nrow(dds), "\n")

#-------------------------------------------------------------------------------
# Variance stabilizing transformation
#-------------------------------------------------------------------------------
dds.norm <- vst(dds)

# Get normalized counts (transposed: samples as rows, genes as columns)
norm.counts <- assay(dds.norm) %>% t()

# Verify dimensions
dim(norm.counts)

#===============================================================================
# STEP 4: NETWORK CONSTRUCTION
#===============================================================================

#-------------------------------------------------------------------------------
# Choose soft-thresholding power
#-------------------------------------------------------------------------------
# Test a range of powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Calculate network topology for each power
sft <- pickSoftThreshold(
  norm.counts,
  powerVector = power,
  networkType = "signed",
  verbose = 5
)

# Extract fit indices
sft.data <- sft$fitIndices
head(sft.data)

# Visualize scale-free topology fit
P1 <- ggplot(sft.data, aes(x = Power, y = SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.02) +
  geom_hline(yintercept = 0.8, color = 'red', linetype = "dashed") +
  labs(x = 'Soft Threshold (Power)', y = 'Scale Free Topology Model Fit (R²)') +
  theme_classic() +
  ggtitle("Scale-Free Fit Index")

P2 <- ggplot(sft.data, aes(x = Power, y = mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 50) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') +
  theme_classic() +
  ggtitle("Mean Connectivity")

# Combine plots
P1 / P2

# Save plot
ggsave("soft_threshold_selection.pdf", width = 8, height = 10)
#we will use power with high r sq but not excessively large (14)

#-------------------------------------------------------------------------------
# Build co-expression network
#-------------------------------------------------------------------------------
# Selected soft power (based on R² > 0.8 with reasonable connectivity)
soft.power <- 14

# Temporarily replace cor function with WGCNA's cor
temp.cor <- cor
cor <- WGCNA::cor

# Convert to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# Build network using blockwise modules
bw.nwt <- blockwiseModules(
  norm.counts,
  maxBlockSize = 8000,           # Maximum genes per block
  corType = "pearson",           # Correlation type
  TOMType = "signed",            # Topological overlap type
  power = soft.power,            # Soft-thresholding power
  mergeCutHeight = 0.25,         # Module merging threshold
  numericLabels = FALSE,         # Use color labels
  randomSeed = 1234,             # For reproducibility
  verbose = 3
)

# Restore original cor function
cor <- temp.cor

#===============================================================================
# STEP 5: MODULE ANALYSIS
#===============================================================================

#-------------------------------------------------------------------------------
# Extract module eigengenes
#-------------------------------------------------------------------------------
module.eigenes <- bw.nwt$MEs
head(module.eigenes)

# Count genes per module
table(bw.nwt$colors)
# Note: Grey module contains unassigned genes

#-------------------------------------------------------------------------------
# Visualize module dendrogram
#-------------------------------------------------------------------------------
# Plot for first block
b <- 1
g <- if (is.list(bw.nwt$blockGenes)) bw.nwt$blockGenes[[b]] else which(bw.nwt$blockGenes == b)

plotDendroAndColors(
  bw.nwt$dendrograms[[b]],
  cbind(bw.nwt$unmergedColors[g], bw.nwt$colors[g]),
  c("Unmerged", "Merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors"
)

#===============================================================================
# STEP 6: MODULE-TRAIT CORRELATION
#===============================================================================

#-------------------------------------------------------------------------------
# Create binary trait variable for herbicide resistance
#-------------------------------------------------------------------------------
traits <- meta %>%
  mutate(resistant.binary = ifelse(grepl('KCTR', group), 1, 0))

#-------------------------------------------------------------------------------
# Calculate module-trait correlations
#-------------------------------------------------------------------------------
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Correlate module eigengenes with traits
module.trait.corr <- cor(module.eigenes, trait.all[, 3:14], use = 'p')
module.trait.corr.pval <- corPvalueStudent(module.trait.corr, nSamples)

# Preview p-values
head(module.trait.corr.pval)

#-------------------------------------------------------------------------------
# Visualize module-trait correlations
#-------------------------------------------------------------------------------
heatmap.data <- cbind(module.eigenes, trait.all[, 3:14])
heatmap.data <- as.data.frame(heatmap.data)

# Create correlation level plot
CorLevelPlot(
  heatmap.data,
  x = names(heatmap.data)[24:29],    # Trait columns
  y = names(heatmap.data)[1:23],     # Module columns
  col = c("#3B4CC0", "#8DB0FE", "#F7F7F7", "#F4A582", "#B40426"),
  main = "Module-Trait Correlations"
)

#===============================================================================
# STEP 7: EXTRACT GENES FROM MODULES OF INTEREST
#===============================================================================

#-------------------------------------------------------------------------------
# Get gene-module mapping
#-------------------------------------------------------------------------------
module.gene.mapping <- as.data.frame(bw.nwt$colors)
colnames(module.gene.mapping) <- "module"

# Extract genes from specific modules (modify module names as needed)
# Example: Extract genes from "royalblue" module
mod.royalblue.genes <- module.gene.mapping %>%
  filter(module == 'royalblue') %>%
  rownames()
mod.royalblue.genes <- as.data.frame(mod.royalblue.genes)
colnames(mod.royalblue.genes) <- "gene_id"

# Example: Extract genes from "lightyellow" module
mod.lightyellow.genes <- module.gene.mapping %>%
  filter(module == 'lightyellow') %>%
  rownames()
mod.lightyellow.genes <- as.data.frame(mod.lightyellow.genes)
colnames(mod.lightyellow.genes) <- "gene_id"

# Export gene lists
write_xlsx(mod.royalblue.genes, path = "results/mod_royalblue_genes.xlsx")
write_xlsx(mod.lightyellow.genes, path = "results/mod_lightyellow_genes.xlsx")
