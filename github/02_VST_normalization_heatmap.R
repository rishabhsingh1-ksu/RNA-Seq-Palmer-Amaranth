#===============================================================================
# Script: 02_VST_normalization_heatmap.R
# Description: Variance stabilizing transformation (VST) of RNA-Seq count data
#              and visualization of gene expression patterns via heatmaps.
#
# Input files:
#   - raw_counts_star.csv: Raw gene count matrix (genes x samples)
#   - metadata1+2+3.csv: Sample metadata with batch and treatment information
#
# Output:
#   - gene_summary1+2+3.xlsx: VST-normalized gene expression summary
#   - Heatmaps for top variable genes and candidate genes
#
# Methods:
#   - Size factor normalization (DESeq2)
#   - Variance stabilizing transformation (VST)
#   - Replicate summarization (mean, SD)
#   - Hierarchical clustering heatmaps
#===============================================================================

#-------------------------------------------------------------------------------
# Setup: Set working directory
#-------------------------------------------------------------------------------
setwd("/path/to/your/working/directory")

#-------------------------------------------------------------------------------
# Load required packages
#-------------------------------------------------------------------------------
library(DESeq2)       # For normalization and VST
library(tidyverse)    # For data manipulation
library(readxl)       # For reading Excel files
library(writexl)      # For writing Excel files
library(pheatmap)     # For heatmap visualization
library(ggplot2)      # For custom plotting
library(matrixStats)  # For row variance calculations

#===============================================================================
# PART 1: VST NORMALIZATION
#===============================================================================

#-------------------------------------------------------------------------------
# Step 1: Load raw count matrix
#-------------------------------------------------------------------------------
counts_df <- as.matrix(read.csv("batch_all_counts.xlsx", header =TRUE))

# Extract gene IDs and convert to matrix
gene_col <- colnames(counts_df)[1]    # First column contains gene IDs
counts_mat <- as.data.frame(counts_df)
rownames(counts_mat) <- counts_mat[[gene_col]]
counts_mat[[gene_col]] <- NULL

# Convert to integer matrix (required for DESeq2)
counts_mat <- as.matrix(counts_mat)
storage.mode(counts_mat) <- "integer"

# Preview count matrix
dim(counts_mat)
counts_mat[1:5, 1:5]

#-------------------------------------------------------------------------------
# Step 2: Load sample metadata
#-------------------------------------------------------------------------------
meta <- read.csv("metadata1+2+3.csv", stringsAsFactors = FALSE)

# Set sample names as rownames
rownames(meta) <- meta$sample

# Convert design variables to factors
meta$group <- factor(meta$group)    # Treatment groups
meta$batch <- factor(meta$batch)    # Batch information

# Preview metadata
head(meta)

#-------------------------------------------------------------------------------
# Step 3: Create DESeq2 object and estimate size factors
#-------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = meta,
  design = ~ group + batch    # Account for batch effects
)

# Estimate size factors (library size normalization)
dds <- estimateSizeFactors(dds)

# View size factors
sizeFactors(dds)

#-------------------------------------------------------------------------------
# Step 4: Variance stabilizing transformation
#-------------------------------------------------------------------------------
# VST stabilizes variance across the range of mean values
# blind = FALSE uses design information to avoid over-correction
vsd <- vst(dds, blind = FALSE)

# Extract VST-transformed values (log2-like scale)
vst_mat <- assay(vsd)

# Preview VST matrix
dim(vst_mat)
vst_mat[1:5, 1:5]

#-------------------------------------------------------------------------------
# Step 5: Summarize expression across replicates
#-------------------------------------------------------------------------------
# Convert to long format and join with metadata
vst_long <- as.data.frame(vst_mat) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr_vst") %>%
  left_join(meta, by = "sample")

# Calculate mean and SD per gene per treatment group
gene_summary <- vst_long %>%
  group_by(gene, group) %>%
  summarise(
    mean_expr = mean(expr_vst, na.rm = TRUE),
    sd_expr = sd(expr_vst, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Preview gene summary
head(gene_summary)

# Export gene summary
write_xlsx(gene_summary, path = "vst_norm_counts.xlsx")

#===============================================================================
# PART 2: HEATMAP VISUALIZATION
#===============================================================================

#-------------------------------------------------------------------------------
# Heatmap 1: Top variable genes across all samples
#-------------------------------------------------------------------------------
# Select top N most variable genes
top_n <- 50
gene_vars <- rowVars(vst_mat)
top_genes <- order(gene_vars, decreasing = TRUE)[1:top_n]
hm_mat <- vst_mat[top_genes, , drop = FALSE]

# Create annotation data frame for columns
ann <- meta %>%
  select(group) %>%
  as.data.frame()
rownames(ann) <- rownames(meta)

# Generate heatmap
pheatmap(
  hm_mat,
  scale = "row",                    # Z-score normalization per gene
  annotation_col = ann,              # Add treatment group annotation
  show_colnames = TRUE,
  show_rownames = FALSE,             # Hide gene names for readability
  main = "Top 50 Variable Genes (VST-normalized)"
)

#-------------------------------------------------------------------------------
# Heatmap 2: Mean expression of candidate genes across treatments
#-------------------------------------------------------------------------------
# Read pre-computed gene summaries with candidate genes
# This requires a gene_summary file with columns: gene_id, group, anno, mean.exp
data <- read_excel("vst_norm_counts", sheet = 2)

# Convert relevant columns to factors
data <- data %>%
  mutate(across(c(gene_id, group, anno), as.factor))

# Remove missing values
data <- na.omit(data)

# Sort annotations alphabetically
data$anno <- factor(data$anno, levels = sort(unique(data$anno)))

# Define color palette (white -> salmon -> red)
col_palette <- c("#F7F7F7", "#F4A582", "#B40426")

# Create candidate gene heatmap
ggplot(data, aes(x = group, 
                  y = factor(anno, levels = sort(unique(anno))),
                  fill = mean.exp)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = col_palette,
    name = "Mean VST\nExpression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  labs(
    title = "Mean Expression of Candidate Genes Across Treatments",
    caption = "KCTR: Resistant; KSS: Susceptible; NT: Non-treated; D: 2,4-D; M: Mesotrione; A: Atrazine; C: Chlorsulfuron; L: Lactofen"
  )

# Save heatmap
ggsave("candidate_gene_heatmap.pdf", width = 10, height = 8, dpi = 300)
ggsave("candidate_gene_heatmap.png", width = 10, height = 8, dpi = 300)

