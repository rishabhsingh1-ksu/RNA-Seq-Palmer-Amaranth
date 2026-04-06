# RNA-Seq Analysis of Herbicide Resistance in *Amaranthus palmeri*

[![DOI]()
[![GEO](https://img.shields.io/badge/GEO-GSE326567-orange.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326567)

## Overview

This repository contains the bioinformatics pipeline and analysis scripts for the RNA-Seq study: **"Deciphering Metabolic Resistance to Multiple Herbicides in *Amaranthus palmeri* via Transcriptome Analysis"** (Singh et al., The Plant Journal).

The study identifies cytochrome P450s (*CYP72A219*, *CYP704B1*), glutathione-*S*-transferases (GSTs), and other genes involved in metabolic resistance to herbicides in Palmer amaranth (*Amaranthus palmeri*).

## Study Design

### Biological System
- **Species**: *Amaranthus palmeri* S. Watson (Palmer amaranth)
- **Populations**: 
  - Resistant: KCTR/KCTR-G2 (Kansas Conservation Tillage Resistant, resistant to 6 herbicide sites of action)
  - Susceptible: KSS (Kansas Susceptible)

### Herbicide Treatments
| HRAC Group | Herbicide | Site of Action | Dose (g ha⁻¹) |
|------------|-----------|----------------|---------------|
| 2 | Chlorsulfuron | ALS-inhibitor | 18 |
| 4 | 2,4-D | Synthetic auxin | 560 |
| 5 | Atrazine | PSII-inhibitor | 2240 |
| 14 | Lactofen | PPO-inhibitor | 175 |
| 27 | Mesotrione | HPPD-inhibitor | 105 |

### Experimental Design
- 3 biological replicates per population
- Leaf samples collected 6 hours after herbicide treatment
- Non-treated controls included for each batch
- RNA sequencing: Paired-end Illumina reads

## Repository Structure

```
RNA-Seq-Palmer-Amaranth/
├── README.md
├── LICENSE
├── scripts/
│   ├── bash/
│   │   ├── 01_run_fastqc_on_directory.sh    # Quality control
│   │   ├── 02_trimmomatic_batch_submit.sh   # Batch trimming submission
│   │   ├── 03_build_STAR_index.sh           # Build genome index
│   │   ├── 04_STAR_align.sh                 # Align reads
│   └── R/
│       ├── 01_DESeq2_differential_expression.R  # DGE analysis
│       ├── 02_VST_normalization_heatmap.R       # Normalization & visualization
│       └── 03_WGCNA_coexpression_analysis.R     # Co-expression network
├── data/
│   ├── raw_conts_star.csv
|   ├── metadata1+2+3.csv
|   └── vst_norm_counts.xlsx
├── docs/
│   └── methods_supplementary.md
└── results/
    └── Discussed in the mansucript
```

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                      ANALYSIS PIPELINE                          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  1. QUALITY CONTROL                                             │
│     • FastQC (v0.11.7) - Raw read quality assessment            │
│     • MultiQC - Aggregate QC reports                            │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  2. READ TRIMMING                                               │
│     • Trimmomatic (v0.39) - Adapter removal & quality trimming  │
│     • Parameters: ILLUMINACLIP, SLIDINGWINDOW:4:15, MINLEN:36   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  3. ALIGNMENT & QUANTIFICATION                                  │
│     • STAR aligner - Map to A. palmeri genome                   │
│     • Reference: Palmer amaranth Hap1 (JBEFMX000000000)         │
│     • --quantModeGeneCounts for read counting                   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  4. DIFFERENTIAL EXPRESSION                                     │
│     • DESeq2 - Pairwise DGE analysis (R vs S)                  │
│     • Criteria: padj < 0.05, |log2FC| ≥ 2                       │
│     • Contrasts: R vs S for each treatment and non-treated                   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  5. CO-EXPRESSION ANALYSIS (WGCNA)                              │
│     • Variance stabilizing transformation (VST)                 │
│     • Signed network, soft-threshold power = 14                 │
│     • Module-trait correlation analysis                         │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  6. VISUALIZATION & ANNOTATION                                  │
│     • Heatmaps of candidate gene expression                     │
│     • Functional annotation of DEGs                             │                                      │
└─────────────────────────────────────────────────────────────────┘
```

## Software Requirements

### Command Line Tools
| Software | Version | Purpose |
|----------|---------|---------|
| FastQC | 0.11.7 | Quality control |
| MultiQC | 1.9+ | QC aggregation |
| Trimmomatic | 0.39 | Read trimming |
| STAR | 2.7+ | Read alignment |
| Java | 1.8+ | Required by Trimmomatic |

### R Packages
| Package | Version | Purpose |
|---------|---------|---------|
| DESeq2 | 1.34+ | Differential expression |
| WGCNA | 1.70+ | Co-expression network |
| tidyverse | 1.3+ | Data manipulation |
| pheatmap | 1.0.12 | Heatmap visualization |
| CorLevelPlot | 0.99+ | Module-trait correlations |
| readxl | 1.4+ | Read Excel files |
| writexl | 1.4+ | Write Excel files |

### Installation

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "WGCNA", "CorLevelPlot"))

# Install CRAN packages
install.packages(c("tidyverse", "pheatmap", "readxl", "writexl", "patchwork"))
```

## Usage

### 1. Quality Control

```bash
# Create output directory
mkdir fastqc_raw

# Run FastQC on all samples
sbatch scripts/bash/01_run_fastqc_on_directory.sh \
    /path/to/raw_reads \
    /path/to/fastqc_raw

# Aggregate with MultiQC
module load MultiQC
cd fastqc_raw
multiqc .
```

### 2. Read Trimming

```bash
# Create output directory
mkdir trimmed_reads

# Run batch trimming
sbatch scripts/bash/03_trimmomatic_batch_submit.sh \
    /path/to/raw_reads \
    /path/to/trimmed_reads

# Verify with post-trimming FastQC
mkdir fastqc_trimmed
sbatch scripts/bash/01_run_fastqc_on_directory.sh \
    /path/to/trimmed_reads \
    /path/to/fastqc_trimmed
```

### 3. Alignment (STAR)

```bash
# Create genome index directory
mkdir star_index

# Build STAR index (one-time setup)
# For 150bp reads, use sjdbOverhang = 149
sbatch scripts/bash/04_build_STAR_index.sh \
    /path/to/star_index \
    /path/to/Apalmeri_genome.fa \
    /path/to/Apalmeri.gtf \
    Apalmeri \
    149

# Align all samples (batch mode)
mkdir aligned
bash scripts/bash/06_STAR_align_batch_submit.sh \
    /path/to/star_index \
    /path/to/trimmed_reads \
    /path/to/aligned
```

### 4. Differential Expression Analysis

```r
# Open R and run the DESeq2 script
source("scripts/R/01_DESeq2_differential_expression.R")
```

### 5. WGCNA Analysis

```r
# Run WGCNA after DESeq2
source("scripts/R/03_WGCNA_coexpression_analysis.R")
```

## Key Results

### Differentially Expressed Genes

| Treatment | DEGs (R vs S) |
|-----------|---------------|
| Chlorsulfuron | 414 |
| 2,4-D | 129 |
| Atrazine | 529 |
| Mesotrione | 152 |
| Lactofen | 688 |

### Candidate Resistance Genes

| Gene | Annotation | Function |
|------|------------|----------|
| *CYP72A219* | Cytochrome P450 | Phase I metabolism |
| *CYP704B1* | Cytochrome P450 | Phase I metabolism |
| GST C-terminal | Glutathione-S-transferase | Phase II conjugation |
| ABCG/ABCC family | ABC transporters | Phase III transport |
| UGT79B7/UGT76F1 | UDP-glycosyltransferases | Phase II conjugation |

## Data Availability

### Raw Sequencing Data
All sequencing data is publicly available at NCBI GEO:
- **Accession**: [GSE326567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326567)

### Reference Genome
- **Palmer amaranth Hap1 genome**: [JBEFMX000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBEFMX000000000)
- Reference: Raiyemo et al., 2025

### Required Input Files (Not Included)
To run these analyses, you will need:

1. **Count matrix**: `raw_counts_STAR.csv`
   - Gene count matrix with genes as rows and samples as columns
   - First column: GeneID

2. **Sample metadata**: `metadata1+2+3.csv`
   - Columns: sample, group, batch
   - group format: `{population}_{treatment}` (e.g., KCTR_NT, KSS_24D)


## Citation

If you use this code or data, please cite:

```bibtex

}
```

## Authors

- **Rishabh Singh** - Kansas State University / University of Illinois Urbana-Champaign ([rs81@illinois.edu](mailto:rs81@illinois.edu))
- **Yaiphabi Kumam** - University of Florida
- **Mohit Mahey** - Michigan State University
- **Eric Patterson** - Michigan State University
- **Sanzhen Liu** - Kansas State University
- **Sarah Lancaster** - Kansas State University
- **Mithila Jugulam** (Corresponding) - Texas A&M AgriLife ([m.jugulam@ag.tamu.edu](mailto:m.jugulam@ag.tamu.edu))

## Acknowledgements

This work was funded by the Kansas Soybean Commission.

Bioinformatics scripts were developed with assistance from Teresa Shippy at the Kansas State University Bioinformatics Center.


## Contact

For questions about the analysis pipeline, please open an issue or contact:
- Rishabh Singh: [rs81@illinois.edu](mailto:rs81@illinois.edu)
- Mithila Jugulam: [m.jugulam@ag.tamu.edu](mailto:m.jugulam@ag.tamu.edu)
