# Supplementary Methods

## Detailed Bioinformatics Pipeline

This document provides additional methodological details for the RNA-Seq analysis pipeline.

---

## 1. Quality Control with FastQC

### Key Parameters Monitored
- **Per Base Sequence Quality**: Should remain close to Phred score 36 throughout read length
- **Adapter Content**: Illumina Universal Adapter typically appears at read 3' ends
- **Sequence Length Distribution**: Post-trimming should show trimmed read lengths
- **GC Content**: Should match expected species GC content (~35-40% for *A. palmeri*)
- **Overrepresented Sequences**: May indicate adapter contamination or rRNA

### Interpretation Guidelines
| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| Per Base Quality | >30 | 20-30 | <20 |
| Adapter Content | <1% | 1-5% | >5% |
| Duplication Level | <20% | 20-50% | >50% |

---

## 2. Read Trimming with Trimmomatic

### Trimmomatic Parameters Explained

```
ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:3:TRUE
```
- `2`: Seed mismatches allowed
- `30`: Palindrome clip threshold  
- `10`: Simple clip threshold
- `3`: Minimum adapter length for detection
- `TRUE`: Keep both reads even if only one adapter found

```
SLIDINGWINDOW:4:15
```
- `4`: Window size (bases)
- `15`: Required average quality within window

```
MINLEN:36
```
- Minimum read length to retain after trimming

### Expected Output
- **Paired files**: Use for downstream analysis (both mates passed QC)
- **Unpaired files**: Orphan reads (mate was discarded) - typically not used

---

## 3. Alignment with STAR

### Reference Genome
- **Genome**: *Amaranthus palmeri* Hap1 assembly
- **Accession**: JBEFMX000000000
- **Reference**: Raiyemo et al., 2025

### Recommended STAR Parameters

```bash
STAR --genomeDir /path/to/index \
     --readFilesIn R1.fastq.gz R2.fastq.gz \
     --readFilesCommand zcat \
     --quantMode GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterMismatchNoverReadLmax 0.04 \
     --outFilterMultimapNmax 20 \
     --alignIntronMin 20 \
     --alignIntronMax 500000 \
     --runThreadN 8
```

### Output Files
- `*ReadsPerGene.out.tab`: Gene-level read counts (use column 2 for unstranded libraries)
- `*Aligned.sortedByCoord.out.bam`: Aligned reads
- `*Log.final.out`: Alignment statistics

---

## 4. Differential Expression with DESeq2

### Statistical Model
```r
design = ~ group    # For simple comparisons
design = ~ batch + group    # When accounting for batch effects
```

### Normalization Method
- **Size factors**: Median-of-ratios method
- Accounts for both library size and RNA composition

### Criteria for Candidate Gene Selection
1. **Statistical significance**: Adjusted p-value (Benjamini-Hochberg) < 0.05
2. **Biological significance**: |log2 fold change| ≥ 2
3. **Expression filter**: Minimum mean normalized counts > 10
4. **Functional relevance**: Associated with herbicide metabolism pathways
5. **Direction**: Upregulated in resistant (R) vs susceptible (S)

### Target Gene Families
- Cytochrome P450s (CYPs) - Phase I metabolism
- Glutathione-S-transferases (GSTs) - Phase II conjugation
- UDP-glycosyltransferases (UGTs) - Phase II conjugation  
- ABC transporters - Phase III transport

---

## 5. WGCNA Analysis

### Pre-processing Steps
1. **Outlier gene removal**: `goodSamplesGenes()` function
2. **Outlier sample removal**: Hierarchical clustering (average linkage)
3. **Low-count filtering**: Remove genes with <10 counts in <75% of samples
4. **Normalization**: Variance stabilizing transformation (VST)

### Network Construction Parameters
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Network type | Signed | Distinguish positive/negative correlations |
| Soft power | 14 | R² > 0.8 with reasonable connectivity |
| TOM type | Signed | Topological overlap matrix |
| Max block size | 8000 | Memory-dependent |
| Merge cut height | 0.25 | Module merging threshold |
| Min module size | 30 | Default, prevents very small modules |

### Soft Power Selection
Choose power where:
- Scale-free topology fit (R²) > 0.8
- Mean connectivity is not too low
- Higher R² is preferred if multiple powers qualify

### Module-Trait Correlation
- **Trait**: Binary resistance variable (KCTR=1, KSS=0)
- **Correlation**: Pearson correlation with module eigengenes
- **Significance**: Student asymptotic p-values

---

## 6. qRT-PCR Validation

### Primer Sequences

| Gene | Direction | Sequence (5'-3') | Amplicon (bp) |
|------|-----------|------------------|---------------|
| CYP72A219 | Forward | GTTTATATCCCTGGTTGGAGGT | 134 |
| CYP72A219 | Reverse | CCTTTGCCCTCTCTTCTCTATT | |
| CYP704B1 | Forward | CAGGAGTTGTTGATGAGGATGA | 95 |
| CYP704B1 | Reverse | GAGCAAAGGAGTTCTCAGGTAG | |

### Reference Gene
- β-tubulin or elongation factor (EF1α) commonly used for *Amaranthus*

---

## Troubleshooting Common Issues

### FastQC Issues
| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| Poor quality at 3' end | Normal for Illumina | Trimmomatic will handle |
| High adapter content | Short inserts | Verify ILLUMINACLIP parameters |
| High duplication | Over-sequencing or low complexity | May need to remove duplicates |

### DESeq2 Issues
| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| Size factor = 0 | Zero counts | Pre-filter low-count genes |
| Convergence warning | Low counts or outliers | Filter genes, check outliers |
| Few DEGs | Low power or high variance | Increase replicates |

### WGCNA Issues
| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| No modules | Power too low/high | Re-evaluate soft threshold |
| Too many grey genes | Power too high | Lower soft threshold |
| Memory errors | Block size too large | Reduce maxBlockSize |

---

## References

1. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data.
2. Bolger, A.M., et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
3. Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21.
4. Love, M.I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.
5. Langfelder, P. & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9(1), 559.
6. Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.
