#!/bin/bash -l
#===============================================================================
# Script: 04_STAR_align.sh
# Description: Align paired-end RNA-Seq reads to the A. palmeri reference genome
#              using STAR aligner. Simultaneously generates gene-level read counts
#              using the --quantMode GeneCounts option.
#
# Usage: sbatch 04_STAR_align.sh <genome_index_dir> <read1.fastq.gz> <read2.fastq.gz> <output_prefix>
#
# Arguments:
#   $1 - Path to STAR genome index directory (created by 04_build_STAR_index.sh)
#   $2 - Path to forward reads (R1) file (.fastq.gz)
#   $3 - Path to reverse reads (R2) file (.fastq.gz)
#   $4 - Output filename prefix (sample name recommended)
#
# Example:
#   sbatch 04_STAR_align.sh \
#       /path/to/star_index \
#       /path/to/trimmed/KCTR1_R1.fastq.gz \
#       /path/to/trimmed/KCTR1_R2.fastq.gz \
#       KCTR1_
#
# Output files (with specified prefix):
#   - *Aligned.out.sam     : Aligned reads in SAM format
#   - *ReadsPerGene.out.tab: Gene-level read counts (use for DESeq2)
#   - *Log.final.out       : Alignment statistics summary
#   - *Log.out             : Detailed run log
#   - *Log.progress.out    : Progress during run
#   - *SJ.out.tab          : Splice junctions detected
#
# ReadsPerGene.out.tab format:
#   Column 1: Gene ID
#   Column 2: Unstranded counts (use this for standard RNA-Seq)
#   Column 3: Forward strand counts (for strand-specific libraries)
#   Column 4: Reverse strand counts (for strand-specific libraries)
#
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM job configuration
#-------------------------------------------------------------------------------
#SBATCH --mem-per-cpu=6G            # Memory per CPU core (96GB total)
#SBATCH --cpus-per-task=16          # Number of CPU cores
#SBATCH --time=1-00:00:00           # Maximum runtime (1 day)
#SBATCH --job-name=STAR_align       # Job name for queue

#-------------------------------------------------------------------------------
# Load required modules
#-------------------------------------------------------------------------------
module load STAR

#-------------------------------------------------------------------------------
# Run STAR alignment with gene counting
#-------------------------------------------------------------------------------
# Parameters explained:
#   --genomeDir          : Path to pre-built STAR index
#   --readFilesIn        : Input FASTQ files (R1 R2 for paired-end)
#   --readFilesCommand   : Command to decompress gzipped files
#   --quantMode GeneCounts: Generate gene-level counts simultaneously
#   --outFileNamePrefix  : Prefix for all output files
#   --runThreadN         : Number of threads (match --cpus-per-task)

STAR --genomeDir $1 \
    --readFilesIn $2 $3 \
    --readFilesCommand gunzip -c \
    --quantMode GeneCounts \
    --outFileNamePrefix $4 \
    --runThreadN 16

#-------------------------------------------------------------------------------
# Post-alignment steps:
# 1. Check Log.final.out for alignment statistics:
#    - Uniquely mapped reads % (should be >70% for good data)
#    - % of reads mapped to multiple loci
#    - % unmapped reads
#
# 2. Combine all *ReadsPerGene.out.tab files into a count matrix:
#    - Extract column 1 (GeneID) and column 2 (unstranded counts)
#    - Merge across all samples for DESeq2 input
#
# 3. Optional: Convert SAM to sorted BAM for visualization:
#    module load SAMtools
#    samtools sort -@ 8 -o sample.sorted.bam sample_Aligned.out.sam
#    samtools index sample.sorted.bam
#-------------------------------------------------------------------------------
