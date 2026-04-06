#!/bin/bash
#===============================================================================
# Script: 02_trimmomatic_paired_end.sh
# Description: Trim paired-end Illumina reads using Trimmomatic. This script:
#              1. Removes adapter sequences (Illumina TruSeq3)
#              2. Performs sliding window quality trimming
#              3. Removes reads shorter than 36 bp
#
# Usage: sbatch 02_trimmomatic_paired_end.sh <forward_reads.fastq.gz> <reverse_reads.fastq.gz> <output_directory>
#
# Arguments:
#   $1 - Path to forward reads (R1) file (.fastq.gz)
#   $2 - Path to reverse reads (R2) file (.fastq.gz)
#   $3 - Path to output directory (must exist before running)
#
# Output files (in specified output directory):
#   - paired_<R1_filename>.fastq.gz   : Paired forward reads (use this)
#   - unpaired_<R1_filename>.fastq.gz : Unpaired forward reads (orphans)
#   - paired_<R2_filename>.fastq.gz   : Paired reverse reads (use this)
#   - unpaired_<R2_filename>.fastq.gz : Unpaired reverse reads (orphans)
#
# Dependencies:
#   - Trimmomatic v0.39
#   - Adapter file: TruSeq3-SE.fa (adjust path as needed)
#
# Trimmomatic Parameters:
#   - ILLUMINACLIP: Remove adapters with 2 seed mismatches, 30 palindrome threshold,
#                   10 simple clip threshold, 3 min adapter length, keepBothReads=TRUE
#   - SLIDINGWINDOW:4:15 - Cut when 4-base average quality drops below 15
#   - MINLEN:36 - Drop reads shorter than 36 bp
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM job configuration
#-------------------------------------------------------------------------------
#SBATCH --job-name=trim-nebnext
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00           # Maximum 1 day runtime

#-------------------------------------------------------------------------------
# Parse command line arguments
#-------------------------------------------------------------------------------
FORWARD=$1                          # Forward reads (R1)
REVERSE=$2                          # Reverse reads (R2)
OUTPUT_DIR=$3                       # Output directory

# Extract base filenames (without path and .fastq.gz extension)
FORWARD_PREFIX=$(basename "$FORWARD" .fastq.gz)
REVERSE_PREFIX=$(basename "$REVERSE" .fastq.gz)

#-------------------------------------------------------------------------------
# Load Trimmomatic module and run
#-------------------------------------------------------------------------------
module load Trimmomatic

# Run Trimmomatic in paired-end mode
# Note: Adjust the adapter file path ($EBROOTTRIMMOMATIC or custom location)
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
    $FORWARD $REVERSE \
    $OUTPUT_DIR/paired_"$FORWARD_PREFIX".fastq.gz \
    $OUTPUT_DIR/unpaired_"$FORWARD_PREFIX".fastq.gz \
    $OUTPUT_DIR/paired_"$REVERSE_PREFIX".fastq.gz \
    $OUTPUT_DIR/unpaired_"$REVERSE_PREFIX".fastq.gz \
    ILLUMINACLIP:/path/to/adapters/TruSeq3-SE.fa:2:30:10:3:TRUE \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

#-------------------------------------------------------------------------------
# Notes:
# - Update ILLUMINACLIP path to match your adapter file location
# - Paired files should be used for downstream analysis
# - Unpaired files contain orphan reads where mate was discarded
# - Run FastQC on trimmed files to verify quality improvement
#-------------------------------------------------------------------------------
