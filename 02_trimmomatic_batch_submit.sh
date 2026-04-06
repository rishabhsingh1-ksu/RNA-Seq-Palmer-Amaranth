#!/bin/bash
#===============================================================================
# Script: 03_trimmomatic_batch_submit.sh
# Description: Batch wrapper script to run Trimmomatic on all paired-end read 
#              files in a directory. Automatically pairs R1 and R2 files and 
#              submits individual trimming jobs for each pair.
#
# Usage: sbatch 02_trimmomatic_batch_submit.sh <path_to_reads_directory> <path_to_output_directory>
#
# Arguments:
#   $1 - Full path to directory containing paired-end reads
#        Files must follow naming convention: *R1_001.fastq.gz and *R2_001.fastq.gz
#   $2 - Full path to output directory (must exist before running)
#
# Output: Trimmed paired and unpaired read files in the output directory
#
# Notes:
#   - Run this script from OUTSIDE the directory containing reads
#   - The output directory must already exist
#   - Script assumes R1/R2 naming convention; modify sed command if different
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM job configuration
#-------------------------------------------------------------------------------
#SBATCH --time=1-00:00:00           # Maximum 1 day runtime

#-------------------------------------------------------------------------------
# Process all R1 files and automatically find corresponding R2 files
#-------------------------------------------------------------------------------
# Loop through all R1 files in the specified directory
for file in "$1"/*R1_001.fastq.gz
do
    # Print current file being processed
    echo "Processing forward reads: $file"
    
    # Generate R2 filename by replacing R1 with R2
    r2=$(echo $file | sed 's/R1/R2/')
    echo "Corresponding reverse reads: $r2"
    
    # Submit trimmomatic job for this read pair
    # Update the path to your trimmomatic script as needed
    sbatch /path/to/scripts/02_trimmomatic_paired_end.sh $file $r2 $2
    # sleep 10
done

#-------------------------------------------------------------------------------
# Example usage:
# 1. Create output directory: mkdir trimmed_reads
# 2. Run from parent directory: 
#    sbatch 03_trimmomatic_batch_submit.sh /path/to/raw_reads /path/to/trimmed_reads
# 3. Monitor jobs: kstat -d 1
# 4. After completion, run FastQC on trimmed reads to verify quality
#-------------------------------------------------------------------------------
