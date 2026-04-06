#!/bin/bash
#===============================================================================
# Script: 01_run_fastqc_on_directory.sh
# Description: Run FastQC quality control on all paired-end FASTQ files in a 
#              directory. FastQC provides graphical overview of sequence quality
#              including per-base quality scores, adapter content, and GC content.
#
# Usage: sbatch 01_run_fastqc_on_directory.sh <path_to_directory_of_readfiles> <path_to_output_dir>
#
# Arguments:
#   $1 - Full path to directory containing .fastq.gz files
#   $2 - Full path to output directory for FastQC reports
#
# Output: HTML and ZIP files for each input FASTQ file
#
# Dependencies:
#   - FastQC v0.11.7 or higher
#   - Java 1.8 or higher
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM job configuration
#-------------------------------------------------------------------------------
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=fastqc
#SBATCH --mem=4G                    # Memory allocation per job
#SBATCH --time=24:00:00             # Maximum runtime (DD-HH:MM:SS format)

#-------------------------------------------------------------------------------
# Load required modules
#-------------------------------------------------------------------------------
module load Java/1.8.0_144
module load FastQC

#-------------------------------------------------------------------------------
# Run FastQC on all gzipped FASTQ files in the input directory
#-------------------------------------------------------------------------------
# Loop through all .fastq.gz files in the specified directory
# -o flag specifies output directory
for file in "$1"/*.fastq.gz; do 
    fastqc -o $2 $file
done

#-------------------------------------------------------------------------------
# Post-processing notes:
# 1. Download .html files to view in web browser
# 2. Pay attention to:
#    - Per Base Sequence Quality (should be close to 36 throughout read length)
#    - Adapter Content (Illumina Universal Adapter typically found at read ends)
#    - Total sequence counts in Basic Statistics
# 3. Use MultiQC to aggregate results: module load MultiQC && multiqc .
#-------------------------------------------------------------------------------
