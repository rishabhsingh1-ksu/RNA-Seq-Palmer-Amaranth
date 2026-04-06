#!/bin/bash -l
#===============================================================================
# Script: 03_build_STAR_index.sh
# Description: Build STAR genome index for Amaranthus palmeri reference genome.
#              This is a one-time setup step required before alignment.
#              The index allows STAR to efficiently map RNA-Seq reads.
#
# Usage: sbatch 04_build_STAR_index.sh <genome_dir> <genome_fasta> <gtf_file> <output_prefix> <read_length-1>
#
# Arguments:
#   $1 - Path to output directory for genome index (must exist)
#   $2 - Path to unzipped reference genome FASTA file (.fa or .fasta)
#   $3 - Path to gene annotation GTF file
#   $4 - Output filename prefix for log files
#   $5 - Read length minus 1 (e.g., for 150bp reads, use 149)
#
# Example:
#   sbatch 04_build_STAR_index.sh \
#       /path/to/star_index \
#       /path/to/Apalmeri_genome.fa \
#       /path/to/Apalmeri.gtf \
#       Apalmeri_index \
#       149
#
# Reference Genome:
#   - Amaranthus palmeri Hap1 assembly
#   - NCBI Accession: JBEFMX000000000
#   - Reference: Raiyemo et al., 2025
#
# Output:
#   - STAR index files in the specified genome directory
#   - Log files with specified prefix
#
# Notes:
#   - This step requires significant memory (~32GB for plant genomes)
#   - Only needs to be run once per genome/annotation combination
#   - Index must be rebuilt if genome or GTF is updated
#
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM job configuration
#-------------------------------------------------------------------------------
#SBATCH --mem-per-cpu=6G            # Memory per CPU core
#SBATCH --cpus-per-task=8           # Number of CPU cores
#SBATCH --time=1-00:00:00           # Maximum runtime (1 day)
#SBATCH --job-name=STAR_index       # Job name for queue

#-------------------------------------------------------------------------------
# Load required modules
#-------------------------------------------------------------------------------
module load STAR

#-------------------------------------------------------------------------------
# Build STAR genome index
#-------------------------------------------------------------------------------
# Parameters explained:
#   --runMode genomeGenerate    : Index building mode
#   --genomeDir                 : Output directory for index files
#   --genomeFastaFiles          : Reference genome FASTA
#   --sjdbGTFfile               : Gene annotation for splice junction database
#   --sjdbOverhang              : Read length - 1 (optimal for splice junction detection)
#   --outFileNamePrefix         : Prefix for output log files
#   --runThreadN                : Number of threads (match --cpus-per-task)
#   --genomeSAindexNbases       : Index parameter (13 works for most plant genomes)
#                                 Use 12-13 for genomes ~500Mb-1Gb

STAR --runMode genomeGenerate \
    --genomeDir $1 \
    --genomeFastaFiles $2 \
    --sjdbGTFfile $3 \
    --sjdbOverhang $5 \
    --outFileNamePrefix $4 \
    --runThreadN 8 \
    --genomeSAindexNbases 13

#-------------------------------------------------------------------------------
# Post-run checklist:
# 1. Verify index files were created in genomeDir:
#    - chrLength.txt, chrNameLength.txt, chrName.txt, chrStart.txt
#    - Genome, genomeParameters.txt, SA, SAindex
#    - sjdbInfo.txt, sjdbList.fromGTF.out.tab, sjdbList.out.tab
# 2. Check Log.out for any warnings or errors
# 3. Note the genome version and GTF version used for reproducibility
#-------------------------------------------------------------------------------
