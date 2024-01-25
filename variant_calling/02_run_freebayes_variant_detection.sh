#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="MxO RAPiD freebayes variant calling"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 1. Alignment to reference genomes using the BWA aligner
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load samtools
module load freebayes

#####################
## Set variables
#####################

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap

# Name of input directory
INPUT=/project/gifvl_vaccinium/cranberryGenotyping/RAPiD_Cranberry_15K/Data/2022
# Name of input directory with FASTQ
FASTQDIR=$INPUT/fastq_files

# Name of the file containing sample names
SAMPLEFILE=$WD/data/mxo_rapid15k_sample_list.txt

# Prefix of the indexed Stevens reference genome
DBPREFIXSTEVENS=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_macrocarpon_Stevens_v1.fasta
# Prefix of the indexed oxycoccos reference genome
DBPREFIXOXY=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_macrocarpon_Stevens_v1.fasta

# Directory to output alignments
ALIGNDIR=$WD/variant_calling/alignment/

# Directory to output variants
VARIANTDIR=$WD/variant_calling/variants/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# List all of the alignment files in the alignment directory
# First find those aligned to STEVENS
ALIGNMENTFILESSTE=$(find $VARIANTDIR -name "*STEVENS_alignment.bam")
# Next find those aligned to OXY
ALIGNMENTFILESOXY=$(find $VARIANTDIR -name "*OXY_alignment.bam")


## Run variant calling for alignment to STEVENS
# Ouput file
OUTPUT=$VARIANTDIR/MXO_STEVENS_ref_variants.vcf
freebayes -f $DBPREFIXSTEVENS $ALIGNMENTFILESSTE > $OUTPUT

## Run variant calling for alignment to OXY
# Ouput file
OUTPUT=$VARIANTDIR/MXO_OXY_ref_variants.vcf
freebayes -f $DBPREFIXOXY $ALIGNMENTFILESOXY > $OUTPUT
