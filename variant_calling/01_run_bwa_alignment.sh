#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="MxO RAPiD BWA alignment"
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
module load bwa
module load samtools

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

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Use the sample file to create a vector of sample names
SAMPLENAMES=$(cut -d \t -f 1 $SAMPLEFILE)

# Iterate over the sample names
for SAMPLE in $SAMPLENAMES; do
  # Create a RG tag
  RG="@RG\tID:$SAMPLE\tSM:$SAMPLE"

  # Find the FASTQ files in the input directory that match the sample name
  SAMPLEFASTQS=$(find $FASTQDIR -name "*$SAMPLE*")

  ## ALIGNMENT TO STEVENS
  # Create the output SAM file name
  OUTPUT=$ALIGNDIR/${SAMPLE}_STEVENS_alignment.bam
  # Run the alignment
  bwa mem -t $NTHREADS -R $RG $DBPREFIXSTEVENS $SAMPLEFASTQS | samtools view -b -o $OUTPUT -

  ## ALIGNMENT TO OXY
  OUTPUT=$ALIGNDIR/${SAMPLE}_OXY_alignment.bam
  # Run the alignment
  bwa mem -t $NTHREADS -R $RG $DBPREFIXOXY $SAMPLEFASTQS | samtools view -b -o $OUTPUT -

done
