#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=08:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="MxO RAPiD variant calling pipeline - Fast QC"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


##
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 1. Quality control of raw reads
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load fastqc
module load cutadapt
module load parallel
module load multiqc

#####################
## Set variables
#####################

# Working directory
WD=/project/gifvl_vaccinium/MxO_GeneticMap

# Name of input directory
INPUT=$WD/data_raw/
# Name of input directory with FASTQ
FASTQDIR=$INPUT/fastq_files

# Name of the file containing sample names
SAMPLEFILE=$WD/data/population_metadata.csv

# Adapter sequence to remove
# 'export' needs to be used in order to put the variable in the global environment
export ADAPTER="AGATCGGAAGAGC"


# Directory to output QC results
QCOUTPUT=$WD/results/variant_calling/qc/
# Directory to output cleaned FASTQ files
CLEANEDFASTQOUTPUT=$WD/results/variant_calling/cleaned_fastq_files/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##########################################
## Run the pipeline - DO NOT EDIT BELOW
##########################################

## Step 1: Quality control of raw reads

# Change working directory
cd $WD

# Create the output directory if it does not exist
if [ ! -d $QCOUTPUT ]; then
  mkdir -p $QCOUTPUT
fi 

# Create a "pre-QC" directory in the QC output directory
PREQCOUTPUT=$QCOUTPUT/pre_qc
if [ ! -d $PREQCOUTPUT ]; then
  mkdir -p $PREQCOUTPUT
fi

# Export FASTQDIR variable for use in parallel
export FASTQDIRUSE=$FASTQDIR

# Use the sample file to create a vector of sample names
dos2unix "$SAMPLEFILE"
SAMPLENAMES=($(cut -d , -f 8 "$SAMPLEFILE" | grep -v '^[[:space:]]*$' | grep -v 'RG_Sample_Code'))


# Iterate over the sample names
for SAMPLE in $SAMPLENAMES; do
  echo -e "\nProcessing sample: $SAMPLE"

  # Find the FASTQ files in the input directory that match the sample name;
  # There will be two files per sample (R1 and R2)
  # Use a wildcard before and after the sample name to capture the full file name
  # (e.g. S1_L001_R1_001.fastq.gz)
  # Store the file names in a variable
  SAMPLEFASTQS=$(find $FASTQDIR -name "*$SAMPLE*")

  # Run FastQC on the raw reads
  fastqc -f fastq -t $NTHREADS -o $PREQCOUTPUT ${SAMPLEFASTQS[@]}

done

# Run multiqc to summarize the FastQC results
multiqc -o $PREQCOUTPUT $PREQCOUTPUT

