#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=08:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="MxO RAPiD variant calling pipeline - read cleaning"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


##
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 2. Read trimming and quality control
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
SAMPLEFILE=$WD/data/mxo_rapid15k_sample_list.txt

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

# Create a "post-QC" directory in the QC output directory
POSTQCOUTPUT=$QCOUTPUT/post_qc
if [ ! -d $POSTQCOUTPUT ]; then
  mkdir -p $POSTQCOUTPUT
fi

# Export FASTQDIR variable for use in parallel
export FASTQDIRUSE=$FASTQDIR

# Use the sample file to create a vector of sample names
SAMPLENAMES=($(cut -d \t -f 1 $SAMPLEFILE))


## Step 1: Read trimming and quality control

# Create the output directory if it does not exist
if [ ! -d $CLEANEDFASTQOUTPUT ]; then
  mkdir -p $CLEANEDFASTQOUTPUT
fi

# Write a function that runs cutadapt in parallel
run_cutadapt() {

  # The first input will be the prefix of the two fastq files
  prefix=$1
  # The second input will be the output directory
  output=$2

  # Use the prefix to find the two FASTQ files
  file1=$(find $FASTQDIRUSE -name "*${prefix}_*_R1_*.fastq.gz")
  file2=$(find $FASTQDIRUSE -name "*${prefix}_*_R2_*.fastq.gz")

  # Run cutadapt in paired-end mode
  file1_basename=$(basename "$file1")
  file2_basename=$(basename "$file2")
  outputfile1=$output/${file1_basename%.fastq.gz}_trim.fastq.gz
  outputfile2=$output/${file2_basename%.fastq.gz}_trim.fastq.gz
  logfile=$output/"${prefix}_cutadapt_log.txt"

  cutadapt -a "$ADAPTER" -A "$ADAPTER" -q 20,20 -m 30 -M 0 --pair-filter both \
  -o $outputfile1 -p $outputfile2 $file1 $file2 > $logfile 2>&1

}

# Export the function
export -f run_cutadapt

# Run the function in parallel
parallel -j $NTHREADS run_cutadapt {} $CLEANEDFASTQOUTPUT ::: ${SAMPLENAMES[@]}


# Step 2: post-trimming quality control

# Iterate over the sample names
for SAMPLE in ${SAMPLENAMES[@]}; do
  echo -e "\nRunning post-trimming fastqc for sample: $SAMPLE"

  # Find the FASTQ files in the input directory that match the sample name;
  # There will be two files per sample (R1 and R2)
  # Use a wildcard before and after the sample name to capture the full file name
  # (e.g. S1_L001_R1_001.fastq.gz)
  # Store the file names in a variable
  SAMPLEFASTQS=$(find $CLEANEDFASTQOUTPUT -name "*$SAMPLE*.fastq.gz")

  # Run FastQC on the raw reads
  fastqc -f fastq -t $NTHREADS -o $POSTQCOUTPUT ${SAMPLEFASTQS[@]}

done

# Run multiqc to summarize the FastQC results
multiqc -o $POSTQCOUTPUT $POSTQCOUTPUT