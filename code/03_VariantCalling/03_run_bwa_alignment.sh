#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --job-name="MxO RAPiD BWA alignment"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## GBS variant calling pipeline
##
## Step 3. Alignment to reference genomes using the BWA aligner
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
WD=/project/gifvl_vaccinium/MxO_GeneticMap
TEMPDIR=/90daydata/gifvl_vaccinium/MxO_GeneticMap

VARDIR=$WD/results/variant_calling

# Path to directory containing the FASTQ files
FASTQDIR=$WD/results/variant_calling/cleaned_fastq_files/


# Prefix of the indexed reference genomes
BLREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/Vmacrocarpon_BenLear_v1-scaffolds.fasta
STREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/V_macrocarpon_Stevens_v1.fasta
OXREFPREFIX=/project/gifvl_vaccinium/MxO_GeneticMap/results/genome_scaffolding_ragtag/Voxycoccos_NJ96-20_v1_ragtag_scaffolded/Voxycoccos_NJ96-20_v1_ragtag_scaffold.fasta


# Directory to output alignment files
# Putting the output in 90Daydata to save space on project directory
ALIGNDIR=$TEMPDIR/variant_calling/alignment/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################
 

## Run the pipeline

# Change working directory
cd $WD

# Create the aligndir if it does not exist
mkdir -p $ALIGNDIR
  

## Get a list of the cleaned FASTQ files
FASTQFILES=$(find $FASTQDIR -name "*trim.fastq.gz")
# Remove the suffix "_R[1,2]_001_trim.fastq.gz" to get the unique sample names
SAMPLES=$(echo "$FASTQFILES" | sed 's/_R[1,2]_001_trim.fastq.gz//g' | sort -u)

# Iterate over those files and align to the BenLear, Stevens, and Oxycoccos reference genomes
# Use paired-end mode for BWA mem
for SAMPLE in $SAMPLES; do
  # Get the FASTQ files for this sample
  fastq1file=${SAMPLE}_R1_001_trim.fastq.gz
  fastq2file=${SAMPLE}_R2_001_trim.fastq.gz

  # Basename of sample
  SAMPLEBASE=$(basename $SAMPLE)

  # Create a RG tag
  RG="@RG\tID:$SAMPLEBASE\tSM:$SAMPLEBASE\tPL:ILLUMINA"

  # Align to the BenLear reference
  REFPREFIX=$BLREFPREFIX
  REFNAME="BenLear"
  OUTPUTBAM=$ALIGNDIR/${SAMPLEBASE}_${REFNAME}_alignment.bam
  
  # Run the alignment in a pipeline
  bwa mem -t $NTHREADS -R $RG $REFPREFIX $fastq1file $fastq2file | \
  samtools fixmate -u -m - - | \
  samtools sort -@ $NTHREADS -o $OUTPUTBAM -
  samtools index $OUTPUTBAM

  # Align to the Stevens reference
  REFPREFIX=$STREFPREFIX
  REFNAME="Stevens"
  OUTPUTBAM=$ALIGNDIR/${SAMPLEBASE}_${REFNAME}_alignment.bam
  # Run the alignment in a pipeline
  bwa mem -t $NTHREADS -R $RG $REFPREFIX $fastq1file $fastq2file | \
  samtools fixmate -u -m - - | \
  samtools sort -@ $NTHREADS -o $OUTPUTBAM -
  samtools index $OUTPUTBAM 

  # Align to the Oxycoccos reference
  REFPREFIX=$OXREFPREFIX
  REFNAME="Oxycoccos"
  OUTPUTBAM=$ALIGNDIR/${SAMPLEBASE}_${REFNAME}_alignment.bam
  # Run the alignment in a pipeline
  bwa mem -t $NTHREADS -R $RG $REFPREFIX $fastq1file $fastq2file | \
  samtools fixmate -u -m - - | \
  samtools sort -@ $NTHREADS -o $OUTPUTBAM -
  samtools index $OUTPUTBAM

done

