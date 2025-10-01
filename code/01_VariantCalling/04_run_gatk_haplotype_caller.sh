#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --job-name="MxO variant calling - haplotype caller"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


## 
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 3. Variant calling with GATK haplotype caller
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load gatk
module load parallel

#####################
## Set variables
#####################

# Working directory
WD=/project/gifvl_vaccinium/MxO_GeneticMap
VARDIR=$WD/results/variant_calling

# Path to directory containing the BAM files
BAMDIR=$VARDIR/alignment/

# Prefix of the indexed reference genomes
BLREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/Vaccinium_macrocarpon_BenLear_v2.fasta
STREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/V_macrocarpon_Stevens_v1.fasta
OXREFPREFIX=//project/gifvl_vaccinium/MxO_GeneticMap/results/genome_scaffolding_ragtag/Voxycoccos_NJ96-20_v1_ragtag_scaffolded/Voxycoccos_NJ96-20_v1_ragtag_scaffold.fasta

# Directory to output vcf files
VARIANTDIR=$VARDIR/variants/gatk/


# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Make subdirectories
mkdir -p $VARIANTDIR
# Make a directory for the haplotype caller output
HAPLODIR=$VARIANTDIR/haplotype_caller
if [ ! -d $HAPLODIR ]; then
  mkdir $HAPLODIR
fi

  
# Create a function to run haplotype caller
# Input: 1) BAM file, 2) output directory, 3) reference
run_haplotype_caller() {
  file=$1
  output=$2
  ref=$3

  # Output filename
  file1=$(basename $file)
  newfile="${file1%.bam}.g.vcf.gz"
  outfile=$output/$newfile

  # Run haplotype caller
  gatk HaplotypeCaller -R $ref -I $file -O $outfile -ERC GVCF

}

# Export the function
export -f run_haplotype_caller

## Run for each reference genome ##

# Subset the "BenLear" alignment files
BENLEARBAMFILES=($(find $BAMDIR -name "*BenLear*.bam"))
# Run the function in parallel
parallel -j $NTHREADS run_haplotype_caller {} $HAPLODIR $BLREFPREFIX ::: ${BENLEARBAMFILES[@]}


# Subset the "Stevens" alignment files
STEVENSBAMFILES=($(find $BAMDIR -name "*Stevens*.bam"))
# Run the function in parallel
parallel -j $NTHREADS run_haplotype_caller {} $HAPLODIR $STREFPREFIX ::: ${STEVENSBAMFILES[@]}


# Subset the "Oxycoccos" alignment files
OXYCOCCOSBAMFILES=($(find $BAMDIR -name "*Oxycoccos*.bam"))
# Run the function in parallel
parallel -j $NTHREADS run_haplotype_caller {} $HAPLODIR $OXREFPREFIX ::: ${OXYCOCCOSBAMFILES[@]}

# End of script