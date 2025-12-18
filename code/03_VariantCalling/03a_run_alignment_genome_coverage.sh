#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --job-name="MxO RAPiD BWA alignment - genome coverage"
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
## Step 3a. Estimate genome coverage using the bam files from the alignment step
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load parallel
module load samtools
module load bedtools

#####################
## Set variables
#####################

# Working directory
WD=/project/gifvl_vaccinium/MxO_GeneticMap
TEMPDIR=/90daydata/gifvl_vaccinium/MxO_GeneticMap

VARDIR=$WD/results/variant_calling

# Path to directory containing the FASTQ files
FASTQDIR=$TEMPDIR/variant_calling/cleaned_fastq_files/


# Prefix of the indexed reference genomes
BLREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/Vmacrocarpon_BenLear_v1-scaffolds.fasta
STREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/V_macrocarpon_Stevens_v1.fasta
OXREFPREFIX=/project/gifvl_vaccinium/MxO_GeneticMap/results/genome_scaffolding_ragtag/Voxycoccos_NJ96-20_v1_ragtag_scaffolded/Voxycoccos_NJ96-20_v1_ragtag_scaffold.fasta


# Directory to output alignment files
# Putting the output in 90Daydata to save space on project directory
ALIGNDIR=$TEMPDIR/variant_calling/alignment/

# Directory to output coverage files
COVERAGEDIR=$VARDIR/genome_coverage/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################
 

## Run the pipeline

# Change working directory
cd $WD

# Create the coverage output directory if it does not exist
mkdir -p $COVERAGEDIR

# Create a function to run genome coverage calculation
run_genome_coverage() {
  local bamfile=$1
  local outputfile=${bamfile%.bam}_coverage.txt

  bedtools genomecov -bga -ibam $bamfile > $outputfile
}




# A list of reference genomes to align to
REFNAMES=("BenLear" "Stevens" "Oxycoccos")

# Iterate over each reference genome name
for REFNAME in ${REFNAMES[@]}; do
    echo -e "Processing reference genome: $REFNAME\n\n"

    # List the BAM files that end in "${REFNAME}_alignment.bam", but exclude the index files
    BAMFILES=($(find $ALIGNDIR -name "*${REFNAME}_alignment.bam" ! -name "*.bai"))

    # Run genome coverage calculation in parallel
    parallel -j $NTHREADS run_genome_coverage {} ::: "${BAMFILES[@]}"

done


# End of script