#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="MxO RAPiD freebayes variant calling"
#SBATCH -p short
#SBATCH -t 48:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 72   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=216G   # maximum memory per node
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
DBPREFIXOXY=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_oxycoccos_physical_v1.fasta

# Directory with input alignment files
ALIGNDIR=$WD/variant_calling/merged_alignment/

# Directory to output variants
VARIANTDIR=$WD/variant_calling/variants/

# Path to BED file containing the RAPiD Flex-Seq probes
PROBEBED=/project/gifvl_vaccinium/cranberryGenotyping/RAPiD_Cranberry_15K/vm_flexseq_probes.bed

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE

# Run Stevens alignment or Oxy?
REF=STEVENS
# REF=OXY


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# List all of the alignment files in the alignment directory
# First find those aligned to STEVENS
ALIGNMENTFILESSTE=$(find $ALIGNDIR -name "*stevens*.bam")
# Next find those aligned to OXY
ALIGNMENTFILESOXY=$(find $ALIGNDIR -name "*oxy*.bam")


if [ $REF = "STEVENS" ];
then

	echo -e "Variant calling using the STEVENS alignment...\n"

	## Run variant calling for alignment to STEVENS
	# Ouput file
	OUTPUT=$VARIANTDIR/mxo_stevens-ref_variants.vcf
	# freebayes -f $DBPREFIXSTEVENS $ALIGNMENTFILESSTE > $OUTPUT
	# Use parallelization
	freebayes-parallel <(fasta_generate_regions.py $DBPREFIXSTEVENS 100000) $SLURM_JOB_CPUS_PER_NODE \
		-f $DBPREFIXSTEVENS --min-alternate-count 10 --use-best-n-alleles 2 -t $PROBEBED $ALIGNMENTFILESSTE > $OUTPUT

elif [ $REF = "OXY" ];
then

	echo -e "Variant calling using the OXY alignment...\n"

	## Run variant calling for alignment to OXY
	# Ouput file
	OUTPUT=$VARIANTDIR/mxo_oxy-ref_variants.vcf
	# freebayes -f $DBPREFIXOXY $ALIGNMENTFILESOXY > $OUTPUT

	# Use parallelization
	freebayes-parallel <(fasta_generate_regions.py $DBPREFIXOXY 100000) $SLURM_JOB_CPUS_PER_NODE \
	  -f $DBPREFIXOXY --min-alternate-count 10 --use-best-n-alleles 2 $ALIGNMENTFILESOXY > $OUTPUT

fi
