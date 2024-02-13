#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="MxO RAPiD freebayes variant filtration"
#SBATCH -p short
#SBATCH -t 02:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 4. Filter variants from freebayes
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load vcftools
module load bcftools

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

# Directory to output alignments
ALIGNDIR=$WD/variant_calling/alignment/

# Directory to output variants
VARIANTDIR=$WD/variant_calling/variants/

# Path to BED file containing the RAPiD Flex-Seq probes
PROBEBED=/project/gifvl_vaccinium/cranberryGenotyping/RAPiD_Cranberry_15K/vm_flexseq_probes.bed

# Parameters for VCFtools filtration
MAXMISSING=0.2
MINDP=10



##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Get the variant files
VARIANTFILESTE=$VARIANTDIR/mxo_stevens-ref_variants.vcf.gz
VARIANTFILEOXY=$VARIANTDIR/mxo_oxy-ref_variants.vcf.gz

#######################
## STEVENS alignment
#######################

# bcftools stats
OUTPUTSTAT=${VARIANTFILESTE%".vcf.gz"}_stats.txt
bcftools stats $VARIANTFILESTE > $OUTPUTSTAT


OUTPUT=${VARIANTFILESTE%".vcf.gz"}_filtered.vcf.gz

vcftools --gzvcf $VARIANTFILESTE \
	--bed $PROBEBED \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--max-missing $MAXMISSING \
	--minDP $MINDP \
	--recode \
	--recode-INFO-all \
	--stdout | gzip -c > $OUTPUT


# bcftools stats
OUTPUTSTAT=${OUTPUT%".vcf.gz"}_stats.txt
bcftools stats $OUTPUT > $OUTPUTSTAT





#######################
## OXY alignment
#######################

# bcftools stats
OUTPUTSTAT=${VARIANTFILEOXY%".vcf.gz"}_stats.txt
bcftools stats $VARIANTFILEOXY > $OUTPUTSTAT

OUTPUT=${VARIANTFILEOXY%".vcf.gz"}_filtered.vcf.gz

vcftools --gzvcf $VARIANTFILEOXY \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--max-missing $MAXMISSING \
	--minDP $MINDP \
	--recode \
	--recode-INFO-all \
	--stdout | gzip -c > $OUTPUT


# bcftools stats
OUTPUTSTAT=${OUTPUT%".vcf.gz"}_stats.txt
bcftools stats $OUTPUT > $OUTPUTSTAT


