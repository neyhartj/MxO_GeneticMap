#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --job-name="MxO variant calling - variant filtration"
#SBATCH -p short
#SBATCH -t 01:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 2   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


## 
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 6. Filter variants from GATK genotype caller
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
WD=/project/gifvl_vaccinium/MxO_GeneticMap
VARDIR=$WD/results/variant_calling
VARIANTDIR=$VARDIR/variants/gatk/
GENOTYPEDIR=$VARIANTDIR/genotype_caller

# Filtering parameters for production SNPs
MAXMISSING=0.20 # Maximum missing data proportion
MINDP=10 # Minimum read depth
MINMAC=1 # Remove monomorphic SNPs
MINQ=40 # Minimum QUAL score


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# The max-missing input for VCFtools is 1 - missing proportion
VCFMAXMISSING=$(echo "1 - $MAXMISSING" | bc)

# List raw.vcf.gz files
VARIANTFILES=$(find $GENOTYPEDIR -name "*_alignment.raw.vcf.gz")

# Iterate over the files
for FILE in $VARIANTFILES; do

	# Run BCFtools stats
	OUTPUTSTAT=${FILE%".raw.vcf.gz"}_raw_stats.txt
	bcftools stats $FILE > $OUTPUTSTAT

	### Filter for production SNPS ##
	OUTPUT=${FILE%".raw.vcf.gz"}_filtered

	vcftools --gzvcf $FILE \
		--remove-indels \
		--min-alleles 2 \
		--max-alleles 2 \
		--max-missing $VCFMAXMISSING \
		--mac $MINMAC \
		--recode \
		--recode-INFO-all \
		--out $OUTPUT

	# Additional filter with bcftools for allele depth and quality
	VARIANTFILE1=${OUTPUT}.recode.vcf
	VARIANTFILE2=${OUTPUT}.vcf
	bcftools view -i "AD[GT] >= $MINDP && F_MISSING <= $MAXMISSING" $VARIANTFILE1 > $VARIANTFILE2

	# Remove the intermediate file
	rm $VARIANTFILE1

	# Zip
	gzip -f $VARIANTFILE2

	# Run stats on the output
	OUTPUTSTAT=${VARIANTFILE2%".vcf"}_stats.txt
	bcftools stats $VARIANTFILE2.gz > $OUTPUTSTAT

done
