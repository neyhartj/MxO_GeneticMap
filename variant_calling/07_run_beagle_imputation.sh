#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=40G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="beagle_imputation"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## Beagle imputation testing
## 
## This script will run Beagle on a set of VCF with random missing
## genotypes using various Beagle parameters. The output will be a set
## of imputed and phased VCF files for comparison with the original.
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load beagle-geno
module load vcftools

## Set variables
# Working directory
WD=/project/gifvl_vaccinium/cranberrygbs/cranberryHistoricalGBS/
# Name of input directory
INPUT=$WD/snps/
# Output directory
OUTPUT=$WD/imputation/beagle_imputation

# Name of input vcf
VCFIN=$INPUT/cranberryGBS_production_snps_resolvedDuplicates_filtered.vcf.gz


# Change working directory
cd $WD 


## Filter

# Filter to retain the germplasm collection, then filter SNPs for excessive 
# missingness (> 50%) and then filter individuals for excessive missingness (> 50%)
vcftools --gzvcf $VCFIN --keep $WD/input/cranberry_gbs_germplasm_collection_individuals.txt \
--recode-INFO-all --recode --stdout | \
vcftools --vcf - --mac 15 --max-missing 0.5 --recode-INFO-all --recode --out $OUTPUT/cranberryGBS_snps_to_impute

# Calculate individual missingness
vcftools --vcf $OUTPUT/cranberryGBS_snps_to_impute.recode.vcf --missing-indv \
--out $OUTPUT/cranberryGBS_germplasm_missingness

# Filter out individuals with missingness > 0.50
awk 'NR > 1 { if ($5 > 0.50) print $1 }' $OUTPUT/cranberryGBS_germplasm_missingness.imiss > $OUTPUT/individual_remove_imputation.txt

vcftools --vcf $OUTPUT/cranberryGBS_snps_to_impute.recode.vcf --remove $OUTPUT/individual_remove_imputation.txt \
--recode-INFO-all --recode --out $OUTPUT/cranberryGBS_snps_to_impute_filtered


## Run beagle

# Run the imputation software
java -Xmx32g -jar /software/7/apps/beagle-geno/5.0/beagle.16May19.351.jar \
gt=$OUTPUT/cranberryGBS_snps_to_impute_filtered.recode.vcf \
out=$OUTPUT/cranberryGBS_germplasm_imputed_snps \
ne=1000 burnin=6 iterations=25

