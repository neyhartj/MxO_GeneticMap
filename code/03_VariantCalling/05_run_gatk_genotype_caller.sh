#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --job-name="MxO variant calling - genotype caller"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


## 
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 4. Variant calling with GATK genotype caller
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load samtools
module load gatk

#####################
## Set variables
#####################

# Working directory
WD=/project/gifvl_vaccinium/MxO_GeneticMap
# Path to variant calling results directory
VARDIR=$WD/results/variant_calling
# Directory of gatk output
VARIANTDIR=$VARDIR/variants/gatk/

# Path to temporary directory on 90Daydata
TEMPDIR=/90daydata/gifvl_vaccinium/MxO_GeneticMap
# Path to .vcf files from haplotype caller on 90Daydata
HAPLODIR=$TEMPDIR/variant_calling/variants/gatk/haplotype_caller

# Prefix of the indexed reference genomes
BLREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/Vaccinium_macrocarpon_BenLear_v2.fasta
STREFPREFIX=/project/gifvl_vaccinium/vaccinium_genomes/V_macrocarpon_Stevens_v1.fasta
OXREFPREFIX=//project/gifvl_vaccinium/MxO_GeneticMap/results/genome_scaffolding_ragtag/Voxycoccos_NJ96-20_v1_ragtag_scaffolded/Voxycoccos_NJ96-20_v1_ragtag_scaffold.fasta

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Make a directory for the genotype caller output
GENOTYPEDIR=$VARIANTDIR/genotype_caller
if [ ! -d $GENOTYPEDIR ]; then
  mkdir $GENOTYPEDIR
fi

  
## Run for each reference genome ##

# Find the GVCF files from the BenLear alignment
BENLEARGVCFFILES=$(for file in $(find $HAPLODIR -name "*BenLear*.g.vcf.gz"); do echo -n "-V $file "; done)
# Combine the GVCF files
BENLEARGVCF=$GENOTYPEDIR/mxo_variant_cohort_BenLear_alignment.g.vcf

gatk --java-options "-Xmx8g" CombineGVCFs -R $BLREFPREFIX $BENLEARGVCFFILES -O $BENLEARGVCF
# Genotype the combined GVCF file
gatk --java-options "-Xmx8g" GenotypeGVCFs -R $BLREFPREFIX -V $BENLEARGVCF -O $GENOTYPEDIR/mxo_variant_cohort_BenLear_alignment.raw.vcf.gz

wait

# Run for the Stevens reference genome

# Find the GVCF files from the Stevens alignment
STEVENSGVCFFILES=($(find $HAPLODIR -name "*Stevens*.g.vcf.gz"))
# Combine the GVCF files
STEVENSGVCF=$GENOTYPEDIR/mxo_variant_cohort_Stevens_alignment.g.vcf
gatk --java-options "-Xmx8g" CombineGVCFs -R $STREFPREFIX $(for file in ${STEVENSGVCFFILES[@]}; do echo -n "-V $file "; done) -O $STEVENSGVCF
# Genotype the combined GVCF file
gatk --java-options "-Xmx8g" GenotypeGVCFs -R $STREFPREFIX -V $STEVENSGVCF -O $GENOTYPEDIR/mxo_variant_cohort_Stevens_alignment.raw.vcf.gz

wait

# Run for the Oxycoccos reference genome
# Find the GVCF files from the Oxycoccos alignment
OXYCOCCOSGVCFFILES=($(find $HAPLODIR -name "*Oxycoccos*.g.vcf.gz"))
# Combine the GVCF files
OXYCOCCOSGVCF=$GENOTYPEDIR/mxo_variant_cohort_Oxycoccos_alignment.g.vcf
gatk --java-options "-Xmx8g" CombineGVCFs -R $OXREFPREFIX $(for file in ${OXYCOCCOSGVCFFILES[@]}; do echo -n "-V $file "; done) -O $OXYCOCCOSGVCF
# Genotype the combined GVCF file
gatk --java-options "-Xmx8g"  GenotypeGVCFs -R $OXREFPREFIX -V $OXYCOCCOSGVCF -O $GENOTYPEDIR/mxo_variant_cohort_Oxycoccos_alignment.raw.vcf.gz

wait

# End of script



