#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="MxO alignment duplicate filtration and merging"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=156G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 2. Mark and remove duplicates in the BAM files and then merge bam files
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load samtools
module load picard

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

# Directory containing alignment files
ALIGNDIR=$WD/variant_calling/alignment/

# Directory to output merged alignment files
MERGEDALIGNDIR=$WD/variant_calling/merged_alignment

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

# Print a message
echo -e "Removing duplicates and merging the $REF alignments...\n"

if [ $REF = "STEVENS" ];
then

    ## STEVEN alignment files

    # List files for the Stevens alignment
    BAMFILES=$(find $ALIGNDIR -name "*STEVENS_alignment.bam")

    # Iterate over the alignment files
    for BAMFILE in $BAMFILES; do

    	# Output file from duplicate marking
    	OUTPUT=${BAMFILE%".bam"}_nodup.bam
    	MARKDUPOUT=${BAMFILE%".bam"}_duplicate_metrics.txt

    	# Mark and remove duplicates
    	java -Xmx100G -jar /software/el9/apps/picard/3.0.0/picard.jar MarkDuplicates \
    		--REMOVE_DUPLICATES true \
    		-I $BAMFILE \
    		-O $OUTPUT \
    		-M $MARKDUPOUT

    done


    # Collect new BAM file names
    NEWBAMFILES=$(find $ALIGNDIR -name "*STEVENS_alignment_nodup.bam")

    # Merge the bam files
    # Sort on coordinates
    # Subset the bam files for only those positions overlapping with the probe BED file
    samtools merge -@ $SLURM_JOB_CPUS_PER_NODE -o - $NEWBAMFILES | \
    	samtools sort -@ $SLURM_JOB_CPUS_PER_NODE -u - | \
    	samtools view -b -o $MERGEDALIGNDIR/mxo_stevens_alignments_merged.bam -@ $SLURM_JOB_CPUS_PER_NODE -L $PROBEBED -

    # Index
    samtools index $MERGEDALIGNDIR/mxo_stevens_alignments_merged.bam

elif [ $REF = "OXY" ];
then

    ## OXY alignment files

    # List files for the Stevens alignment
    BAMFILES=$(find $ALIGNDIR -name "*OXY_alignment.bam")

    # Iterate over the alignment files
    for BAMFILE in $BAMFILES; do

    	# Output file from duplicate marking
    	OUTPUT=${BAMFILE%".bam"}_nodup.bam
    	MARKDUPOUT=${BAMFILE%".bam"}_duplicate_metrics.txt

    	# Mark and remove duplicates
    	java -Xmx100G -jar /software/el9/apps/picard/3.0.0/picard.jar MarkDuplicates \
    		--REMOVE_DUPLICATES true \
    		-I $BAMFILE \
    		-O $OUTPUT \
    		-M $MARKDUPOUT

    done


    # Collect new BAM file names
    NEWBAMFILES=$(find $ALIGNDIR -name "*OXY_alignment_nodup.bam")

    # Merge the bam files
    # Sort on coordinates
    # Subset the bam files for only those positions overlapping with the probe BED file
    samtools merge -@ $SLURM_JOB_CPUS_PER_NODE -o - $NEWBAMFILES | \
    	samtools sort -@ $SLURM_JOB_CPUS_PER_NODE -u - | \
    	samtools view -b -o $MERGEDALIGNDIR/mxo_oxy_alignments_merged.bam -@ $SLURM_JOB_CPUS_PER_NODE -

    # Index
    samtools index $MERGEDALIGNDIR/mxo_oxy_alignments_merged.bam

fi

