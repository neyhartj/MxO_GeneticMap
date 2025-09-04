#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=168:00:00  # walltime limit (HH:MM:SS)
##SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=768G   # maximum memory per node
##SBATCH --mem=128G   # maximum memory per node
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --job-name="reassemble oxycoccos genome - hifiasm"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change the working directory
PROJ=/project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap
SUBDIR=$PROJ/analysis/GenomeAssembly

cd $SUBDIR

# Load modules
module load miniconda
module load minimap2
module load hifiasm
module load bbtools
module load bwa


# Load the environment
source activate wgs_env

NTHREADS=$SLURM_JOB_CPUS_PER_NODE



##############################################################################################


# Set input variable
FASTQ=$PROJ/data_raw/Oxycoccos_WGS_data/SRR14876345.fastq.gz


# Filter the ONT reads
OUTPUT=oxycoccos_ONT_SRR14876345_filtered.fastq
filtlong --min_length 1000 --keep_percent 90 $FASTQ > $OUTPUT

# Use nanoplot to check for quality
NanoPlot --fastq oxycoccos_ONT_SRR14876345_filtered.fastq -o oxycoccos_ONT_SRR14876345_filtered_nanoplot_out

INPUT=$OUTPUT

# Run preliminary assembly using hifiasm
# hifiasm -t $NTHREADS --ont -o oxycoccos_ONT.asm $INPUT

# Convert GFA to FASTA
# SAMPLE="oxycoccos_ONT"
# awk '/^S/{print ">"$2;print $3}' "$SAMPLE".asm.bp.p_ctg.gfa > "$SAMPLE".asm.bp.p_ctg.fasta

# Run stats
# stats.sh -Xmx256g in="$SAMPLE".asm.bp.p_ctg.fasta out="$SAMPLE".asm.bp.p_ctg.fasta.stats

# Test with FLYE
flye --nano-raw $INPUT -g 500m -o oxycoccos_ONT_flye_out -t $NTHREADS

# Run stats
stats.sh -Xmx256g in=oxycoccos_ONT_flye_out/assembly.fasta out=oxycoccos_ONT_flye_out/assembly_fasta.out

## Polish with pilon

# Index the assembly; if needed

