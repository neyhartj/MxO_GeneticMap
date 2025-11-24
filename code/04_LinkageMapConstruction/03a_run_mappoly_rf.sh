#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="MxO RAPiD variant calling pipeline - Fast QC"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


##
## MxO Genetic Map
##
## Variant calling pipeline
##
## Step 1. Quality control of raw reads
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load r

#####################
## Set variables
#####################

# Working directory
WD=/project/gifvl_vaccinium/MxO_GeneticMap

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##########################################
## Run the pipeline - DO NOT EDIT BELOW
##########################################

## Step 1: Quality control of raw reads

# Change working directory
cd $WD

# Run the R script to calculate recombination fractions
Rscript code/04_LinkageMapConstruction/03a_run_mappoly_rf.R