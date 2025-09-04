#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --job-name="scaffold the oxycoccos genome with ragtag"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# 
# Use RagTag to scaffold the oxycoccos genome assembly using the Vaccinium macrocarpon genome as a reference
# 

# Change the working directory
PROJ=/project/gifvl_vaccinium/MxO_GeneticMap

# Create output directories
# Main output directory first
RESULTS=${PROJ}/results
# Subdirectory for this analysis
SUBDIR=genome_scaffolding_ragtag

# Path to the reference fasta file
GENOME_DIR=/project/gifvl_vaccinium/vaccinium_genomes/
REF=${GENOME_DIR}/V_macrocarpon_Stevens_v1.fasta

# Path to the query fasta file
# This is the assembly to be scaffolded
QUERY=${GENOME_DIR}/Voxycoccos_NJ96-20_v1-scaffolds.fasta

# Name of the output of the RagTag scaffolding
ASSEMBLYNAME=Voxycoccos_NJ96-20_v1_ragtag_scaffolded

NTHREADS=$SLURM_JOB_CPUS_PER_NODE

##############################################################################################
## Generally do not edit below this line
##############################################################################################


# Load modules
module load miniconda
module load minimap2

# Load the environment
source activate wgs_env


# Change to the project directory
cd $PROJ

mkdir -p $RESULTS/$SUBDIR

# Output path
OUTPUT=$RESULTS/$SUBDIR/$ASSEMBLYNAME

## Scaffold the assembly with RagTag
ragtag.py scaffold -o $OUTPUT -t $NTHREADS $REF $QUERY

# End of script