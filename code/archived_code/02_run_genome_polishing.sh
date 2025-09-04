#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

##SBATCH --time=168:00:00  # walltime limit (HH:MM:SS)
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=32G   # maximum memory per node
##SBATCH --mem=128G   # maximum memory per node
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --job-name="polish genome assembly"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change the working directory
PROJ=/project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap
SUBDIR=$PROJ/analysis/GenomeAssembly

cd $SUBDIR

# Load modules
module load bwa
module load trimmomatic


# Load the environment
source activate wgs_env

NTHREADS=$SLURM_JOB_CPUS_PER_NODE



##############################################################################################

# Set the input fastq
FORWARD=$PROJ/data_raw/Oxycoccos_WGS_data/SRR14876344_1.fastq.gz
REVERSE=$PROJ/data_raw/Oxycoccos_WGS_data/SRR14876344_2.fastq.gz

trimmomatic PE -threads $NTHREADS 
