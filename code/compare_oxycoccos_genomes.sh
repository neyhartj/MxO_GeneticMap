#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=08:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="compare oxycoccos genoems"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change the working directory
cd /project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap 

# Load modules
module load miniconda

# Load the environment
source activate wgs_env


# Set variables
genome1=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_oxycoccos_physical_v1.fasta
genome2=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_microcarpum_physical_v1.fasta


# Run the script
dnadiff $genome1 $genome2
