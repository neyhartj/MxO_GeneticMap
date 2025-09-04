#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=08:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=24   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="align stevens oxycoccos genomes"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# 
# Align the ragtag-rescaffolded oxycoccos genome to the Stevens reference and identify structural variants with SYRI
# Also produce a dotplot of the two genomes with mummerplot
# 

# Change the working directory
PROJ=/project/gifvl_vaccinium/MxO_GeneticMap

# Create output directories
# Main output directory first
RESULTS=${PROJ}/results
# Subdirectory for this analysis
SUBDIR=genome_alignment_syri

# Path to the reference fasta file
GENOME_DIR=/project/gifvl_vaccinium/vaccinium_genomes/
REF=${GENOME_DIR}/V_macrocarpon_Stevens_v1.fasta

# Path to the ragtag-rescaffolded oxycoccos fasta file
QUERY=$RESULTS/genome_scaffolding_ragtag/Voxycoccos_NJ96-20_v1_ragtag_scaffolded/Voxycoccos_NJ96-20_v1_ragtag_scaffold.fasta

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

# Define outdir
OUTDIR=$RESULTS/$SUBDIR

# Define a subdir for nucmer results
NUCMER_OUT=$OUTDIR/stevens_oxy_nucmer
mkdir -p $NUCMER_OUT
# Align the oxycoccos genome to the stevens using stevens as the reference
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT $REF $QUERY

# Filter for 1-1 alignments with at least 90% identity and at least 100 bp long
delta-filter -m -i 90 -l 100 $NUCMER_OUT.delta > ${NUCMER_OUT}_i90_l100.delta

# Produce a dotplot
mummerplot --png --large --layout -p ${NUCMER_OUT}_i90_l100_mumerplot ${NUCMER_OUT}_i90_l100.delta

# Convert to tab of coordinates
show-coords -THrd ${NUCMER_OUT}_i90_l100.delta > ${NUCMER_OUT}_i90_l100.coords

# Run SYRI
# Define a subdir for syri results
SYRI_OUT=$OUTDIR/stevens_oxy_syri
mkdir -p $SYRI_OUT
SYRI_PREFIX=${SYRI_OUT}_i90_l100_syri_

syri -r $REF -q $QUERY \
	-c ${NUCMER_OUT}_i90_l100.coords -d ${NUCMER_OUT}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG \
	-k \
	--nc 12 \
	--prefix $SYRI_PREFIX
# 


# # Visualize
# echo -e "#file\tname\ttags\n${newstevens2}\tstevens\n${oxyassembly}\toxycoccos" > ${OUTDIR}_genomes.txt
# plotsr --sr ${SYRI_PREFIX}.out --genomes ${OUTDIR}_genomes.txt -o ${SYRI_PREFIX}_plot.png



