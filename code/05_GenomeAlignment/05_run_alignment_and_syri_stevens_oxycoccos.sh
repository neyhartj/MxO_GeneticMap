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
NUCMER_OUT=$OUTDIR/nucmer
mkdir -p $NUCMER_OUT
NUCMER_PREFIX=${NUCMER_OUT}/stevens_oxy_nucmer
# Align the oxycoccos genome to the stevens using stevens as the reference
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_PREFIX $REF $QUERY

# Filter for 1-1 alignments with at least 90% identity and at least 100 bp long
delta-filter -m -i 90 -l 100 $NUCMER_PREFIX.delta > ${NUCMER_PREFIX}_i90_l100.delta

# Produce a dotplot
MUMMERPLOT_OUT=$OUTDIR/mummerplot
mkdir -p $MUMMERPLOT_OUT
MUMMERPLOT_PREFIX=${MUMMERPLOT_OUT}/stevens_oxy_mummerplot
mummerplot --png --large --layout -p ${MUMMERPLOT_PREFIX} ${NUCMER_PREFIX}_i90_l100.delta

# Convert to tab of coordinates
show-coords -THrd ${NUCMER_PREFIX}_i90_l100.delta > ${NUCMER_PREFIX}_i90_l100.coords

# Run SYRI
# Define a subdir for syri results
SYRI_OUT=$OUTDIR/syri
mkdir -p $SYRI_OUT
SYRI_PREFIX=stevens_oxy_syri

syri -r $REF -q $QUERY \
	-c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG \
	-k \
	--nc 12 \
	--dir $SYRI_OUT \
	--prefix $SYRI_PREFIX
# 


# # Visualize
# echo -e "#file\tname\ttags\n${newstevens2}\tstevens\n${oxyassembly}\toxycoccos" > ${OUTDIR}_genomes.txt
# plotsr --sr ${SYRI_PREFIX}.out --genomes ${OUTDIR}_genomes.txt -o ${SYRI_PREFIX}_plot.png



