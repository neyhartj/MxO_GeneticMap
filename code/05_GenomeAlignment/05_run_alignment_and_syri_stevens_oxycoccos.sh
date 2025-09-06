#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=02:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=8G   # maximum memory per node
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


# Modify the REF and QUERY fasta files to make sure the chromosome names are identical
# This is important for SYRI to work properly
# The REF fasta chromosome name follow the pattern "chr[##]_Vaccinium_macrocarpon_Stevens_v1". We will rename them to "chr[##]"
# The QUERY fasta chromosome names follow the pattern "chr[##]_Vaccinium_macrocarpon_Stevens_v1_RagTag". We will rename them to "chr[##]"
# For both genomes, we will only keep the chromosomes that start with "chr" (i.e. we will remove any contigs that do not belong to a chromosome)
# Create modified versions of the fasta files in the outdir
NEWGENOMEDIR=$RESULTS/modified_genomes
mkdir -p $NEWGENOMEDIR

# Name of the modified reference fasta
newstevens=${NEWGENOMEDIR}/V_macrocarpon_Stevens_v1_renamed.fasta
# Rename the chromosomes in the reference fasta and only keep those that start with "chr"
awk '
/^>/ {
	if ($0 ~ /^>chr[0-9]+_Vaccinium_macrocarpon_Stevens_v1/) {
		match($0, /^>chr([0-9]+)_Vaccinium_macrocarpon_Stevens_v1/, arr)
		num = arr[1]
		if (length(num) == 1) num = "0" num
		print ">chr" num
		keep=1
	} else {
		keep=0
	}
	next
}
keep { print }
' $REF > $newstevens

# Name of the modified query fasta
oxyassembly=${NEWGENOMEDIR}/Voxycoccos_NJ96-20_v1_ragtag_scaffolded_renamed.fasta
# Rename the chromosomes in the query fasta
awk '
/^>/ {
	if ($0 ~ /^>chr[0-9]+_Vaccinium_macrocarpon_Stevens_v1_RagTag/) {
		match($0, /^>chr([0-9]+)_Vaccinium_macrocarpon_Stevens_v1_RagTag/, arr)
		num = arr[1]
		if (length(num) == 1) num = "0" num
		print ">chr" num
		keep=1
	} else {
		keep=0
	}
	next
}
keep { print }
' $QUERY > $oxyassembly

# Define a subdir for nucmer results
NUCMER_OUT=$OUTDIR/nucmer
mkdir -p $NUCMER_OUT
NUCMER_PREFIX=${NUCMER_OUT}/stevens_oxy_nucmer
# Align the oxycoccos genome to the stevens using stevens as the reference; use the renamed versions of the fasta files
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_PREFIX $newstevens $oxyassembly

# Filter for 1-1 alignments with at least 90% identity and at least 100 bp long
delta-filter -1 -l 1000 $NUCMER_PREFIX.delta > ${NUCMER_PREFIX}_i90_l100.delta

# # Produce a dotplot
# MUMMERPLOT_OUT=$OUTDIR/mummerplot
# mkdir -p $MUMMERPLOT_OUT
# MUMMERPLOT_PREFIX=${MUMMERPLOT_OUT}/stevens_oxy_mummerplot
# mummerplot --png --large --layout -p ${MUMMERPLOT_PREFIX} ${NUCMER_PREFIX}_i90_l100.delta

# Convert to tab of coordinates
show-coords -THrd ${NUCMER_PREFIX}_i90_l100.delta > ${NUCMER_PREFIX}_i90_l100.coords

# Run SYRI
# Define a subdir for syri results
SYRI_OUT=results/syri
mkdir -p $SYRI_OUT
SYRI_PREFIX=stevens_oxy_

# Run SYRI using the filtered nucmer alignments and the renamed fasta files
syri -r $newstevens \
	-q $oxyassembly \
	-c ${NUCMER_PREFIX}_i90_l100.coords \
	-d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log \
	--log DEBUG \
	-k \
	--nc $NTHREADS \
	--dir $SYRI_OUT \
	--prefix $SYRI_PREFIX
# 


# Visualize
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${newstevens}\tstevens\n${oxyassembly}\toxycoccos" > $GENOME_FILE

# Output of syri
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out


plotsr --sr $SYRI_OUT_FILE \
	--genomes $GENOME_FILE \
	-o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

# End of script

