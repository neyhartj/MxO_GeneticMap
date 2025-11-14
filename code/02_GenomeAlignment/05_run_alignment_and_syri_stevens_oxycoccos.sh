#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=04:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12   # X processor core(s) per node X 2 threads per core
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
SUBDIR=genome_alignment

# Path to the directory containing genomes
GENOME_DIR=/project/gifvl_vaccinium/vaccinium_genomes/

# Path to the Ben Lear reference fasta file.
BEN_LEAR_REF=${GENOME_DIR}/Vaccinium_macrocarpon_BenLear_v2.fasta

# Path to the Stevens fasta file
STEVENS_REF=${GENOME_DIR}/V_macrocarpon_Stevens_v1.fasta

# Path to the ragtag-rescaffolded oxycoccos fasta file
OXY_REF=$RESULTS/genome_scaffolding_ragtag/Voxycoccos_NJ96-20_v1_ragtag_scaffolded/Voxycoccos_NJ96-20_v1_ragtag_scaffold.fasta

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


# Modify all genomes to make sure the chromosome names are identical
# This is important for SYRI to work properly
# The REF fasta chromosome name follow the pattern "chr[##]_Vaccinium_macrocarpon_Stevens_v1". We will rename them to "chr[##]"
# The QUERY fasta chromosome names follow the pattern "chr[##]_Vaccinium_macrocarpon_Stevens_v1_RagTag". We will rename them to "chr[##]"
# For both genomes, we will only keep the chromosomes that start with "chr" (i.e. we will remove any contigs that do not belong to a chromosome)
# Create modified versions of the fasta files in the outdir
NEWGENOMEDIR=$RESULTS/modified_genomes
mkdir -p $NEWGENOMEDIR

# Name of the modified ben lear fasta
newbenlear=${NEWGENOMEDIR}/Vaccinium_macrocarpon_BenLear_v2_renamed.fasta
# Use bioawk to only keep chromosomes that start with "chr" and take the reverse complement of chromosomes 1, 4, 5, 6, and 9
bioawk -c fastx '$name ~ /^chr/ { if ($name=="chr01" || $name=="chr04" || $name=="chr05") { print ">"$name"\n"revcomp($seq) } else { print ">"$name"\n"$seq } }' $BEN_LEAR_REF > $newbenlear

# Name of the modified stevens fasta
newstevens=${NEWGENOMEDIR}/V_macrocarpon_Stevens_v1_renamed.fasta
# Use bioawk to only keep the chromosomes that start with "chr", remove the suffix '_Vaccinium_macrocarpon_Stevens_v1', and rename them to "chr[##]"
# Do no modify the orientation of any chromosomes in Stevens
bioawk -c fastx '$name ~ /^chr/ { sub(/_.*/, "", $name); printf ">chr%02d\n%s\n", substr($name, 4), $seq }' $STEVENS_REF > $newstevens

# Name of the modified query fasta
oxyassembly=${NEWGENOMEDIR}/Voxycoccos_NJ96-20_v1_ragtag_scaffolded_renamed.fasta
bioawk -c fastx '$name ~ /^chr/ { sub(/_.*/, "", $name); printf ">chr%02d\n%s\n", substr($name, 4), $seq }' $OXY_REF > $oxyassembly

# Use bioawk to print the names and lengths of the chromosomes in each modified fasta file
echo -e "\nChromosome lengths in the modified Ben Lear reference fasta:"
bioawk -c fastx '$name ~ /^chr/ { print $name "\t" length($seq) }' $newbenlear
echo -e "\nChromosome lengths in the modified Stevens reference fasta:"
bioawk -c fastx '$name ~ /^chr/ { print $name "\t" length($seq) }' $newstevens
echo -e "\nChromosome lengths in the modified Oxycoccos query fasta:"
bioawk -c fastx '$name ~ /^chr/ { print $name "\t" length($seq) }' $oxyassembly
echo -e "\n"

# Define a subdir for nucmer results
NUCMER_OUT=$OUTDIR/nucmer
mkdir -p $NUCMER_OUT

# Align Oxy to Stevens with nucmer
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/stevens_oxy_nucmer $newstevens $oxyassembly
# Align Stevens to Oxy with nucmer (reverse direction)
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/oxy_stevens_nucmer $oxyassembly $newstevens

# Align Stevens to Ben Lear with nucmer
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/stevens_benlear_nucmer $newstevens $newbenlear
# Align Ben Lear to Stevens with nucmer (reverse direction)
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/benlear_stevens_nucmer $newbenlear $newstevens	

# Align Oxy to Ben Lear with nucmer
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/oxy_benlear_nucmer $oxyassembly $newbenlear
# Align Ben Lear to Oxy with nucmer (reverse direction)
nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/benlear_oxy_nucmer $newbenlear $oxyassembly

# Find all delta files
delta_files=($NUCMER_OUT/*.delta)
# Iterate over delta files and filter them
for delta_file in "${delta_files[@]}"; do
	echo "Processing $delta_file"
	# Get the prefix (file name without path and extension)
	prefix=$(basename "$delta_file" .delta)
	# Filter for 1-1 alignments at least 1000 bp long
	delta-filter -1 -l 1000 "$delta_file" > "${NUCMER_OUT}/${prefix}_i90_l100.delta"

	# Create a tab-delimited coords file
	show-coords -THrd "${NUCMER_OUT}/${prefix}_i90_l100.delta" > "${NUCMER_OUT}/${prefix}_i90_l100.coords"	
done

# Run SYRI
# Define a subdir for syri results
SYRI_OUT=results/$SUBDIR/syri
mkdir -p $SYRI_OUT

# Run SYRI for the Stevens vs Oxycoccos alignment
SYRI_PREFIX=stevens_oxy_
NUCMER_PREFIX=$NUCMER_OUT/stevens_oxy_nucmer
syri -r $newstevens -q $oxyassembly -c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX

# Plot with plotsr
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${newstevens}\tstevens\n${oxyassembly}\toxycoccos" > $GENOME_FILE
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png


# Run SYRI for the Oxycoccos vs Stevens alignment
SYRI_PREFIX=oxy_stevens_
NUCMER_PREFIX=$NUCMER_OUT/oxy_stevens_nucmer
syri -r $oxyassembly -q $newstevens -c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX	

# Plot with plotsr
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${oxyassembly}\toxycoccos\n${newstevens}\tstevens" > $GENOME_FILE
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

# Run SYRI for the Stevens vs Ben Lear alignment
SYRI_PREFIX=stevens_benlear_
NUCMER_PREFIX=$NUCMER_OUT/stevens_benlear_nucmer
syri -r $newstevens -q $newbenlear -c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX

# Plot with plotsr
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${newstevens}\tstevens\n${newbenlear}\tbenlear" > $GENOME_FILE
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

# Run SYRI for the Ben Lear vs Stevens alignment
SYRI_PREFIX=benlear_stevens_
NUCMER_PREFIX=$NUCMER_OUT/benlear_stevens_nucmer
syri -r $newbenlear -q $newstevens -c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX

# Plot with plotsr
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${newbenlear}\tbenlear\n${newstevens}\tstevens" > $GENOME_FILE
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

# Run SYRI for the Oxycoccos vs Ben Lear alignment
SYRI_PREFIX=oxy_benlear_
NUCMER_PREFIX=$NUCMER_OUT/oxy_benlear_nucmer
syri -r $oxyassembly -q $newbenlear -c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX	

# Plot with plotsr
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${oxyassembly}\toxycoccos\n${newbenlear}\tbenlear" > $GENOME_FILE
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

# Run SYRI for the Ben Lear vs Oxycoccos alignment
SYRI_PREFIX=benlear_oxy_
NUCMER_PREFIX=$NUCMER_OUT/benlear_oxy_nucmer
syri -r $newbenlear -q $oxyassembly -c ${NUCMER_PREFIX}_i90_l100.coords -d ${NUCMER_PREFIX}_i90_l100.delta \
	--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX	

# Plot with plotsr
GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
echo -e "#file\tname\ttags\n${newbenlear}\tbenlear\n${oxyassembly}\toxycoccos" > $GENOME_FILE
SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png


# End of script