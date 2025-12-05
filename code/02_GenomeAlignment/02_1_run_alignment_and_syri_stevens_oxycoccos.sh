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
BEN_LEAR_REF=${GENOME_DIR}/Vmacrocarpon_BenLear_v1-scaffolds.fasta

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
# Simply add _renamed to the file name; use the filename variable above to create the new file name
newbenlear=${NEWGENOMEDIR}/$(basename ${BEN_LEAR_REF%".fasta"} )_renamed.fasta
# # Use bioawk to only keep chromosomes that start with "chr" and take the reverse complement of chromosomes 1, 4, 5, 6, and 9
# bioawk -c fastx '$name ~ /^chr/ { if ($name=="chr01" || $name=="chr04" || $name=="chr05") { print ">"$name"\n"revcomp($seq) } else { print ">"$name"\n"$seq } }' $BEN_LEAR_REF > $newbenlear
# The ben lear chromosomes are named "Vmac_chr[##]", so just remove Vmac_ prefix
bioawk -c fastx '$name ~ /^Vmac_chr/ { sub(/^Vmac_/, "", $name); printf ">%s\n%s\n", $name, $seq }' $BEN_LEAR_REF > $newbenlear

# Name of the modified stevens fasta
newstevens=${NEWGENOMEDIR}/$(basename ${STEVENS_REF%".fasta"} )_renamed.fasta
# Use bioawk to only keep the chromosomes that start with "chr", remove the suffix '_Vaccinium_macrocarpon_Stevens_v1', and rename them to "chr[##]"
# Do no modify the orientation of any chromosomes in Stevens
bioawk -c fastx '$name ~ /^chr/ { sub(/_.*/, "", $name); printf ">chr%02d\n%s\n", substr($name, 4), $seq }' $STEVENS_REF > $newstevens

# Name of the modified query fasta
oxyassembly=${NEWGENOMEDIR}/Voxycoccos_NJ96-20_v1_ragtag_scaffolded_renamed.fasta
bioawk -c fastx '$name ~ /^chr/ { sub(/_.*/, "", $name); printf ">chr%02d\n%s\n", substr($name, 4), $seq }' $OXY_REF > $oxyassembly

# Use bioawk to print the names and lengths of the chromosomes in each modified fasta file
echo -e "\nChromosome lengths in the modified Ben Lear reference fasta:"
bioawk -c fastx '{ print $name "\t" length($seq) }' $newbenlear
echo -e "\nChromosome lengths in the modified Stevens reference fasta:"
bioawk -c fastx '{ print $name "\t" length($seq) }' $newstevens
echo -e "\nChromosome lengths in the modified Oxycoccos query fasta:"
bioawk -c fastx '{ print $name "\t" length($seq) }' $oxyassembly
echo -e "\n"


## Run nucmer to align the genomes

# Define a subdir for nucmer results
NUCMER_OUT=$OUTDIR/nucmer
mkdir -p $NUCMER_OUT

# Define a list of genome pairs to align; iterate over this list
# The list should include the following:
# 1) Stevens vs Oxycoccos
# 2) Oxycoccos vs Stevens
# 3) Stevens vs Ben Lear
# 4) Ben Lear vs Stevens
# 5) Oxycoccos vs Ben Lear
# 6) Ben Lear vs Oxycoccos
ref_genome_list=("$newstevens" "$oxyassembly" "$newstevens" "$newbenlear" "$oxyassembly" "$newbenlear")
query_genome_list=("$oxyassembly" "$newstevens" "$newbenlear" "$newstevens" "$newbenlear" "$oxyassembly")
output_prefix_list=("stevens_oxy_nucmer" "oxy_stevens_nucmer" "stevens_benlear_nucmer" "benlear_stevens_nucmer" "oxy_benlear_nucmer" "benlear_oxy_nucmer")

# Iterate over the genome pairs and run nucmer
for i in ${!ref_genome_list[@]}; do
	ref_genome=${ref_genome_list[$i]}
	query_genome=${query_genome_list[$i]}
	output_prefix=${output_prefix_list[$i]}

	# Run nucmer
	nucmer --maxmatch --noextend -c 500 -b 500 -l 100 -p $NUCMER_OUT/$output_prefix $ref_genome $query_genome

done

# Filter the nucmer alignments and create coords files
# Filter for 1-1 alignments at least 10000 bp long
delta_files=($NUCMER_OUT/*.delta)
# Do not include "filter" delta files from previous runs
delta_files=(${delta_files[@]:0})
delta_files=(${delta_files[@]//*_filter.delta*})
# Iterate over delta files and filter them
for delta_file in "${delta_files[@]}"; do
	echo "Processing $delta_file"
	# Get the prefix (file name without path and extension)
	prefix=$(basename "$delta_file" .delta)
	# Filter for 1-1 alignments at least 10000 bp long
	delta-filter -1 -l 10000 "$delta_file" > "${NUCMER_OUT}/${prefix}_filter.delta"

	# Create a tab-delimited coords file
	show-coords -THrd "${NUCMER_OUT}/${prefix}_filter.delta" > "${NUCMER_OUT}/${prefix}_filter.coords"	
done

# Run SYRI for the nucmer alignments
# Define a subdir for syri results
SYRI_OUT=results/$SUBDIR/syri
mkdir -p $SYRI_OUT

# Iterate over the genome pairs again to run SYRI
# The list should include the following:
# 1) Stevens vs Oxycoccos
# 2) Oxycoccos vs Stevens
# 3) Stevens vs Ben Lear
# 4) Ben Lear vs Stevens
# 5) Oxycoccos vs Ben Lear
# 6) Ben Lear vs Oxycoccos

# Use the reference and query genome lists defined above
for i in ${!ref_genome_list[@]}; do
	ref_genome=${ref_genome_list[$i]}
	query_genome=${query_genome_list[$i]}
	output_prefix=${output_prefix_list[$i]}

	# Get the filtered delta and coords file names
	NUCMER_PREFIX=$NUCMER_OUT/$output_prefix
	# Define SYRI prefix
	SYRI_PREFIX=${output_prefix}_

	# Delta file
	DELTA_FILE=${NUCMER_PREFIX}_filter.delta
	# Coords file
	COORDS_FILE=${NUCMER_PREFIX}_filter.coords

	# Run SYRI
	syri -r $ref_genome -q $query_genome -c $COORDS_FILE -d $DELTA_FILE \
		--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX

	# Plot with plotsr
	GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
	# Get the base names of the reference and query genomes
	ref_name=$(basename $ref_genome)
	query_name=$(basename $query_genome)
	echo -e "#file\tname\ttags\n${ref_genome}\t${ref_name}\n${query_genome}\t${query_name}" > $GENOME_FILE
	SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
	plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

done




## Use minimap2 to for whole genome alignment (for comparison)
MINIMAP_OUT=$OUTDIR/minimap2
mkdir -p $MINIMAP_OUT

# Iterate over the genome pairs and run minimap2; use the list defined above
# create new output prefixes for minimap2
output_prefix_list=("stevens_oxy_minimap2" "oxy_stevens_minimap2" "stevens_benlear_minimap2" "benlear_stevens_minimap2" "oxy_benlear_minimap2" "benlear_oxy_minimap2")

# Iterate over the genome pairs again to run minimap2
for i in ${!ref_genome_list[@]}; do
	ref_genome=${ref_genome_list[$i]}
	query_genome=${query_genome_list[$i]}
	output_prefix=${output_prefix_list[$i]}

	# Align with minimap2
	minimap2 -ax asm5 -eqx -t $NTHREADS $ref_genome $query_genome | samtools sort -o $MINIMAP_OUT/${output_prefix}.bam
	samtools index $MINIMAP_OUT/${output_prefix}.bam

done

# Run SYRI on the minimap2 alignments
for i in ${!ref_genome_list[@]}; do
	ref_genome=${ref_genome_list[$i]}
	query_genome=${query_genome_list[$i]}
	output_prefix=${output_prefix_list[$i]}

	# BAM file from minimap2
	BAM_FILE=$MINIMAP_OUT/${output_prefix}.bam
	# Define SYRI prefix
	SYRI_PREFIX=${output_prefix}_

	# Run SYRI
	syri -r $ref_genome -q $query_genome -c $BAM_FILE \
		--lf ${SYRI_PREFIX}.log --log DEBUG -k --nc $NTHREADS --dir $SYRI_OUT --prefix $SYRI_PREFIX

	# Plot with plotsr
	GENOME_FILE=${SYRI_OUT}/${SYRI_PREFIX}genomes.txt
	# Get the base names of the reference and query genomes
	ref_name=$(basename $ref_genome)
	query_name=$(basename $query_genome)
	echo -e "#file\tname\ttags\n${ref_genome}\t${ref_name}\n${query_genome}\t${query_name}" > $GENOME_FILE
	SYRI_OUT_FILE=${SYRI_OUT}/${SYRI_PREFIX}syri.out
	plotsr --sr $SYRI_OUT_FILE --genomes $GENOME_FILE -o ${SYRI_OUT}/${SYRI_PREFIX}plot.png

done


# End of script