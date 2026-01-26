#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gifvl_vaccinium
#SBATCH --time=08:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="align stevens oxycoccos genomes - diversity"
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
NEWGENOMEDIR=$RESULTS/modified_genomes

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
module load bcftools

# Load the environment
source activate wgs_env


# Change to the project directory
cd $PROJ

mkdir -p $RESULTS/$SUBDIR

# Define outdir
OUTDIR=$RESULTS/$SUBDIR

# Directory with output from minimap2
MINIMAP_OUT=$OUTDIR/minimap2


# Name of the modified ben lear fasta
newbenlear=${NEWGENOMEDIR}/$(basename ${BEN_LEAR_REF%".fasta"} )_renamed.fasta
# Name of the modified stevens fasta
newstevens=${NEWGENOMEDIR}/$(basename ${STEVENS_REF%".fasta"} )_renamed.fasta
# Name of the modified query fasta
oxyassembly=${NEWGENOMEDIR}/Voxycoccos_NJ96-20_v1_ragtag_scaffolded_renamed.fasta


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
output_prefix_list=("stevens_oxy_minimap2" "oxy_stevens_minimap2" "stevens_benlear_minimap2" "benlear_stevens_minimap2" "oxy_benlear_minimap2" "benlear_oxy_minimap2")


# Use the reference and query genome lists defined above
for i in ${!ref_genome_list[@]}; do
	ref_genome=${ref_genome_list[$i]}
	query_genome=${query_genome_list[$i]}
	output_prefix=${output_prefix_list[$i]}

	# Print message
	printf "\n\nRunning mpileup and pixy for minimpa output of %s vs %s\n" "$(basename $ref_genome)" "$(basename $query_genome)"

	# Get the .bam file with this output prefix
	BAM_FILE=$MINIMAP_OUT/${output_prefix}.bam

	# Generate mpileup using bcftools and output to a compressed vcf file; use multi-threading
	bcftools mpileup -f $ref_genome $BAM_FILE --threads $SLURM_JOB_CPUS_PER_NODE | \
	bcftools call -m -Oz -o $MINIMAP_OUT/temp_variants.vcf.gz --threads $SLURM_JOB_CPUS_PER_NODE

	# Index the vcf file
	bcftools index $MINIMAP_OUT/temp_variants.vcf.gz

	# Create a txt file with the sample name (the bam file name) to use for reheadering
	echo -e "query" > $MINIMAP_OUT/new_sample_names.txt

	# Use bcftools reheader to change the sample name to the bam file name
	bcftools reheader -s $MINIMAP_OUT/new_sample_names.txt $MINIMAP_OUT/temp_variants.vcf.gz > $MINIMAP_OUT/temp_variants_renamed.vcf.gz
	bcftools index $MINIMAP_OUT/temp_variants_renamed.vcf.gz


	# Create a dummy VCF for the reference genome to include in pixy analysis
	# This uses the 'header' from your query VCF but sets genotypes to 0/0
	bcftools view -h $MINIMAP_OUT/temp_variants.vcf.gz > $MINIMAP_OUT/header.txt

	# Create a version of the VCF where every site is 0/0 (matching the reference)
	bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t.\t.\tGT\t0/0\n' $MINIMAP_OUT/temp_variants.vcf.gz | \
	cat $MINIMAP_OUT/header.txt - | bcftools view -Oz -o $MINIMAP_OUT/dummy.vcf.gz

	# New sample name file
	echo -e "reference" > $MINIMAP_OUT/new_sample_names.txt
	# Reheader the dummy VCF to have sample name 'reference'
	bcftools reheader -s $MINIMAP_OUT/new_sample_names.txt $MINIMAP_OUT/dummy.vcf.gz > $MINIMAP_OUT/dummy_renamed.vcf.gz
	# Index
	bcftools index $MINIMAP_OUT/dummy_renamed.vcf.gz

	# Combine the real VCF and the dummy reference VCF
	bcftools merge $MINIMAP_OUT/dummy_renamed.vcf.gz $MINIMAP_OUT/temp_variants_renamed.vcf.gz -Oz -o $MINIMAP_OUT/${output_prefix}_variants.vcf.gz
	bcftools index $MINIMAP_OUT/${output_prefix}_variants.vcf.gz

	# Remove dummy and temp
	rm $MINIMAP_OUT/dummy.vcf.gz* $MINIMAP_OUT/dummy_renamed.vcf.gz* $MINIMAP_OUT/temp_variants.vcf.gz* $MINIMAP_OUT/temp_variants_renamed.vcf.gz* $MINIMAP_OUT/header.txt
	# Create a populations file for pixy
	# The VCF will contain a single sample with the name output_prefix.bam
	echo -e "query\tp1\nreference\tp2" > $MINIMAP_OUT/populations.txt

	# Run pixy to calculate nucleotide diversity (pi) and d_xy in 10 kb windows
	pixy --vcf $MINIMAP_OUT/${output_prefix}_variants.vcf.gz \
	--populations $MINIMAP_OUT/populations.txt \
	--window_size 10000 \
	--n_cores $SLURM_JOB_CPUS_PER_NODE \
	--output_folder $MINIMAP_OUT \
	--output_prefix pixy_${output_prefix} \
	--stats pi dxy

done

# End of script

