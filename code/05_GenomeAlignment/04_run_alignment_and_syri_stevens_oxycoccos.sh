#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

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


# Change the working directory
cd /project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap/analysis/GenomeComparison

# Load modules
module load miniconda
module load minimap2

# Load the environment
source activate wgs_env


# Set variables
# Stevens reference
stevens=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_macrocarpon_Stevens_v1.fasta
# Oxycoccos assembly
oxyassembly=/project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap/analysis/GenomeAssembly/oxycoccos_ont_flye_allmaps_assembly_chr.fasta
oxyassembly=/project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap/analysis/GenomeAssembly/Vaccinium_oxycoccos_Kawash2022_assembly_db_allmaps_scaffolding/Vaccinium_oxycoccos_Kawash2022_assembly_db_allmaps_scaffolding_assembly_chr.fasta

OUTDIR=Stevens_vs_voxy_kawash2022_rescaffolded


NTHREADS=$SLURM_JOB_CPUS_PER_NODE



#####################################################################################################################


# Change the headers in the stevens reference and remove any non-chromosome contigs
newstevens=Vaccinium_macrocarpon_Stevens_v1_renamed.fasta
awk '{
if ($0 ~ /^>chr[0-9]+$/) {
match($0, /^>chr([0-9]+)$/, m)
if (m[1] + 0 < 10) {
print ">chr0" m[1]
} else {
print $0
}
} else {
print $0
}
}' $stevens > $newstevens

newstevens2=Vaccinium_macrocarpon_Stevens_v1_renamed_chr.fasta
awk '/^>/ {keep = ($0 ~ /^>chr/)} keep' $newstevens > $newstevens2

# Remove the intermediate file
rm $newstevens


##

# Align the oxycoccos genome to the stevens using stevens as the reference
nucmer --maxmatch -c 500 -b 500 -l 100 -p $OUTDIR $newstevens2 $oxyassembly

# Filter
delta-filter -m -i 90 -l 100 $OUTDIR.delta > ${OUTDIR}_i90_l100.delta

# Convert to tab of coordinates
show-coords -THrd ${OUTDIR}_i90_l100.delta > ${OUTDIR}_i90_l100.coords

# Run SYRI
prefix=${OUTDIR}_i90_l100_syri_

syri -r $newstevens2 -q $oxyassembly \
	-c ${OUTDIR}_i90_l100.coords -d ${OUTDIR}_i90_l100.delta \
	--lf ${prefix}.log --log DEBUG \
	-k \
	--nc 12 \
	--prefix $prefix

#awk '$11=="INV" {print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}' ${prefix}.syri.out > ${prefix}.inversions.tsv
#awk '$11=="TRANS" || $11=="INVTR" {print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}' ${prefix}.syri.out > ${prefix}.translocations.tsv

# Visualize
echo -e "#file\tname\ttags\n${newstevens2}\tstevens\n${oxyassembly}\toxycoccos" > ${OUTDIR}_genomes.txt
plotsr --sr ${prefix}_syri.out --genomes ${OUTDIR}_genomes.txt -o ${prefix}_plot.png



