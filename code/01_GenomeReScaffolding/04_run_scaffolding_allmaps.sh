#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # X processor core(s) per node X 2 threads per core
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --job-name="scaffold genome with allmaps"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change the working directory
PROJ=/project/gifvl_vaccinium/cranberryGenotyping/MxO_GeneticMap
SUBDIR=$PROJ/analysis/GenomeAssembly

cd $SUBDIR

# Load modules
module load miniconda
module load blast+
module load r

NTHREADS=$SLURM_JOB_CPUS_PER_NODE

# Path to the draft assembly
# INPUTASSEMBLY=oxycoccos_ONT_flye_out/assembly.fasta
INPUTASSEMBLY=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_oxycoccos_physical_v1.fasta

# Path to the marker contextual sequences
MARKERFASTA=/project/gifvl_vaccinium/cranberryGenotyping/CranberryConsensusLinkageMap/output/blast_output/unique_scaffold_map_positions_flanking_sequence.fasta

ASSEMBLYNAME=Vaccinium_oxycoccos_Kawash2022_assembly_db

##############################################################################################

# ASSEMBLYDIR=oxycoccos_ONT_flye_out
ASSEMBLYDIR=$SUBDIR

OUTDIR=${ASSEMBLYNAME}_allmaps_scaffolding

### BLAST marker sequences to the draft assembly

# Make a blast DB
makeblastdb -in $INPUTASSEMBLY -dbtype nucl -out $ASSEMBLYDIR/$ASSEMBLYNAME

# Run blast
blastn -query $MARKERFASTA -db $ASSEMBLYDIR/$ASSEMBLYNAME -outfmt "6 qseqid sseqid sstart send pident length evalue bitscore" -evalue 1e-10 -perc_identity 95 -num_threads $NTHREADS > $ASSEMBLYDIR/${ASSEMBLYNAME}_blast.tsv

# Subset the blast output with an R script and create the input file for allmaps
Rscript 03a_merge_map_with_blast.R $ASSEMBLYDIR/${ASSEMBLYNAME}_blast.tsv composite_linkage_map.csv ${OUTDIR}/allmaps_input.csv


# Load the environment
source activate wgs_env

# Convert the CSV files to a bed file
python -m jcvi.assembly.allmaps merge ${OUTDIR}/allmaps_input.csv -o ${OUTDIR}/allmaps_input.bed

# Run scaffolding using allmaps
python -m jcvi.assembly.allmaps path ${OUTDIR}/allmaps_input.bed $INPUTASSEMBLY

# Rename the new assembly
mv ${OUTDIR}/allmaps_input.chr.fasta ${OUTDIR}/${OUTDIR}_assembly_chr.fasta

# Run stats on the new assembly
stats.sh -Xmx256g in=${OUTDIR}/${OUTDIR}_assembly_chr.fasta out=${OUTDIR}/${OUTDIR}_assembly_chr_fasta_stats.out
