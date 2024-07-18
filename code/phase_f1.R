## MxO Genetic Map
##
## Phase the F1 individuals
## Construct a genetic map
##
##

library(tidyverse)
library(vcfR)
library(mappoly)
library(polyBreedR)
library(onemap)


# Load population data ----------------------------------------------------

pop_data <- read_csv(file = file.path(data_dir, "mxo_rapid15k_sample_information.csv"))

# Subset the two parents and the F1
indiv_names <- c("CNJ98-325-33", "NJ96-20", "STEVENS")
pop_data1 <- subset(pop_data, geno_id %in% indiv_names)


# Load the vcf data -------------------------------------------------------

# steven alignment
trio_stevens_ref_vcf <- read.vcfR(file = "data_raw/filtered_stevens_output.recode.vcf")
# oxy alignment
trio_oxy_ref_vcf <- read.vcfR(file = "data_raw/filtered_oxy_output.recode.vcf")

# Read in the entire VCF
f2_vcf_stevens <- read.vcfR(file = "data_raw/mxo_stevens-ref_variants_filtered.vcf.gz")
# Rename the VCF
colnames(f2_vcf_stevens@gt) <- c("FORMAT", pop_data$geno_id[match(x = colnames(f2_vcf_stevens@gt), table = pop_data$RG_Sample_Code, nomatch = 0)])
# Remove CNJ03-79 invididuals
idx_remove <- grep(pattern = "CNJ03-79", x = colnames(f2_vcf_stevens@gt))
f2_vcf_stevens <- f2_vcf_stevens[,-idx_remove]

f2_vcf_oxy <- read.vcfR(file = "data_raw/mxo_oxy-ref_variants_filtered.vcf.gz")
# Rename the VCF
colnames(f2_vcf_oxy@gt) <- c("FORMAT", pop_data$geno_id[match(x = colnames(f2_vcf_oxy@gt), table = pop_data$RG_Sample_Code, nomatch = 0)])
# Remove CNJ03-79 invididuals
idx_remove <- grep(pattern = "CNJ03-79", x = colnames(f2_vcf_oxy@gt))
f2_vcf_oxy <- f2_vcf_oxy[,-idx_remove]



# Phase the F1 using unambiguous markers ----------------------------------

trio_stevens_ref_geno <- extract.gt(x = f2_vcf_stevens[,c("FORMAT", pop_data1$geno_id)], element = "GT")

# First find heterozygous markers in the F1
# Set all others to NA
which_het <- which(trio_stevens_ref_geno[,"CNJ98-325-33"] == "0/1")
trio_stevens_ref_geno[-which_het,"CNJ98-325-33"] <- NA




f1_geno_mat <- trio_stevens_ref_geno[,"CNJ98-325-33", drop = FALSE]
f1_geno_mat[f1_geno_mat[,] == "0/0", ] <- NA




# Unambiguous markers are the following:
# 1. Contrasting homozygotes in the parents
# 2. Hets in the F1

# Find contrasting homozygotes in the parents
parents_stevens_ref_geno <- trio_stevens_ref_geno[,c("NJ96-20", "STEVENS")]

which_contrasting_homos <- parents_stevens_ref_geno[,1] == "0/0" & parents_stevens_ref_geno[,2] == "1/1"
# Which markers have OXY with 1/1 (homozygous alternate)? These would imply that stevens would be homo ref
#






# Create an empty matrix to store haplotype information
# 1 = stevens
# 2 = oxy
# f1_hap_matrix <- matrix(as.numeric(NA)




# Use onemap to construct the genetic map ---------------------------------




vcf_f2 <- onemap_read_vcfR(vcfR.object = f2_vcf_stevens, cross = "f2 intercross",
                           parent1 = indiv_names[2], parent2 = indiv_names[3], f1 = indiv_names[1],
                           only_biallelic = TRUE, verbose = TRUE)




# Use the mappoly vignette to construct a genetic map ---------------------

# Convert the VCF to a dosage file
vcf2csv(vcf.file = )


dat <- read_geno




















