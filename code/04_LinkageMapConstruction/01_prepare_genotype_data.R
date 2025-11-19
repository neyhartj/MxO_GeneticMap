## MxO Genetic map
##
## Prepare the genotypic data
##
## This code will filter the VCF based on the following parameters:
## 1. Allele depth
## 2. Missingness
##
## The genotypes will be converted for use in OneMap or R/qtl
##

# Load packages
library(tidyverse)
library(vcfR)
library(polyBreedR)
library(onemap)
library(berryBreedR)
library(mappoly)

pop_metadata <- read_csv("data/population_metadata.csv")

cran_dir <- find_cran_dir()


# Minimum allele depth to retain a genotype call
min_AD <- 5
# Maximum missingness per SNP
max_snp_missing <- 0.20


# Names of the parents
parents <- c("STEVENS", "NJ98-20")
# Names of the F1s
f1s <- c("CNJ98-325-33", "CNJ98-309-19")




# Process the Ben Lear reference SNPs -------------------------------------

reference_genome <- "BenLear"

# VCF file to use
vcf_file <- file.path(results_dir, "/variant_calling/variants/gatk/genotype_caller/mxo_variant_cohort_BenLear_alignment_filtered.vcf.gz")

# Load the file
vcf_in <- read.vcfR(file = vcf_file)

##
## Subset genotypes
##
clones_subset <- c("STEVENS", "CNJ98-325-33", "NJ96-20", "CNJ98-309-19",
                   subset(pop_metadata, category == "S1", individual, drop = TRUE))

length(clones_subset)
ind_keep <- intersect(clones_subset, colnames(vcf_in@gt))
length(ind_keep)
setdiff(clones_subset, ind_keep)

vcf_in1 <- vcf_in[,c("FORMAT", ind_keep)]
vcf_in1

##
## Subset chromosomes
##
idx <- which(vcf_in1@fix[,"CHROM"] != "Vmac_Unknown")
vcf_in1 <- vcf_in1[idx, ]
vcf_in1


##
## Subset markers - remove monomorphic
##
gt <- extract.gt(x = vcf_in1, element = "GT")
ds <- GT2DS(GT = gt, n.core = 12)

sdx <- apply(X = ds, MARGIN = 1, FUN = sd, na.rm = TRUE)
idx <- which(sdx > 0)

# Subset
vcf_in1 <- vcf_in1[idx, ]
vcf_in1


## Filter out (i.e. replace with NA) genotype calls with insufficient allele depth

# Extract the AD information
ad <- extract.gt(x = vcf_in1, element = "AD")
# Extract dp information
dp <- extract.gt(x = vcf_in1, element = "DP")
# Extract the GT information; make a copy
gt1 <- gt <- extract.gt(x = vcf_in1, element = "GT")


# Matrix of logicals for hets
is_het <- gt == "0/1" | gt == "1/0"

# Split to ref and alt
ad_ref <- ADsplit(AD = ad, ALT = FALSE)
ad_alt <- ADsplit(AD = ad, ALT = TRUE)

# Apply a minimum depth filter to the AD data strings
# First for homozygous
which_min_ad_hom <- !is_het & (ad_ref >= min_AD | ad_alt >= min_AD)
which_min_ad_het <- is_het & (ad_ref >= min_AD & ad_alt >= min_AD)
which_min_ad <- which_min_ad_hom | which_min_ad_het


# If a sample is NA or does not meet the AD threshold, set as missing (./.)
if (any(is.na(gt1))) {
  gt1[is.na(gt1)] <- "./."
}

if (any(!which_min_ad)) {
  gt1[!which_min_ad] <- "./."
}



## Filter on missingness

# Calculate SNP missingness
snp_missing <- rowMeans(gt1 == "./.")
hist(snp_missing)
# subset
idx <- which(snp_missing <= max_snp_missing)
gt2 <- gt1[idx, , drop = FALSE]
dim(gt2)

ad1 <- ad[idx,]
dp1 <- dp[idx,]

fixed <- getFIX(vcf_in1)
fixed1 <- fixed[idx, ]
info <- getINFO(vcf_in1)
info1 <- info[idx]








# Load the VCF file -------------------------------------------------------

# List the VCF files to use
vcf_file1 <- list.files(path = "data_raw/", pattern = "*.vcf.gz", full.names = TRUE)
vcf_file2 <- file.path(cran_dir, "Genotyping/2023/RAPiD/Data/cran1571_BenLearv1_snps_alias_renamed.vcf.gz")
vcf_files <- c(vcf_file1, vcf_file2)

# Iterate over VCF files
for (file in vcf_files) {
  vcf_in <- read.vcfR(file = file)

  # Determine the reference genome
  if (any(grepl(pattern = "ben", x = vcf_in@meta, ignore.case = T))) {
    reference_genome <- "BenLear"
  } else if (any(grepl(pattern = "stev", x = vcf_in@meta, ignore.case = T))) {
    reference_genome <- "Stevens"
  }

  # Replace sample names with geno ids
  gt_colnames <- colnames(vcf_in@gt)
  idx <- match(x = gt_colnames, table = pop_metadata$RG_Sample_Code, nomatch = 0)
  gt_colnames <- c("FORMAT", pop_metadata$individual[idx])
  vcf_in1 <- vcf_in
  colnames(vcf_in1@gt) <- gt_colnames

# # Subset the VCF for only the relevant clones
# clones_subset <- c("ST", "CNJ98-325-33_rep1", "CNJ98-325-33_rep2", "NJ96-20_rep1", "NJ96-20_rep2", "CNJ98-309-19",
#                    subset(pop_metadata, category == "S1", individual, drop = TRUE))
clones_subset <- c("STEVENS", "CNJ98-325-33", "NJ96-20", "CNJ98-309-19",
                   subset(pop_metadata, category == "S1", individual, drop = TRUE))

ind_keep <- intersect(clones_subset, colnames(vcf_in1@gt))
length(ind_keep)
setdiff(clones_subset, ind_keep)

vcf_in1 <- vcf_in1[,c("FORMAT", ind_keep)]
vcf_in1

# # Check duplicates and retain the one with the lowest missing data
# dup_ind <- grep(pattern = "_rep", x = clones_subset, value = TRUE)
# dup_clones_uniq <- unique(sub(pattern = "_rep.*", replacement = "", x = dup_clones))
# for (dup in dup_clones_uniq) {
#   dup_ind_i <- grep(pattern = dup, x = dup_ind, value = TRUE)
#   dup_gt <- vcf_in1@gt[,dup_ind_i]
#   dup_gt_miss <- colMeans(apply(X = dup_gt, MARGIN = 2, FUN = grepl, pattern = "\\./\\."))
#   dup_ind_i_drop <- dup_ind_i[which.max(dup_gt_miss)]
#   ind_keep <- setdiff(ind_keep, dup_ind_i_drop)
# }
#
# length(ind_keep)
# vcf_in1 <- vcf_in[,c("FORMAT", ind_keep)]
# vcf_in1

# # Rename
# ind_names <- colnames(vcf_in1@gt)[-1]
# ind_names <- sub(pattern = "_rep.*", replacement = "", x = ind_names)
# ind_names <- sub(pattern = "ST", replacement = "STEVENS", x = ind_names, ignore.case = FALSE)
# colnames(vcf_in1@gt)[-1] <- ind_names


# Filter the VCF -----------------------------------------------

# Remove monomorphic markers
gt <- extract.gt(x = vcf_in1, element = "GT")
ds <- GT2DS(GT = gt, n.core = 12)

sdx <- apply(X = ds, MARGIN = 1, FUN = sd, na.rm = TRUE)
idx <- which(sdx > 0)

# Subset
vcf_in1 <- vcf_in1[idx, ]
vcf_in1


## Filter out (i.e. replace with NA) genotype calls with insufficient allele depth

# Extract the AD information
ad <- extract.gt(x = vcf_in1, element = "AD")
# Extract dp information
dp <- extract.gt(x = vcf_in1, element = "DP")
# Extract the GT information; make a copy
gt1 <- gt <- extract.gt(x = vcf_in1, element = "GT")


# Matrix of logicals for hets
is_het <- gt == "0/1" | gt == "1/0"

# Split to ref and alt
ad_ref <- ADsplit(AD = ad, ALT = FALSE)
ad_alt <- ADsplit(AD = ad, ALT = TRUE)

# Apply a minimum depth filter to the AD data strings
# First for homozygous
which_min_ad_hom <- !is_het & (ad_ref >= min_AD | ad_alt >= min_AD)
which_min_ad_het <- is_het & (ad_ref >= min_AD & ad_alt >= min_AD)
which_min_ad <- which_min_ad_hom | which_min_ad_het


# If a sample is NA or does not meet the AD threshold, set as missing (./.)
if (any(is.na(gt1))) {
  gt1[is.na(gt1)] <- "./."
}

if (any(!which_min_ad)) {
  gt1[!which_min_ad] <- "./."
}



## Filter on missingness

# Calculate SNP missingness
snp_missing <- rowMeans(gt1 == "./.")
hist(snp_missing)
# subset
idx <- which(snp_missing <= max_snp_missing)
gt2 <- gt1[idx, , drop = FALSE]
dim(gt2)

ad1 <- ad[idx,]
dp1 <- dp[idx,]

fixed <- getFIX(vcf_in1)
fixed1 <- fixed[idx, ]
info <- getINFO(vcf_in1)
info1 <- info[idx]

# Remove SNPs that are not on chromosomes -----------------------------------

max_pos_char <- max(nchar(fixed1[,"POS"]))

# Fix the chromosome names
if (reference_genome == "BenLear") {
  fixed1[,"CHROM"] <- sub(pattern = "Vmac_", replacement = "", x = fixed1[,"CHROM"])
} else {
  fixed1[,"CHROM"] <- paste0("chr", str_pad(string = parse_number(fixed1[,"CHROM"]), width = 2, side = "left", pad = "0"))
}

fixed1[,"ID"] <- paste0(fixed1[,"CHROM"], "_", str_pad(string = fixed1[,"POS"], width = max_pos_char, side = "left", pad = "0"))


idx <- grep(pattern = "^chr[0-9]{2}$", x = fixed1[,"CHROM"])
fixed2 <- fixed1[idx, ]
info2 <- info1[idx]

gt3 <- gt2[idx, ]
ad3 <- ad1[idx, ]
dp3 <- dp1[idx, ]

dim(gt3)

# Rename the entries ------------------------------------------------------

vcf_sample_names <- colnames(gt3)
new_sample_names <- update_alias(x = vcf_sample_names, alias = as.data.frame(select(pop_metadata, individual, RG_Sample_Code)))

colnames(gt3) <- new_sample_names
colnames(ad3) <- new_sample_names
colnames(dp3) <- new_sample_names




# Filter on expected genotype frequencies in the parents and F1 ------------------------------------------------

# First remove SNPs where STEVENS is homozygous alternate
#
# DO NOT DO THIS IF THE REFERENCE IS BEN LEAR OR OXYCOCCOS

if (reference_genome == "BenLear") {
  stevens_hom_alt <- integer(0)
} else {
  stevens_hom_alt <- which(gt3[,"STEVENS"] == "1/1")
}



# Next remove SNPs where the F1s are not hets
f1_not_het <- gt3[, f1s] != "0/1" & gt3[, f1s] != "0|1"
f1_not_het <- which(apply(X = f1_not_het, MARGIN = 1, FUN = all))

idx_rm <- union(stevens_hom_alt, f1_not_het)
length(idx_rm)

# Filter
gt4 <- gt3[-idx_rm, ]
ad4 <- ad3[-idx_rm, ]
dp4 <- dp3[-idx_rm, ]
fixed4 <- fixed2[-idx_rm, ]

# # Subset only known chromosomes
# idx <- grep(pattern = "^chr[0-9]{2}$", x = fixed3[,"CHROM"])
# gt4 <- gt4[idx, ]
# ad4 <- ad4[idx, ]
# dp4 <- dp4[idx, ]
# fixed4 <- fixed3[idx, ]


dim(gt4)

# Save a new VCF ----------------------------------------------------------

# Recalculate info
DP <- apply(ad4, 2, function(x) {
  sapply(strsplit(x, split = ",", fixed = T), function(u) {
    sum(as.integer(u))
  })
})
NS <- rowSums(DP > 0, na.rm = TRUE)
DP.AVG <- rowMeans(DP, na.rm = TRUE)
alt <- apply(ad4, 1, function(x) { sum(as.integer(sub(pattern = "([0-9]{1,})(,)([0-9]{1,})", replacement = "\\3", x = x)), na.rm = TRUE) })
AF <- round(alt/rowSums(DP, na.rm = TRUE), 3)
AF[is.na(AF)] <- "."
info_toprint <- polyBreedR:::make_info(cbind(NS = NS, DP.AVG = round(DP.AVG, 1), AF = AF))

fixed_toprint = cbind(fixed4, INFO = info_toprint)
filename_out <- "data/mxo_benlear-ref_variants_processed.vcf.gz"
write_vcf(filename = filename_out, fixed = fixed_toprint, geno = list(GT = gt4, AD = ad4, DP = DP))




# Create vcf objects for each family --------------------------------------

vcf_filename <- filename_out
# vcf_in2 <- read.vcfR(file = vcf_filename)
vcf_in2 <- read.vcfR(file = vcf_filename)
vcf_in2

# split pop metadata by the F1 individual
pop_by_family <- pop_metadata %>%
  split(.$S0_name) %>%
  purrr::map(~select(.x, individual, parent1, parent2, S0_name) %>% unlist() %>% unique() %>% intersect(., colnames(vcf_in2@gt)))

# For each family, create a onemap object; then convert it for use in R/qtl
for (family in pop_by_family) {
  vcf_in_family <- vcf_in2[, c("FORMAT", family)]
  f1 <- intersect(family, f1s)

  family_gt <- extract.gt(x = vcf_in_family, element = "GT")
  # Convert phased to unphased gt
  family_gt <- sub(pattern = "\\|", replacement = "/", x = family_gt)
  # Extract the GT data for the parents and the F1
  parents_f1 <- family_gt[,c("STEVENS", "NJ96-20", f1)]

  # Impute
  parents_f1_1 <- parents_f1
  # Any F1 genotype is a het if Stevens is 0/0 and Oxy is 1/1 or vice versa
  stevens_00_oxy_11 <- parents_f1_1[,"STEVENS"] == "0/0" & parents_f1_1[,"NJ96-20"] == "1/1"
  sum(stevens_00_oxy_11, na.rm = TRUE)
  parents_f1_1[stevens_00_oxy_11, f1] <- "0/1"

  # Use homozygous alternate for stevens if the reference genome is not stevens
  if (reference_genome != "Stevens") {
    stevens_11_oxy_00 <- parents_f1_1[,"STEVENS"] == "1/1" & parents_f1_1[,"NJ96-20"] == "0/0"
    sum(stevens_11_oxy_00, na.rm = TRUE)
    parents_f1_1[stevens_11_oxy_00, f1] <- "0/1"
    dim(parents_f1_1)
  } else if (reference_genome == "Stevens") {
    # If Stevens, set all Stevens genotypes that are not 0/0 to missing
    stevens_01 <- parents_f1_1[,"STEVENS"] == "0/1"
    parents_f1_1[stevens_01, "STEVENS"] <- "./."
  }

  # Remove missing and non-hets in the f1
  parents_f1_2 <- parents_f1_1[!is.na(parents_f1_1[, f1]), ]
  parents_f1_3 <- parents_f1_2[parents_f1_2[, f1] == "0/1", ]
  dim(parents_f1_3)

  # Attempt to reconstruct the parent haplotypes using the F1
  # i.e. switch heterozygous genotypes in the parent to homozygous depending on the genotypes
  # of the F1
  parents_f1_4 <- cbind(parents_f1_3, STEVENS_HOM = as.character(NA), "NJ96-20_HOM" = as.character(NA))
  parents_f1_4[parents_f1_4[,"STEVENS"] == "0/0", "STEVENS_HOM"] <- "0/0"
  parents_f1_4[parents_f1_4[,"NJ96-20"] == "1/1", "NJ96-20_HOM"] <- "1/1"
  parents_f1_4[parents_f1_4[,"STEVENS"] == "1/1", "STEVENS_HOM"] <- "1/1"
  parents_f1_4[parents_f1_4[,"NJ96-20"] == "0/0", "NJ96-20_HOM"] <- "0/0"

  parents_f1_4[parents_f1_4[,"STEVENS"] == "0/1" & parents_f1_4[,"NJ96-20"] == "1/1", "STEVENS_HOM"] <- "0/0"
  parents_f1_4[parents_f1_4[,"STEVENS"] == "0/1" & parents_f1_4[,"NJ96-20"] == "0/0", "STEVENS_HOM"] <- "1/1"

  parents_f1_4[parents_f1_4[,"STEVENS"] == "1/1" & parents_f1_4[,"NJ96-20"] == "0/1", "NJ96-20_HOM"] <- "0/0"
  parents_f1_4[parents_f1_4[,"STEVENS"] == "0/0" & parents_f1_4[,"NJ96-20"] == "0/1", "NJ96-20_HOM"] <- "1/1"

  # Remove any markers that are not STEVENS == 0/0 and OXY == 1/1
  idx <- union(which(parents_f1_4[,"STEVENS_HOM"] == "0/0" & parents_f1_4[,"NJ96-20_HOM"] == "1/1"),
               which(parents_f1_4[,"STEVENS_HOM"] == "1/1" & parents_f1_4[,"NJ96-20_HOM"] == "0/0"))
  parents_f1_5 <- parents_f1_4[idx, ]
  dim(parents_f1_5)

  # Replace these parent GTs in the vcf
  family_gt1 <- cbind(family_gt[row.names(parents_f1_5), setdiff(colnames(family_gt), colnames(parents_f1_5))], parents_f1_5)
  # Remove STEVENS and NJ96-20
  family_gt2 <- family_gt1[, setdiff(colnames(family_gt1), c("STEVENS", "NJ96-20"))]
  vcf_in_family1 <- vcf_in_family
  vcf_in_family1@gt <- cbind(FORMAT = "GT", family_gt2)
  vcf_in_family1@fix <- vcf_in_family1@fix[vcf_in_family1@fix[, "ID"] %in% row.names(family_gt2), ]

  onemap_in <- onemap_read_vcfR(vcfR.object = vcf_in_family1, cross = "f2 intercross",
                                parent1 = "STEVENS_HOM", parent2 = "NJ96-20_HOM", f1 = f1, verbose = TRUE)

  # Output files for r/qtl
  geno_mat <- onemap_in$geno
  geno_mat[geno_mat == 0] <- NA
  geno_mat[geno_mat == 1] <- "A"
  geno_mat[geno_mat == 2] <- "H"
  geno_mat[geno_mat == 3] <- "B"

  qtl_geno <- as.data.frame(geno_mat)
  qtl_geno <- rbind(as.integer(sub(pattern = "chr", replacement = "", x = onemap_in$CHROM)), onemap_in$POS / 1e6, qtl_geno)
  qtl_geno <- cbind(id = row.names(qtl_geno), qtl_geno)
  qtl_geno$id[1:2] <- ""
  # Save
  write.csv(x = qtl_geno, file = paste0("data/mxo_rqtl_geno_F1-", f1, ".csv"), quote = FALSE, row.names = FALSE)

  qtl_pheno <- data.frame(id = row.names(geno_mat), fake = rnorm(length(row.names(geno_mat))))
  write.csv(x = qtl_pheno, file = paste0("data/mxo_rqtl_pheno_F1-", f1, ".csv"), quote = FALSE, row.names = FALSE)

}







