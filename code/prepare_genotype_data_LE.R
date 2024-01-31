## MxO S1 Mapping
##
## Genotyping
##
## Code to filter the VCF file
##

## vcftools --vcf Final_DP10_Corrected_VAC_155005_SNPs.vcf --remove-indels --min-alleles 2 \
## --max-alleles 2 --keep mxo_vcf_sample_list.txt --minQ 40 --recode --recode-INFO-all \
## --out Final_DP10_Corrected_VAC_155005_SNPs_biallelic_mxos1_qual



# Load packages
library(tidyverse)
library(vcfR)
library(snps)
library(statgenIBD)
library(onemap)

# Load population information
pop_metadata <- read_csv("Data/population_metadata.csv")


# Load the VCF file -------------------------------------------------------

vcf_in <- read.vcfR(file = "C:/Users/Lindsay.Erndwein/OneDrive - USDA/Documents/R/hybrid_oxy-ref_rapid_merge.vcf.tar.gz")

# Get snp info
snp_info <- vcf_in@fix %>%
  as.data.frame() %>%
  mutate(chrom_num = parse_number(CHROM),
         chrom_num = str_pad(chrom_num, 2, "left", "0"),
         marker = paste0("S", chrom_num, "_", POS)) %>%
  select(marker, chrom = CHROM, pos = POS, ref = REF, alt = ALT) %>%
  mutate_if(is.character, parse_guess) %>%
  mutate(chrom = factor(chrom, levels = paste0("chr", 1:12)))



# Get genotype calls
gt <- vcf_in@gt

# Function to convert to numeric
convert_gt <- function(x) {
  x1 <- sapply(strsplit(x, ":"), "[[", 1)
  x2 <- x1
  x2[x2 == "./."] <- NA
  x2[x2 == "0/0"] <- 2
  x2[x2 == "0/1"] <- 1
  x2[x2 == "1/1"] <- 0
  as.numeric(x2)
}

gt2 <- apply(X = gt[,-1, drop = FALSE], MARGIN = 2, FUN = convert_gt)
rownames(gt2) <- snp_info$marker
colnames(gt2) <- pop_metadata$individual[match(x = colnames(gt2), table = pop_metadata$RG_Sample_Code)]

geno_mat <- t(gt2)

# Save
filename <- file.path("Results/mxos1_15k_recoded_genotypes.RData")
save("geno_mat", "snp_info", file = filename)

if (file.exists(filename)) load(filename)



# Missingness
missingness <- calc_missing(x = geno_mat - 1)
snp_miss <- colMeans(is.na(geno_mat))
geno_miss <- rowMeans(is.na(geno_mat))

# Filter on missingness
geno_mat1 <- geno_mat[geno_miss < 0.20, snp_miss < 0.20]

# Filter out completely monomorphic SNPs
maf <- calc_maf(x = geno_mat1 - 1, check.matrix = TRUE)
geno_mat2 <- geno_mat1[,maf > 0, drop = FALSE]


## Filter out completely correlated markers
# First split by chromosome
snp_info1 <- snp_info %>%
  subset(marker %in% colnames(geno_mat2))

snps_per_chrom <- snp_info1 %>%
  split(.$chrom)

# Maximum SNP distance to separate into haplogroups
snp_haplo_max_dist <- 300
# Maximum SNP similarity to remove snps
max_snp_sim <- 0.95

# Group SNPs by haplotypes
snp_haplos_per_chrom <- snps_per_chrom %>%
  purrr::map_df(~{
    # SNP distance
    x_pos <- setNames(.x$pos, .x$marker)
    x_clust <- hclust(dist(x_pos))
    x_group <- cutree(tree = x_clust, h = snp_haplo_max_dist)
    # Assign group to snps
    mutate(.x, group = paste0(unique(.x$chrom), "_", x_group))
  })

# Filter SNPs within haplotypes for LD and maf
pruned_snps_per_haplo <- snp_haplos_per_chrom %>%
  distinct() %>%
  split(.$group) %>%
  purrr::map_df(~{
    if (nrow(.x) == 1) {
      .x
    } else {
      # Calculate pairwise similarity between snps
      x_mat <- geno_mat2[,.x$marker]
      x_pairs <- combn(x = .x$marker, m = 2)
      x_sim <- apply(X = x_pairs, MARGIN = 2, FUN = function(snps) mean(x_mat[,snps[1]] == x_mat[,snps[2]], na.rm = T))
      x_sim_df <- data.frame(marker1 = x_pairs[1,], marker2 = x_pairs[2,], similarity = x_sim) %>%
        spread(marker2, similarity) %>%
        column_to_rownames("marker1") %>%
        as.matrix()
      x_sim_df1 <- rbind(cbind(NA, x_sim_df), NA)
      dimnames(x_sim_df1) <- replicate(2, c(row.names(x_sim_df), tail(colnames(x_sim_df1), 1)), simplify = FALSE)
      diag(x_sim_df1) <- 1
      x_sim_df2 <- as.matrix(Matrix::forceSymmetric(x_sim_df1))

      # Prune
      x_pruned <- prune_LD(x = x_mat - 1, sim.mat = x_sim_df2, r2.max = max_snp_sim)
      # Return the snp info for pruned snps
      subset(.x, marker %in% colnames(x_pruned))

    }

  })


# Edit snp_info
snp_info2 <- pruned_snps_per_haplo %>%
  arrange(chrom, pos)

# Subset geno_mat
geno_mat3 <- geno_mat2[, snp_info2$marker, drop = FALSE]

# MAF
snp_maf <- calc_maf(x = geno_mat3 - 1)
# Filter on MAF
geno_mat4 <- geno_mat3[, snp_maf >= (15 / nrow(geno_mat3)), drop = FALSE]

# Subset these markers from the vcf and save a new vcf
vcf_out <- vcf_in
vcf_out@fix[,"ID"] <- paste0("S", str_pad(parse_number(vcf_out@fix[,"CHROM"]), 2, "left", "0"), "_", vcf_out@fix[,"POS"])
# Index of markers to keep
marker_idx_keep <- which(vcf_out@fix[,"ID"] %in% colnames(geno_mat4))
vcf_out@fix <- vcf_out@fix[marker_idx_keep, ]

gt1 <- gt[marker_idx_keep, c("FORMAT", subset(pop_metadata, individual %in% c(offspring_names, row.names(parents2)), RG_Sample_Code, drop = TRUE))]
colnames(gt1) <- c("FORMAT", pop_metadata$individual[match(colnames(gt1)[-1], table = pop_metadata$RG_Sample_Code)])

vcf_out@gt <- gt1

# Save the vcf
write.vcf(x = vcf_out, file = "Data/mxo_offspring_parnet_geno.vcf")


# Assign alleles to parents -----------------------------------------------

# Subset the parents
parents <- geno_mat4[c("STEVENS", "NJ96-20"), , drop = FALSE]
apply(X = parents, MARGIN = 1, FUN = table)

# Remove NA
parents1 <- parents[, colMeans(is.na(parents)) == 0, drop = FALSE]
# Subset markers where one parent is 2 and the other is 0
which_seg_mar <- (parents1[1,] == 2 & parents1[2,] == 0) | (parents1[1,] == 0 & parents1[2,] == 2)

parents2 <- parents1[, which_seg_mar]


## Assign A B or H to markers in the offspring
# A = stevens
# B = NJ96-20
# H = het

offspring_names <- intersect(row.names(geno_mat4), subset(pop_metadata, category == "S1" & S0_name == "CNJ98-325-33", individual, drop = TRUE))
offspring1 <- geno_mat4[offspring_names, colnames(parents1)]
offspring2 <- offspring1[, colnames(parents2)]

# Iterate over snps in the offspring
offspring_parent_alleles <- offspring2
for (j in seq_len(ncol(offspring_parent_alleles))) {
  snp_j <- as.character(offspring2[,j])
  # 0 or 2 for A allele?
  a_allele <- parents2["STEVENS",j]
  b_allele <- setdiff(c(0, 2), a_allele)
  snp_j[snp_j == a_allele] <- "A"
  snp_j[snp_j == "1"] <- "H"
  snp_j[snp_j == b_allele] <- "B"
  offspring_parent_alleles[,j] <- snp_j
}

snp_info4 <- snp_info2 %>%
  subset(marker %in% colnames(offspring_parent_alleles)) %>%
  arrange(chrom, pos)

offspring_parent_alleles <- offspring_parent_alleles[,snp_info4$marker]


# Save this data in CSV format for qtl
data_to_save <- rbind(
  c("pheno", colnames(offspring_parent_alleles)),
  c("", as.numeric(snp_info4$chrom)),
  c("", as.numeric(snp_info4$pos)),
  cbind(seq_len(nrow(offspring_parent_alleles)), offspring_parent_alleles)
)

write_csv(x = as.data.frame(data_to_save), file = "Data/mxo_offspring_geno_qtl.csv", col_names = FALSE)


# Load genotypes in QTL and estimate recombination freqs ------------------

# Read data to qtl
mxo_cross <- read.cross(format = "csv", file = "Data/mxo_offspring_geno_qtl.csv")

# Calculate segregation pattern
mxo_cross_seg <- geno.table(cross = mxo_cross) %>%
  rownames_to_column("marker") %>%
  left_join(., snp_info4)




# A function to measure the recombination fraction between adjacent markers per chromosom
est_rf_adj <- function(cross, chr) {
  # Get marker names per chromosome
  marker_names <- markernames(cross = cross, chr = chr)
  # Empty matrix
  marker_rf <- matrix(NA, nrow = length(marker_names), ncol = length(marker_names),
                      dimnames = list(marker_names, marker_names))

  # Iterate over i to length(marker_names)-1
  for (i in seq_len(length(marker_names) - 1)) {
    ii <- i + 1
    markers_pull <- marker_names[c(i, ii)]
    markers_pull_rf <- est.rf(cross = pull.markers(cross = cross, markers = markers_pull))$rf
    marker_rf[markers_pull[1], markers_pull[2]] <- markers_pull_rf[2,1]
  }

  # Convert to map distance
  marker_d <- qtl:::imf.h(marker_rf)


}






sim_cross <- sim.cross(map = sim.map(len = rep(100, 12), n.mar = 50))
sim_cross_rf <- est.rf(cross = sim_cross)
# Calculate haldane mapping distance
# sim_cross_rf_map <-

# Calculate rf between first two markers
sim_cross_rf_ij <- est.rf(pull.markers(cross = sim_cross, markers = markernames(sim_cross)[1:2]))
sim_cross_d_ij <- qtl::imf.h(sim_cross_rf_ij$rf[2,1])





# Read in using onemap ----------------------------------------------------


onemap_data_in <- onemap_read_vcfR()













snp_af1 <- snp_af[colnames(parents2)]

inner_join(snp_info, data.frame(marker = names(snp_af1), ref_freq = snp_af1)) %>%
  ggplot(aes(x = pos, y = ref_freq)) +
  geom_point() +
  facet_wrap(~ chrom)



# K matrix
K <- rrBLUP::A.mat(X = geno_mat2 - 1, min.MAF = 0, max.missing = 1, impute.method = "mean", return.imputed = TRUE)

# PCA
K_pca <- prcomp(x = K$A)
K_pca_df <- broom::tidy(K_pca) %>%
  mutate(PC = paste0("PC", PC))

K_pca_df %>%
  subset(PC %in% c("PC1", "PC2")) %>%
  spread(PC, value) %>%
  mutate(parent = ifelse(str_detect(row, "CNJ"), NA, row)) %>%
  ggplot(aes(x = PC1, y = PC2, color = parent)) +
  geom_point()

