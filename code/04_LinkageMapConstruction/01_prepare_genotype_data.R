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
parents <- subset(pop_metadata, category == "parent", individual, drop = TRUE)
parents
# Names of the F1s
f1s <- subset(pop_metadata, category == "S0", individual, drop = TRUE)
f1s



# Process SNP data ----------------------

# List VCF files
vcf_file_list <- list.files(path = results_dir, pattern = "vcf.gz", full.names = TRUE, recursive = TRUE)

# Iterate
for (vcf_file in vcf_file_list) {

  if (grepl(pattern = "BenLear", x = basename(vcf_file))) {
    reference_genome <-  "BenLear"
  } else if (grepl(pattern = "Oxycoccos", x = basename(vcf_file))) {
    reference_genome <-  "Oxycoccos"
  } else if (grepl(pattern = "Stevens", x = basename(vcf_file))) {
    reference_genome <-  "Stevens"
  }

  # Load the file
  vcf_in <- read.vcfR(file = vcf_file)


  # Name of genotypes
  gt_names <- colnames(vcf_in@gt)
  gt_names <- sub(pattern = "RAPiD-Genomics_F310_", replacement = "", x = gt_names)
  gt_names <- sub(pattern = "_i5-.*", replacement = "", x = gt_names)
  idx <- match(x = gt_names, table = pop_metadata$RG_Sample_Code, nomatch = 0)
  gt_names[-1] <- pop_metadata$individual[idx]

  colnames(vcf_in@gt) <- gt_names


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
  idx <- which(vcf_in1@fix[,"CHROM"] != "Unknown")
  vcf_in1 <- vcf_in1[idx, ]
  vcf_in1

  ##
  ## Sort on chromosome and position
  ##
  vcf_in1 <- vcf_in1[order(vcf_in1@fix[,"CHROM"], as.numeric(vcf_in1@fix[,"POS"])), ]
  unique(vcf_in1@fix[,"CHROM"])


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

  fixed2 <- cbind(fixed1, INFO = info1)

  # Create a newer vcf object
  vcf_in2 <- vcf_in1
  vcf_in2@gt <- gt
  vcf_in2@fix <- fixed2




  ## Remove SNPs that are not on chromosomes -----------------------------------

  max_pos_char <- max(nchar(fixed2[,"POS"]))

  # Fix the chromosome names
  if (reference_genome == "BenLear") {
    fixed2[,"CHROM"] <- sub(pattern = "Vmac_", replacement = "", x = fixed2[,"CHROM"])
  } else {
    fixed2[,"CHROM"] <- paste0("chr", str_pad(string = parse_number(fixed2[,"CHROM"]), width = 2, side = "left", pad = "0"))
  }

  fixed2[,"ID"] <- paste0(fixed2[,"CHROM"], "_", str_pad(string = fixed2[,"POS"], width = max_pos_char, side = "left", pad = "0"))


  idx <- grep(pattern = "^chr[0-9]{2}$", x = fixed2[,"CHROM"])
  fixed2 <- fixed2[idx, ]
  info2 <- info1[idx]

  gt3 <- gt2[idx, ]
  ad3 <- ad1[idx, ]
  dp3 <- dp1[idx, ]

  dim(gt3)



  ## Rename the entries ------------------------------------------------------

  vcf_sample_names <- colnames(gt3)
  new_sample_names <- update_alias(x = vcf_sample_names, alias = as.data.frame(select(pop_metadata, individual, RG_Sample_Code)))

  colnames(gt3) <- new_sample_names
  colnames(ad3) <- new_sample_names
  colnames(dp3) <- new_sample_names




  ## Filter on expected genotype frequencies in the parents and F1 ------------------------------------------------

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
  fixed3 <- fixed2[-idx_rm, ]

  dim(gt4)



  ## Save a new VCF ----------------------------------------------------------

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

  fixed4 <- fixed3[, which(colnames(fixed3) != "INFO")]
  fixed_toprint = cbind(fixed4, INFO = info_toprint)
  filename_out <- sub(pattern = "X", replacement = reference_genome, "data/mxo_variant_cohort_X_alignment_variants_processed.vcf.gz", ignore.case = FALSE)
  write_vcf(filename = filename_out, fixed = fixed_toprint, geno = list(GT = gt4, AD = ad4, DP = DP))

}




# Create r/qtl objects for each family --------------------------------------

for (reference_genome in c("BenLear", "Stevens", "Oxycoccos")) {

  pattern <- sub(pattern = "X", replacement = reference_genome, x = ".*X.*.vcf.gz")

  vcf_filename <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
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

    ####
    #### Impute ####
    ####
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
    # Reorder
    parents_f1_5 <- parents_f1_5[order(row.names(parents_f1_5)), ]
    dim(parents_f1_5)

    # Create a data.frame for import to mappoly
    progeny_gt <- family_gt[row.names(parents_f1_5), setdiff(colnames(family_gt), colnames(parents_f1_5))]
    family_gt1 <- cbind(progeny_gt, parents_f1_5)
    # Remove STEVENS and NJ96-20
    family_gt2 <- family_gt1[, setdiff(colnames(family_gt1), c("STEVENS", "NJ96-20"))]

    # Unique marker names
    mars <- row.names(family_gt2)
    # Ordered GT
    family_gt2 <- family_gt2[mars, ]
    family_fix2 <- vcf_in_family@fix[match(x = mars, table = vcf_in_family@fix[,"ID"]), ]

    vcf_in_family1 <- vcf_in_family
    vcf_in_family1@gt <- cbind(FORMAT = "GT", family_gt2)
    vcf_in_family1@fix <- family_fix2

    # # Convert to a geno DF for mappoly
    # ds <- GT2DS(GT = family_gt2, n.core = 12)
    # geno_df <- data.frame(marker = row.names(ds), parent1 = ds[,f1], parent2 = ds[, f1],
    #                       chrom = family_fix2[,"CHROM"], position = as.numeric(family_fix2[,"POS"]), row.names = NULL)
    # geno_df <- cbind(geno_df, ds[, colnames(progeny_gt)])
    # write_csv(x = geno_df, file = "data/tmp.csv")
    #
    #
    # dat <- mappoly::read_geno_csv(file.in = "data/tmp.csv", ploidy = 2)
    # print(dat, detailed = TRUE)
    #
    # dat1 <- filter_missing(input.data = dat, type = "marker", filter.thres = 0.05, inter = TRUE)
    #
    # dat1 <- filter_missing(input.data = dat1, type = "individual", filter.thres = 0.05, inter = TRUE)
    #
    # pval.bonf <- 0.05/dat1$n.mrk
    # mrks.chi.filt <- filter_segregation(dat1, chisq.pval.thres =  pval.bonf, inter = TRUE)
    # seq.init <- make_seq_mappoly(mrks.chi.filt)
    #
    # all.rf.pairwise <- est_pairwise_rf2(input.seq = seq.init, ncpus = 12)
    #

    onemap_in <- onemap_read_vcfR(vcfR.object = vcf_in_family1, cross = "f2 intercross",
                                  parent1 = "STEVENS_HOM", parent2 = "NJ96-20_HOM", f1 = f1, verbose = TRUE)

    # Output files for r/qtl
    geno_mat <- onemap_in$geno
    geno_mat[geno_mat == 0] <- NA

    # # Get frequencies of B allele
    # geno_mat_b_freq <- colMeans(geno_mat - 1, na.rm = TRUE) / 2
    # plot(geno_mat_b_freq)

    geno_mat[geno_mat == 1] <- "A"
    geno_mat[geno_mat == 2] <- "H"
    geno_mat[geno_mat == 3] <- "B"

    qtl_geno <- as.data.frame(geno_mat)
    qtl_geno <- rbind(as.integer(sub(pattern = "chr", replacement = "", x = onemap_in$CHROM)), onemap_in$POS / 1e6, qtl_geno)
    qtl_geno <- cbind(id = row.names(qtl_geno), qtl_geno)
    qtl_geno$id[1:2] <- ""
    # Save
    write.csv(x = qtl_geno, file = paste0("data/mxo_rqtl_geno_F1-", f1, "_", reference_genome, ".csv"), quote = FALSE, row.names = FALSE)

    qtl_pheno <- data.frame(id = row.names(geno_mat), fake = rnorm(length(row.names(geno_mat))))
    write.csv(x = qtl_pheno, file = paste0("data/mxo_rqtl_pheno_F1-", f1, "_", reference_genome, ".csv"), quote = FALSE, row.names = FALSE)

  }

}







