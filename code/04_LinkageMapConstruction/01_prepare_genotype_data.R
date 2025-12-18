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

# get the F1 / S0 names
s0_names <- subset(pop_metadata, category == "S0", individual, drop = TRUE)


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
  vcf_in1 <- vcf_in
  vcf_in1

  # Name of genotypes
  gt_names <- colnames(vcf_in1@gt)
  gt_names <- sub(pattern = "RAPiD-Genomics_F310_", replacement = "", x = gt_names)
  gt_names <- sub(pattern = "HC_RAPiD-Genomics_F368_", replacement = "", x = gt_names)
  gt_names <- sub(pattern = "_i5-.*", replacement = "", x = gt_names)
  idx <- match(x = gt_names, table = pop_metadata$RG_Sample_Code, nomatch = 0)
  gt_names[-1] <- pop_metadata$individual[idx]

  colnames(vcf_in1@gt) <- gt_names


  for (s0 in s0_names) {

    # Identify the clones to subset for this family
    clones_subset <- subset(pop_metadata, S0_name == s0)
    parent1 <- unique(clones_subset$parent1)
    parent2 <- unique(clones_subset$parent2)
    f2s <- clones_subset$individual
    clones_subset <- c(parent1, parent2, s0, f2s)

    ind_keep <- intersect(clones_subset, colnames(vcf_in1@gt))

    # Print
    cat("\nParents kept in the vcf file:", intersect(ind_keep, c(parent1, parent2)))
    cat("F1 kept in the vcf file:", intersect(ind_keep, s0))
    cat("Number of F2s kept in the vcf file:", length(intersect(ind_keep, f2s)), "/", length(f2s))

    vcf_in2 <- vcf_in1[,c("FORMAT", ind_keep)]
    vcf_in2


    ##
    ## Subset chromosomes
    ##
    idx <- grep(pattern = "Unknown", x = vcf_in2@fix[,"CHROM"], invert = TRUE)
    vcf_in2 <- vcf_in2[idx, ]
    vcf_in2

    ##
    ## Sort on chromosome and position
    ##
    vcf_in2 <- vcf_in2[order(vcf_in2@fix[,"CHROM"], as.numeric(vcf_in2@fix[,"POS"])), ]
    vcf_in2

    ##
    ## Subset markers - remove monomorphic
    ##
    gt <- extract.gt(x = vcf_in2, element = "GT")
    ds <- GT2DS(GT = gt, n.core = 12)

    sdx <- apply(X = ds, MARGIN = 1, FUN = sd, na.rm = TRUE)
    idx <- which(sdx > 0)

    # Subset
    vcf_in2 <- vcf_in2[idx, ]
    vcf_in2


    ## Filter out (i.e. replace with NA) genotype calls with insufficient allele depth

    # Extract the AD information
    ad <- extract.gt(x = vcf_in2, element = "AD")
    # Extract dp information
    dp <- extract.gt(x = vcf_in2, element = "DP")
    # Extract the GT information; make a copy
    gt1 <- gt <- extract.gt(x = vcf_in2, element = "GT")


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

    fixed <- getFIX(vcf_in2)
    fixed1 <- fixed[idx, ]
    info <- getINFO(vcf_in2)
    info1 <- info[idx]

    fixed2 <- cbind(fixed1, INFO = info1)

    # Create a newer vcf object
    vcf_in3 <- vcf_in2
    vcf_in3@gt <- gt
    vcf_in3@fix <- fixed2
    vcf_in3



    ## Remove SNPs that are not on chromosomes

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



    ## Rename the entries

    vcf_sample_names <- colnames(gt3)
    new_sample_names <- update_alias(x = vcf_sample_names, alias = as.data.frame(select(pop_metadata, individual, RG_Sample_Code)))

    colnames(gt3) <- new_sample_names
    colnames(ad3) <- new_sample_names
    colnames(dp3) <- new_sample_names




    # Next remove SNPs where the F1s are not hets
    f1_het <- gt3[, s0] == "0/1" | gt3[, s0] == "0|1"
    idx <- which(f1_het)
    length(idx)

    # Filter
    gt4 <- gt3[idx, ]
    ad4 <- ad3[idx, ]
    dp4 <- dp3[idx, ]
    fixed3 <- fixed2[idx, ]

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

    filename_out <- paste0("data/mxo_variants_processed_", reference_genome, "_family_", s0, ".vcf.gz")
    write_vcf(filename = filename_out, fixed = fixed_toprint, geno = list(GT = gt4, AD = ad4, DP = DP))

  }

}




# Create r/qtl objects for each family --------------------------------------

# List the VCF files
vcf_file_list <- list.files(path = data_dir, pattern = "mxo_variants_processed", full.names = TRUE)

# Iterate over these files
for (vcf_file in vcf_file_list) {

  # Which reference genome
  if (grepl(pattern = "BenLear", x = basename(vcf_file))) {
    reference_genome <- "BenLear"
  } else if (grepl(pattern = "Stevens", x = basename(vcf_file))) {
    reference_genome <- "Stevens"
  } else {
    reference_genome <- "Oxyccocos"
  }

  # Which family?
  if (grepl(pattern = s0_names[1], x = basename(vcf_file))) {
    s0 <- s0_names[1]
  } else {
    s0 <- s0_names[2]
  }

  vcf_in_family <- read.vcfR(file = vcf_file)
  vcf_in_family


  f1 <- s0
  parents <- unique(unlist(subset(pop_metadata, S0_name == s0, c(parent1, parent2))))
  parent1 <- parents[1]
  parent2 <- parents[2]

  # split pop metadata by the F1 individual
  family <- pop_metadata %>%
    subset(S0_name == s0) %>%
    select(individual, parent1, parent2, S0_name) %>%
    unlist() %>%
    unique()

  family_gt <- extract.gt(x = vcf_in_family, element = "GT")
  # Convert phased to unphased gt
  family_gt <- sub(pattern = "\\|", replacement = "/", x = family_gt)
  # Extract the GT data for the parents and the F1
  parents_f1 <- family_gt[,c(parents, f1)]



  ####
  #### Impute ####
  ####

  parents_f1_1 <- parents_f1
  # Any F1 genotype is a het if parent1 is 0/0 and parent2 is 1/1 or vice versa
  parent1_00_parent2_11 <- parents_f1_1[,parent1] == "0/0" & parents_f1_1[,parent2] == "1/1"
  sum(parent1_00_parent2_11, na.rm = TRUE)
  parents_f1_1[parent1_00_parent2_11, f1] <- "0/1"
  dim(parents_f1_1)

  # Use homozygous alternate for stevens if the reference genome is not stevens
  if (reference_genome != "Stevens") {
    parent1_11_parent2_00 <- parents_f1_1[,parent1] == "1/1" & parents_f1_1[,parent2] == "0/0"
    # sum(parent1_11_parent2_00, na.rm = TRUE)
    parents_f1_1[parent1_11_parent2_00, f1] <- "0/1"
    # dim(parents_f1_1)
  } else if (reference_genome == "Stevens" & "STEVENS" %in% parents) {
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
  parent1_hom <- paste0(parent1, "_HOM")
  parent2_hom <- paste0(parent2, "_HOM")

  parents_f1_4 <- cbind(parents_f1_3, hom1 = as.character(NA), hom2 = as.character(NA))
  colnames(parents_f1_4)[-1:-3] <- c(parent1_hom, parent2_hom)
  parents_f1_4[parents_f1_4[,parent1] == "0/0", parent1_hom] <- "0/0"
  parents_f1_4[parents_f1_4[,parent2] == "1/1", parent2_hom] <- "1/1"
  parents_f1_4[parents_f1_4[,parent1] == "1/1", parent1_hom] <- "1/1"
  parents_f1_4[parents_f1_4[,parent2] == "0/0", parent2_hom] <- "0/0"

  parents_f1_4[parents_f1_4[,parent1] == "0/1" & parents_f1_4[,parent2] == "1/1", parent1_hom] <- "0/0"
  parents_f1_4[parents_f1_4[,parent1] == "0/1" & parents_f1_4[,parent2] == "0/0", parent1_hom] <- "1/1"

  parents_f1_4[parents_f1_4[,parent1] == "1/1" & parents_f1_4[,parent2] == "0/1", parent2_hom] <- "0/0"
  parents_f1_4[parents_f1_4[,parent1] == "0/0" & parents_f1_4[,parent2] == "0/1", parent2_hom] <- "1/1"

  # Remove any markers that are not parent1 == 0/0 and parent2 == 1/1
  idx <- union(which(parents_f1_4[,parent1_hom] == "0/0" & parents_f1_4[,parent2_hom] == "1/1"),
               which(parents_f1_4[,parent1_hom] == "1/1" & parents_f1_4[,parent2_hom] == "0/0"))
  parents_f1_5 <- parents_f1_4[idx, ]
  # Reorder
  parents_f1_5 <- parents_f1_5[order(row.names(parents_f1_5)), ]
  dim(parents_f1_5)

  # Create a data.frame for import
  progeny_gt <- family_gt[row.names(parents_f1_5), setdiff(colnames(family_gt), colnames(parents_f1_5))]
  family_gt1 <- cbind(progeny_gt, parents_f1_5)
  # Remove STEVENS and NJ96-20
  family_gt2 <- family_gt1[, setdiff(colnames(family_gt1), c(parent1, parent2))]

  # Unique marker names
  mars <- row.names(family_gt2)
  # Ordered GT
  family_gt2 <- family_gt2[mars, ]
  family_fix2 <- vcf_in_family@fix[match(x = mars, table = vcf_in_family@fix[,"ID"]), ]

  vcf_in_family1 <- vcf_in_family
  vcf_in_family1@gt <- cbind(FORMAT = "GT", family_gt2)
  vcf_in_family1@fix <- family_fix2
  vcf_in_family1

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
                                parent1 = parent1_hom, parent2 = parent2_hom, f1 = f1, verbose = TRUE)

  # Output files for r/qtl
  geno_mat <- onemap_in$geno
  geno_mat[geno_mat == 0] <- NA

  # # Get frequencies of B allele
  # geno_mat_b_freq <- colMeans(geno_mat - 1, na.rm = TRUE) / 2
  # plot(geno_mat_b_freq)

  # Assign Macro vs Oxy alleles
  if (parent1 %in% c("BEN_LEAR", "STEVENS")) {
    macro_num <- 1
    oxy_num <- 3

  } else {
    macro_num <- 3
    oxy_num <- 1

  }
  macro_allele_code <- "M"
  oxy_allele_code <- "O"

  geno_mat[geno_mat == macro_num] <- macro_allele_code
  geno_mat[geno_mat == 2] <- "H"
  geno_mat[geno_mat == oxy_num] <- oxy_allele_code

  qtl_geno <- as.data.frame(geno_mat)
  qtl_geno <- rbind(as.integer(sub(pattern = "chr", replacement = "", x = onemap_in$CHROM)), onemap_in$POS / 1e6, qtl_geno)
  qtl_geno <- cbind(id = row.names(qtl_geno), qtl_geno)
  qtl_geno$id[1:2] <- ""
  # Save
  write.csv(x = qtl_geno, file = paste0("data/mxo_rqtl_geno_F1-", f1, "_", reference_genome, ".csv"), quote = FALSE, row.names = FALSE)

  qtl_pheno <- data.frame(id = row.names(geno_mat), fake = rnorm(length(row.names(geno_mat))))
  write.csv(x = qtl_pheno, file = paste0("data/mxo_rqtl_pheno_F1-", f1, "_", reference_genome, ".csv"), quote = FALSE, row.names = FALSE)


}







