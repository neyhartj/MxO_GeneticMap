# MxO Genetic Map
#
# Construct a linkage map
#
# This script uses R/qtl to construct a interspecific linkage map
#

library(tidyverse)
library(mappoly)

pop_metadata <- read_csv("data/population_metadata.csv")

# Colors
col <- c("slateblue", "violetred", "green3")

# # Load the VCF of the data
# vcf_in <- read.vcfR(file = "data/cran1571_BenLearv1_snps_alias.vcf.gz")

max_mar_miss <- 0.05
max_ind_miss <- 0.10
# Percent similarity to call genotypes duplicates
min_geno_dup_sim <- 0.95


# Use mappoly -------------------------------------------------------------

## Family 1 ----------------------------------------------------------------

### Convert the QTL files for mappoly ---------------------------------------

# List the QTL files
qtl_files <- list.files(path = "data", pattern = "rqtl", full.names = TRUE)

# First start with the larger F2 family
cross1_file <- grep(pattern = "geno_F1-CNJ98-325-33_BenLear", x = qtl_files, value = TRUE)

# Read the CSV file
cross1_geno_df <- read.csv(file = cross1_file)

# Extract the genotype matrix
geno_mat <- cross1_geno_df[-1:-2, -1]
# Add row names
row.names(geno_mat) <- cross1_geno_df$id[-1:-2]
# Convert to a matrix
geno_mat <- as.matrix(geno_mat)
# Convert to numeric
geno_mat[geno_mat == "B"] <- 0
geno_mat[geno_mat == "H"] <- 1
geno_mat[geno_mat == "A"] <- 2
geno_mat <- apply(geno_mat, 1, as.numeric)
row.names(geno_mat) <- colnames(cross1_geno_df)[-1]

# Format for mappoly
geno_mappoly <- data.frame(marker = row.names(geno_mat),
                           STEVENS = 1, `NJ96-20` = 1,
                           chrom = paste0("chr_", str_pad(string = unlist(cross1_geno_df[1,-1]), width = 2, side = "left", pad = "0")),
                           position = as.numeric(unlist(cross1_geno_df[2,-1])) * 1e6,
                           row.names = NULL)
geno_mappoly <- cbind(geno_mappoly, geno_mat)

# Save
write_csv(x = geno_mappoly, file = file.path(data_dir, "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear.csv"))



### Read data into mappoly --------------------------------------------------

dat <- read_geno_csv(file.in = file.path(data_dir, "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear.csv"), ploidy = 2)

# Filter on missing marker data
dat2 <- filter_missing(dat, type = "marker", filter.thres = max_mar_miss, inter = FALSE)

# Filter on missing individual data
dat3 <- filter_missing(input.data = dat2, type = "individual",
                       filter.thres = max_ind_miss, inter = FALSE)

dat3_seg_dist <- data.frame(marker = names(dat3$chisq.pval), chrom = dat3$chrom, position = dat3$genome.pos,
                            chisq_pvalue = dat3$chisq.pval, row.names = NULL)

## Remove markers that exceed the cutoff for segregation distortion unless they seem biologically relevant
# Try to do this using moving average
dat3_seg_dist_smooth <- dat3_seg_dist %>%
  split(.$chrom) %>%
  map_dfr(~{
    .x$neg_log_p <- -log10(.x$chisq_pvalue)
    fit <- loess(neg_log_p ~ position, data = .x, span = 0.2)
    pred <- predict(object = fit, newdata = .x, se = TRUE)
    .x$fitted <- pred$fit
    .x$se <- pred$se.fit
    .x
  })



marker_keep <- subset(dat3_seg_dist_smooth, neg_log_p <= fitted + 12, marker, drop = TRUE)
markers_drop <- setdiff(dat3_seg_dist_smooth$marker, marker_keep)
chisq.pval.thresh <- 10^-max(dat3_seg_dist_smooth$fitted + 12)

# Create the object to filter markers
mrks_chi_filt <- list(keep = marker_keep, exclude = marker_drop, chisq.pval.thres = chisq.pval.thresh,
                      data.name = "dat3")
class(mrks_chi_filt) <- "mappoly.chitest.seq"


# Select markers
seq_filt <- make_seq_mappoly(mrks_chi_filt)
seq_filt

# Filter redundant markers
seq_red <- elim_redundant(input.seq = seq_filt)
seq_red

# Create the sequence of ordered markers
seq_init <- make_seq_mappoly(seq_red)
plot(seq_init)

# Run two-point analysis
n_cores = parallel::detectCores() - 1
tpt <- est_pairwise_rf(seq_init, ncpus = n_cores)

# Save this
save("tpt", file = file.path(results_dir, "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear_tpt.RData"))



