# MxO Genetic Map
#
# Construct a linkage map
#
# This script uses R/qtl to construct a interspecific linkage map
#

library(tidyverse)
library(qtl)
library(vcfR)
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

# Use R/qtl ------------------------------------------


## Family 1 ----------------------------------------------------------------



### Load data into qtl ------------------------------------------------------

# List the QTL files
qtl_files <- list.files(path = "data", pattern = "rqtl", full.names = FALSE)

# First start with the larger F2 family
cross1_files <- grep(pattern = "CNJ98-325-33_BenLear", x = qtl_files, value = TRUE)

# Read into qtl
cross1 <- read.cross(format = "csvs", dir = "data", genfile = cross1_files[1],
                     phefile = cross1_files[2], crosstype = "f2", estimate.map = FALSE)

summary(cross1)

# The map is in Mbp; pull it
pmap <- pull.map(cross1)
pmap_df <- map2table(map = pmap) %>%
  mutate(pos = pos * 1e6) %>%
  rownames_to_column("marker")


### Data filtering ----------------------------------------------------------

# Missing data
plotMissing(x = cross1)

cross1_1 <- cross1

# Calculate missing data proportion by marker
p_miss_mar <- 1 - (ntyped(cross1_1, "mar") / nind(cross1_1))
hist(p_miss_mar)

# Remove markers with high missing data (> 10%)
markers_drop <- names(p_miss_mar)[p_miss_mar > max_mar_miss]
length(markers_drop)
cross1_1 <- drop.markers(cross = cross1_1, markers = markers_drop)

# Calculate missing data proportion by individual
p_miss_ind <- 1 - (ntyped(cross1_1) / sum(nmar(cross1_1)))
hist(p_miss_ind)

# Remove markers with high missing data (> 10%)
ind_keep <- names(p_miss_ind)[p_miss_ind <= max_ind_miss]
nind(cross1_1) - length(ind_keep)
cross1_1 <- subset(x = cross1_1, ind = ind_keep)

summary(cross1_1)


# Potential duplicates
#
cross1_1_cg <- comparegeno(cross1_1)
hist(cross1_1_cg[lower.tri(cross1_1_cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cross1_1_cg[lower.tri(cross1_1_cg)])

cross1_1_cg_dup <- which(cross1_1_cg > min_geno_dup_sim, arr=TRUE)
cross1_1_cg_dup <- cross1_1_cg_dup[cross1_1_cg_dup[,1] < cross1_1_cg_dup[,2],]
cross1_1_cg_dup

cross3 <- subset(cross1_1, ind = -cross1_1_cg_dup[,2])
cross3

## Markers

# Remove markers that are not segregating
cross3_gt_table <- geno.table(cross = cross3)
markers_no_seg <- which(apply(X = cross3_gt_table[,3:5], MARGIN = 1, FUN = function(x) any(x == 0)))
cross4 <- drop.markers(cross = cross3, markers = names(markers_no_seg))
cross4

# calculate tests for seg distortion
cross4_gt_table <- geno.table(cross = cross4)
# Get distorted markers
(bonf_thresh <- 0.05 / sum(nmar(cross4)))
cross4_seg_dist <- rownames(cross4_gt_table[cross4_gt_table$P.value < bonf_thresh, ])
length(cross4_seg_dist)
# Set these aside
cross4$seg.distortion <- cross4_seg_dist



# Identify duplicate markers
cross4_dup_mark <- findDupMarkers(cross = cross4, exact.only = FALSE)
# Unique duplicate markers
cross4_uniq_dup_mark <- unlist(cross4_dup_mark, use.names = FALSE)
length(cross4_uniq_dup_mark)
cross4$duplicate <- cross4_uniq_dup_mark


# Drop distorted and duplicate markers
marker_drop <- union(cross4_seg_dist, cross4_uniq_dup_mark)
length(marker_drop)

cross5 <- drop.markers(cross = cross4, markers = marker_drop)
cross5



# Plot segregation distortion and genotype frequencies
cross4_gt_table1 <- cross4_gt_table %>%
  mutate(neg_log_pvalue_seg_dist = -log10(P.value)) %>%
  mutate(nGenotypes = AA + BB + AB,
         nAlleles = nGenotypes * 2,
         pStevens = ((2 * AA) + AB) / nAlleles,
         pOxy = ((2 * BB) + AB) / nAlleles,
         pStevensHom = AA / nGenotypes,
         pHet = AB / nGenotypes,
         pOxyHom = BB / nGenotypes) %>%
  mutate(chr = fct_inseq(chr)) %>%
  # Add map information
  rownames_to_column("marker") %>%
  left_join(., pmap_df)

# Plot
g_seg_dist <- cross4_gt_table1 %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(aes(y = neg_log_pvalue_seg_dist)) +
  facet_grid(~ chr, scales = "free_x") +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = pretty) +
  theme_classic()
g_seg_dist


# h lines for expected frequencies
hline_df <- data.frame(type = fct_inorder(c("Expected P(Hom)", "Expected P(Het)")), y = c(0.5, 0.25))

g_geno_freq <- cross4_gt_table1 %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_hline(data = hline_df, aes(yintercept = y, lty = type)) +
  geom_smooth(aes(y = pStevensHom, color = "P(AA) - Stevens"), method = "loess", span = 0.05, se = FALSE) +
  geom_smooth(aes(y = pOxyHom, color = "P(BB) - Oxycoccos"), method = "loess", span = 0.05, se = FALSE) +
  geom_smooth(aes(y = pHet, color = "P(AB) - Heterozygote"), method = "loess", span = 0.05, se = FALSE) +
  facet_wrap(~ chr, ncol = 4, scales = "free_x") +
  scale_color_manual(values = col, name = NULL) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "Genotype frequency", breaks = pretty, limits = c(0, 1)) +
  scale_linetype_discrete(name = NULL) +
  neyhart::theme_manhattan()
g_geno_freq
ggsave(filename = "mxo_genotype_frequency_CNJ98-325-33.jpg", plot = g_geno_freq, width = 8, height = 6, dpi = 300, path = fig_dir)





### Construct the genetic map -----------------------------------------------

# Estimate recombination fraction
cross5 <- est.rf(cross = cross5)

# Plot the recombination fraction
plotRF(x = cross5)

# Pull the RF matrix and the LOD matrix
cross5_lod <- pull.rf(cross = cross5, what = "lod")
cross5_rf <- pull.rf(cross = cross5, what = "rf")

# Histogram of adjacent marker RF
rf <- pull.rf(cross = cross5, what = "rf")

# Estimate the map based on the physical order
cross5_map <- est.map(cross = cross5, error.prob = 0.03, map.function = "kosambi", verbose = TRUE, n.cluster = 4)


# Use the ripple function to test local order and refine the marker order
cross5_rippled <- ripple(cross = cross5, chr = 12, window = 4, error.prob = 0.01, method = "likelihood")
summary(cross5_rippled)

# Form linkage groups
cross5_lg <- formLinkageGroups(cross = cross5, max.rf = 0.5, min.lod = 25)
table(cross5_lg)














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

print(dat,detailed = T)

plot(dat)

# Filter on missing marker data
dat2 <- filter_missing(dat, type = "marker", filter.thres = 0.05)

# Filter on missing individual data
dat3 <- filter_missing(input.data = dat2, type = "individual",
                       filter.thres = 0.05, inter = TRUE)

dat3_seg_dist <- data.frame(marker = names(dat3$chisq.pval), chrom = dat3$chrom, position = dat3$genome.pos,
                            chisq_pvalue = dat3$chisq.pval, row.names = NULL)

# Plot segregation distortion
chrom_colors <- ggsci::pal_aaas()(2)
chrom_colors <- rep(chrom_colors, length.out = 12)

g_seg_dist <- dat3_seg_dist %>%
  ggplot(aes(x = position / 1e6)) +
  geom_point(aes(y = -log10(chisq_pvalue), color = chrom), size = 1) +
  geom_hline(yintercept = -log10(10^-15)) +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")
g_seg_dist
ggsave(filename = "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear_segregation_distortion.jpg", plot = g_seg_dist, path = fig_dir,
       width = 10, height = 5, dpi = 300)


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


dat3_seg_dist_smooth %>%
  subset(chrom == "chr_12") %>%
  ggplot(aes(x = position / 1e6)) +
  geom_point(aes(y = neg_log_p), size = 0.5, color = "grey50") +
  # geom_ribbon(aes(ymin = fitted - (3 * se), ymax = fitted + (10 * se)), fill = alpha("white", 0), color = "slateblue") +
  geom_line(aes(y = fitted + 12), lwd = 1, color = "slateblue") +
  facet_wrap(~ chrom, scales = "free_x", switch = "x", nrow = 3) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")


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



