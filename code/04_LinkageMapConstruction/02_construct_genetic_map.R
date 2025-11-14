# MxO Genetic Map
#
# Construct a linkage map
#
# This script uses R/qtl to construct a interspecific linkage map
#

library(tidyverse)
library(qtl)
# library(ASMap)
library(vcfR)

pop_metadata <- read_csv("data/population_metadata.csv")

# Colors
col <- c("slateblue", "violetred", "green3")

# # Load the VCF of the data
# vcf_in <- read.vcfR(file = "data/cran1571_BenLearv1_snps_alias.vcf.gz")

max_mar_miss <- 0.05
max_ind_miss <- 0.10
# Percent similarity to call genotypes duplicates
min_geno_dup_sim <- 0.95

# Construct the map for family 1 ------------------------------------------

## Load data into qtl ------------------------------------------------------

# List the QTL files
qtl_files <- list.files(path = "data", pattern = "rqtl", full.names = FALSE)

# First start with the larger F2 family
cross1_files <- grep(pattern = "CNJ98-325-33", x = qtl_files, value = TRUE)

# Read into qtl
cross1 <- read.cross(format = "csvs", dir = "data", genfile = cross1_files[1],
                     phefile = cross1_files[2], crosstype = "f2", estimate.map = FALSE)

summary(cross1)

# The map is in Mbp; pull it
pmap <- pull.map(cross1)
pmap_df <- map2table(map = pmap) %>%
  mutate(pos = pos * 1e6) %>%
  rownames_to_column("marker")


## Data filtering ----------------------------------------------------------

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


g_geno_freq <- cross4_gt_table1 %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_hline(yintercept = 0.50, lty = 1) +
  geom_hline(yintercept = 0.25, lty = 2) +
  geom_smooth(aes(y = pStevensHom, color = "P(AA) - Stevens"), method = "loess", span = 0.05, se = FALSE) +
  geom_smooth(aes(y = pOxyHom, color = "P(BB) - Oxycoccos"), method = "loess", span = 0.05, se = FALSE) +
  geom_smooth(aes(y = pHet, color = "P(AB) - Heterozygote"), method = "loess", span = 0.05, se = FALSE) +
  facet_wrap(~ chr, ncol = 4, scales = "free_x") +
  scale_color_manual(values = col, name = NULL) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "Genotype frequency", breaks = pretty, limits = c(0, 1)) +
  neyhart::theme_manhattan()
g_geno_freq
ggsave(filename = "output/mxo_genotype_frequency_CNJ98-325-33.jpg", plot = g_geno_freq, width = 8, height = 6, dpi = 300)





## Construct the genetic map -----------------------------------------------

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

# Run the genetic mapping function
map_cross4 <- mstmap(object = cross4, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, id = "id")





