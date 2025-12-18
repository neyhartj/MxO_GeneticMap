# MxO Genetic Map
#
# Construct a linkage map
#
# This script uses R/qtl to construct a interspecific linkage map
#

library(tidyverse)
library(qtl)
library(vcfR)
library(GenomicBreeding)
# library(mappoly)

pop_metadata <- read_csv("data/population_metadata.csv")

# Colors
col <- c("slateblue", "violetred", "green3")

# # Load the VCF of the data
# vcf_in <- read.vcfR(file = "data/cran1571_BenLearv1_snps_alias.vcf.gz")

max_mar_miss <- 0.05
max_ind_miss <- 0.10
max_snp_r2 <- 0.95^2
# Percent similarity to call genotypes duplicates
min_geno_dup_sim <- 0.98
min_maf <- 0.05


# List the QTL files
qtl_files <- list.files(path = "data", pattern = "rqtl", full.names = TRUE)
qtl_files <- subset(qtl_files, grepl(pattern = ".csv", x = qtl_files))



# Ad hoc functions --------------------------------------------------------

# A function for interpolating genetic positions
interpolate_genetic_pos <- function(target, genetic.map) {
  # x: df of all markers to include, chr, pos
  # genetic.map: df of marker, chr, pos, and cM

  # Order
  target <- target[order(target$chr, target$pos), ]

  # unique chromosomes
  chrs <- unique(genetic.map$chr)

  # Iterate over chromosomes
  inter_map_list <- list()
  for (i in seq_along(chrs)) {
    chr_i <- chrs[[i]]
    # Subset maps
    target_i <- target[target$chr == chr_i, ]
    genetic.map_i <- genetic.map[genetic.map$chr == chr_i, ]

    # Fill in the genetic map for 0
    if (min(genetic.map_i$pos != 0)) {
      genetic.map_i <- rbind(data.frame(marker = "dummy", chr = chr_i, cM = 0, pos = 0), genetic.map_i)
    }

    # Run interpolation with the genetic and physical map
    target_cm <- approx(
      x = genetic.map_i$pos,
      y = genetic.map_i$cM,
      xout = target_i$pos,
      rule = 2
    )$y

    # Add cm to the target df
    target_i$cM <- target_cm
    inter_map_list[[chr_i]] <- target_i
  }

  # Rbind
  inter_map <- do.call(rbind, inter_map_list)
  row.names(inter_map) <- NULL
  # Return the map
  return(inter_map)

}

aggregate_matrix <- function(M, fact, idx = NULL){
  if (is.null(idx)) {
    idx <- list(seq(1, ncol(M)))
  }
  Rlist <- lapply(idx, function(id) {
    iid <- seq(1, length(id), by = fact)
    iid <- cbind(iid, c(iid[-1]-1, length(id)))
    R <- matrix(NA, nrow(iid), nrow(iid))
    for(i in 1:(nrow(iid)-1)){
      for(j in (i+1):nrow(iid)){
        R[j,i] <-  R[i,j] <- mean(M[iid[i,1]:iid[i,2], iid[j,1]:iid[j,2]], na.rm = TRUE)
      }
    }
    R
  })
  if (length(Rlist) == 1) {
    return(Rlist[[1]])
  } else {
    Rlist
  }
}




# Family 1 / BenLear genome -----------------------------------------------

prefix <- "F1-CNJ98-325-33_BenLear"


## Use R/qtl ------------------------------------------

### Load data into qtl ------------------------------------------------------

# First start with the larger F2 family
cross1_files <- basename(grep(pattern = "CNJ98-325-33_BenLear", x = qtl_files, value = TRUE))

# Read into qtl
cross1 <- read.cross(format = "csvs", dir = "data", genfile = cross1_files[1],
                     phefile = cross1_files[2], crosstype = "f2", estimate.map = FALSE,
                     genotypes = c("M", "H", "O", "D", "C"), alleles = c("M", "O"))

summary(cross1)

# The map is in Mbp; pull it
pmap <- pull.map(cross1)
pmap_df <- map2table(map = pmap) %>%
  mutate(pos = pos * 1e6) %>%
  rownames_to_column("marker")


### Data filtering ----------------------------------------------------------

# Missing data
plotMissing(x = cross1)

cross2 <- cross1

# Filter markers by missingness
mar_missing <- nmissing(cross = cross2, what = "mar") / nind(cross2)
mar_drop <- markernames(cross2)[mar_missing > max_mar_miss]
cross2 <- drop.markers(cross = cross2, markers = mar_drop)
cross2

# Filter individuals by missingness
ind_missing <- nmissing(cross = cross2, what = "ind") / sum(nmar(cross2))
ind_keep <- ind_missing <= max_ind_miss
cross2 <- subset(x = cross2, ind = ind_keep)
cross2

# Remove markers with low MAF
geno_mat <- pull.geno(cross2) - 1
# compute MAF
af <- colMeans(geno_mat, na.rm = TRUE) / 2
plot(af)
maf <- pmin(af, 1 - af)
hist(maf)

# Remove markers with < 0.05 MAF
(min_maf * nrow(geno_mat))
mar_drop <- colnames(geno_mat)[maf < min_maf]
cross2 <- drop.markers(cross = cross2, markers = mar_drop)

# Compute marker genotype frequencies
gt_table <- geno.table(cross = cross2)

# calculate tests for seg distortion
cross2_seg_dist <- gt_table %>%
  rownames_to_column("marker") %>%
  left_join(., pmap_df) %>%
  rename(chisq_pvalue = P.value) %>%
  mutate(chr = fct_inseq(chr))

# Plot segregation distortion
chrom_colors <- ggsci::pal_aaas()(2)
chrom_colors <- rep(chrom_colors, length.out = 12)

g_seg_dist <- cross2_seg_dist %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(aes(y = -log10(chisq_pvalue), color = chr), size = 1) +
  geom_hline(yintercept = -log10(10^-15)) +
  facet_grid(~ chr, scales = "free_x", switch = "x") +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")
g_seg_dist
ggsave(filename = paste0("mxo_rqtl_", prefix, "_segregation_distortion.jpg"), plot = g_seg_dist, path = fig_dir,
       width = 10, height = 5, dpi = 300)


# Plot genotype frequencies
gt_table1 <- gt_table %>%
  mutate(neg_log_pvalue_seg_dist = -log10(P.value)) %>%
  mutate(nGenotypes = MM + OO + MO,
         nAlleles = nGenotypes * 2,
         pStevens = ((2 * MM) + MO) / nAlleles,
         pOxy = ((2 * OO) + MO) / nAlleles,
         pStevensHom = MM / nGenotypes,
         pHet = MO / nGenotypes,
         pOxyHom = OO / nGenotypes) %>%
  mutate(chr = fct_inseq(chr)) %>%
  # Add map information
  rownames_to_column("marker") %>%
  left_join(., pmap_df)

# h lines for expected frequencies
hline_df <- data.frame(type = fct_inorder(c("Expected P(Hom)", "Expected P(Het)")), y = c(0.25, 0.50))

smooth_span <- 0.1
g_geno_freq <- gt_table1 %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_hline(data = hline_df, aes(yintercept = y, lty = type)) +
  geom_smooth(aes(y = pStevensHom, color = "P(MM) - Stevens"), method = "loess", span = smooth_span, se = FALSE) +
  geom_smooth(aes(y = pOxyHom, color = "P(OO) - Oxycoccos"), method = "loess", span = smooth_span, se = FALSE) +
  geom_smooth(aes(y = pHet, color = "P(MO) - Heterozygote"), method = "loess", span = smooth_span, se = FALSE) +
  facet_wrap(~ chr, ncol = 4, scales = "free_x") +
  scale_color_manual(values = col, name = NULL) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "Genotype frequency", breaks = pretty, limits = c(0, 1)) +
  scale_linetype_discrete(name = NULL) +
  neyhart::theme_manhattan()
g_geno_freq
ggsave(filename = paste0("mxo_genotype_freq_", prefix, ".jpg"), plot = g_geno_freq, width = 8, height = 6, dpi = 300, path = fig_dir)




## Remove markers that exceed the cutoff for segregation distortion unless they seem biologically relevant
# Try to do this using smoothing
cross2_seg_dist_smooth <- cross2_seg_dist %>%
  split(.$chr) %>%
  map_dfr(~{
    .x$neg_log_p <- -log10(.x$chisq_pvalue)
    fit <- loess(neg_log_p ~ pos, data = .x, span = 0.15)
    pred <- predict(object = fit, newdata = .x, se = TRUE)
    .x$fitted <- pred$fit
    .x$se <- pred$se.fit
    .x
  })

fit_offset <- 5

cross2_seg_dist_smooth %>%
  subset(chr == 12) %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(aes(y = neg_log_p), size = 0.5, color = "grey50") +
  # geom_ribbon(aes(ymin = fitted - (3 * se), ymax = fitted + (10 * se)), fill = alpha("white", 0), color = "slateblue") +
  geom_line(aes(y = fitted + fit_offset), lwd = 1, color = "slateblue") +
  facet_wrap(~ chr, scales = "free_x", nrow = 3) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")



marker_keep <- subset(cross2_seg_dist_smooth, neg_log_p <= fitted + fit_offset, marker, drop = TRUE)
markers_drop_segdist <- setdiff(cross2_seg_dist_smooth$marker, marker_keep)
chisq.pval.thresh <- 10^-max(cross2_seg_dist_smooth$fitted + fit_offset)

# Remove these markers from the cross object
cross3 <- drop.markers(cross = cross2, markers_drop_segdist)
cross3$seg.dist <- markers_drop_segdist


# Identify redundant markers
cross3_dup_mark <- findDupMarkers(cross = cross3, exact.only = FALSE)
# Unique duplicate markers
cross3_uniq_dup_mark <- unlist(cross3_dup_mark, use.names = FALSE)
length(cross3_uniq_dup_mark)
if (is.null(cross3_uniq_dup_mark)) {
  cross3_uniq_dup_mark <- character(0)
}
cross3$duplicate <- cross3_uniq_dup_mark
# Remove the duplicate markers
cross3 <- drop.markers(cross = cross3, markers = cross3_uniq_dup_mark)
cross3

## Remove highly correlated markers
# Pull out the genotype matrix
geno_mat <- pull.geno(cross = cross3) - 1
# Impute with the mean
geno_mat <- apply(X = geno_mat, MARGIN = 2, FUN = function(snp) {
  snp[is.na(snp)] <- mean(snp, na.rm = TRUE)
  snp
})


# # Compute the correlation across all markers
# all_geno_mat_r <- cor(geno_mat)
# # Aggregate and plot
# fct <- 5
# all_geno_mat_r_agg <- aggregate_matrix(M = all_geno_mat_r, fact = fct)
# diag(all_geno_mat_r_agg) <- 1
# # Get chromosome breakpoints
# mar_by_chrom <- nmar(cross3)
# chrom_idx <- NULL
# for (i in seq_along(mar_by_chrom)) {
#   if (is.null(chrom_idx)) {
#     start <- 0 + 1
#   } else {
#     start <- max(chrom_idx$idx) + 1
#   }
#   end <- start + mar_by_chrom[i] - 1
#   end <- start + ceiling(((end - start) / fct))
#   chrom_idx1_i <- data.frame(chrom = i, idx = seq(start, end))
#   chrom_idx <- rbind(chrom_idx, chrom_idx1_i)
# }
# chrom_idx$idx <- as.character(chrom_idx$idx)
#
#
# all_geno_mat_r_agg_df <- as.data.frame(all_geno_mat_r_agg) %>%
#   rownames_to_column("idx1") %>%
#   gather(idx2, corr, -idx1) %>%
#   mutate(idx2 = parse_number(idx2),
#          idx2 = as.character(idx2)) %>%
#   left_join(., chrom_idx, by = c("idx1" = "idx")) %>%
#   left_join(., chrom_idx, by = c("idx2" = "idx")) %>%
#   mutate_at(vars(idx1, idx2), fct_inseq) %>%
#   mutate(chrom.x = fct_inseq(as.character(chrom.x)),
#          chrom.y = fct_inseq(as.character(chrom.x)),
#          chrom.y = fct_rev(chrom.y))
# g_r <- all_geno_mat_r_agg_df %>%
#   # subset(chrom.x %in% 1:2 & chrom.y %in% 1:2) %>%
#   ggplot(aes(x = idx1, y = idx2, fill = corr)) +
#   geom_tile() +
#   # facet_grid(chrom.y ~ chrom.x, scales = "free") +
#   scale_fill_gradient2(low = "slateblue", mid = "white", high = "firebrick") +
#   theme(axis.ticks = element_blank(), axis.text = element_blank())
#
# ggsave(filename = paste0("mxo_", prefix, "_marker_ld.jpg"),plot = g_r, path = fig_dir, height = 8, width = 8, dpi = 300)


# Split by chromosome
chrs <- chrnames(cross3)
mar_keep_list <- list()
for (chr in chrs) {
  chr_mars <- markernames(cross = cross3, chr = chr)
  geno_mat_chr <- geno_mat[,chr_mars, drop = FALSE]
  # Compute the correlation
  geno_mat_r2 <- cor(geno_mat_chr)^2
  # Prune based on LD
  geno_mat_pruned <- prune_LD(geno_mat_chr, cor.mat = geno_mat_r2, r2.max = max_snp_r2, check.matrix = FALSE)
  # Store the markers to be kept
  mar_keep_list[[chr]] <- colnames(geno_mat_pruned)
}

# Vector of markers to keep
mar_keep <- unlist(mar_keep_list)
length(mar_keep)
# Vector of markers that were removed
ld_pruned_mars <- setdiff(markernames(cross3), mar_keep)

# Store this information in a new cross
cross4 <- drop.markers(cross = cross3, markers = ld_pruned_mars)
cross4$ld.pruned <- ld_pruned_mars
cross4


# Identify and remove potential duplicate individuals
#
cross4_cg <- comparegeno(cross4)
hist(cross4_cg[lower.tri(cross4_cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cross4_cg[lower.tri(cross4_cg)])


cross4_cg_dup <- which(cross4_cg > 0.95, arr=TRUE)
cross4_cg_dup <- cross4_cg_dup[cross4_cg_dup[,1] < cross4_cg_dup[,2],]
cross4_cg_dup

if (nrow(cross4_cg_dup) > 0) {
  cross4 <- subset(cross4, ind = -cross4_cg_dup[,2])
}

cross4







### Construct the genetic map -----------------------------------------------

# Estimate recombination fraction
cross5 <- est.rf(cross = cross4)

# Check alleles
(allele_check <- checkAlleles(cross = cross5, threshold = 3))

# ## Plot LOD score versus rf
# rf <- pull.rf(cross5)
# lod <- pull.rf(cross5, what="lod")
# plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

# Plot the recombination fraction
plotRF(x = cross5)


# Switch alleles at markers
cross6 <- switchAlleles(cross = cross5, markers = as.character(allele_check$marker))
# Re-run est.RF
cross6 <- est.rf(cross = cross6)

# re-plot the recombination fraction
plotRF(x = cross6)

(allele_check <- checkAlleles(cross = cross6, threshold = 3))


# Look at double XOs
cross6_xo <- countXO(cross = cross6)
summary(cross6_xo)
plot(cross6_xo)

# Remove troublesome individuals and re-estimate RF
cross7 <- subset(cross6, ind = (cross6_xo < 250))
c(cross7 = nind(cross7), cross6 = nind(cross6))
# Re-estimate recombination fractions
cross7 <- est.rf(cross = cross7)
plotRF(cross7)

# Look at overall error rates and log-likelihood
loglik <- err <- seq(0.01, 0.10, by = 0.01)
for (i in seq_along(err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(cross7, error.prob=err[i], map.function = "kosambi", verbose = TRUE, n.cluster = 12)
  loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)

# Plot error rate and loglikelihood
plot(err, lod, xlab="Genotyping error rate",
     ylab=expression(paste(log[10], " likelihood")))
# Which error rate maximizes the MLE?
(error_rate_opt <- err[which.max(lod)])


# Estimate the map based on the physical order
cross7_map <- est.map(cross = cross7, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE, n.cluster = 12)

# New cross object
cross8 <- cross7

plotMap(x = cross7_map)


# Maximum cM length of a map
max_cm <- 150
# Minimum LDiff to drop a marker
min_ldiff <- 15


# Several issues; try running dropone to find problematic markers
# Start with chromosome 8
dropone <- droponemarker(cross = cross8, chr = 8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Now look at chromosome 3
dropone <- droponemarker(cross = cross8, chr = 3, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 3, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Now look at chromosome 5
dropone <- droponemarker(cross = cross8, chr = 5, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 5, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Now look at chromosome 12, but more cautiously
dropone <- droponemarker(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Run 12 again
dropone <- droponemarker(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= 5])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Run 12 a third time
dropone <- droponemarker(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= 5])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Re-estimate the whole map
new_map <- est.map(cross = cross8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE, n.cluster = 12)
plot(new_map)

# Identify all the markers that were dropped and store them
dropped_mar <- setdiff(markernames(cross7), markernames(cross8))
cross8$duplicate <- cross5$duplicate
cross8$seg.dist <- cross5$seg.dist
cross8$ld.pruned <- cross5$ld.pruned
cross8$dropped.mar <- dropped_mar


## Add the map to the cross object
cross9 <- replace.map(cross = cross8, map = new_map)
# summarize the map
summaryMap(object = cross9)
# Export the map
map_toprint <- pull.map(cross = cross9, as.table = TRUE) %>%
  rownames_to_column("marker") %>%
  rename(cM = pos) %>%
  left_join(., pmap_df)

# Save
write_csv(x = map_toprint, file = file.path(results_dir, paste0("mxo_rqtl_", prefix, "_base_map.csv")))




### Add markers that were previously removed

# Vector of markers to be added back in to the map
# markers_to_add <- c(cross9$duplicate, cross9$ld.pruned)
# markers_to_add <- cross9$duplicate
markers_to_add <- cross9$ld.pruned
markers_to_add <- sort(unique(markers_to_add))
length(markers_to_add)

# Vector of all markers in the map
all_markers <- sort(unique(c(map_toprint$marker, markers_to_add)))
length(all_markers)

# A data frame of markers to be interpolated
all_mar <- pmap_df %>%
  subset(marker %in% all_markers)

# Run interpolation
cross_interp_map <- interpolate_genetic_pos(target = all_mar, genetic.map = map_toprint)


# Plot phyical vs genetic
cross_interp_map %>%
  ggplot(aes(x = pos, y = cM)) +
  geom_point() +
  facet_wrap(~ chr)



# Save
write_csv(x = cross_interp_map, file = file.path(results_dir, paste0("mxo_rqtl_", prefix, "_interpolated_map.csv")))





### Augment the cross object and export to qtl2 --------------------------------

# Create a new cross object with this new map
cross10_map <- cross_interp_map %>%
  as.data.frame() %>%
  column_to_rownames("marker") %>%
  select(chr, cM) %>%
  table2map()
# Remove markers with duplicated positions
cross10_map1 <- map2table(cross10_map)
cross10_map1 <- table2map(cross10_map1[!duplicated(cross10_map1), ])

cross10 <- drop.markers(cross = cross2, markers = setdiff(markernames(cross2), row.names(map2table(cross10_map1))))
cross10 <- replace.map(cross = cross10, map = cross10_map1)
cross10 <- subset(x = cross10, ind = as.character(cross9$pheno$id))
cross10

# Convert to qtl2
cross_qtl <- cross10
cross_qtl2 <- qtl2::convert2cross2(cross = cross10)
cross_qtl2$pmap <- cross_interp_map %>%
  subset(marker %in% markernames(cross_qtl)) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("marker") %>%
  select(chr, pos) %>%
  table2map()

# Save
save("cross_qtl", "cross_qtl2", file = file.path(data_dir, paste0("mxo_rqtl2_", prefix, "_cross_object.RData")))



# Family 2 / BenLear genome -----------------------------------------------

prefix <- "F1-CNJ98-309-19_BenLear"


# First start with the larger F2 family
cross1_files <- basename(grep(pattern = prefix, x = qtl_files, value = TRUE))

# Read into qtl
cross1 <- read.cross(format = "csvs", dir = "data", genfile = cross1_files[1],
                     phefile = cross1_files[2], crosstype = "f2", estimate.map = FALSE,
                     genotypes = c("M", "H", "O", "D", "C"), alleles = c("M", "O"))
summary(cross1)

# The map is in Mbp; pull it
pmap <- pull.map(cross1)
pmap_df <- map2table(map = pmap) %>%
  mutate(pos = pos * 1e6) %>%
  rownames_to_column("marker")


### Data filtering ----------------------------------------------------------

# Missing data
plotMissing(x = cross1)

cross2 <- cross1

# Filter markers by missingness
mar_missing <- nmissing(cross = cross2, what = "mar") / nind(cross2)
mar_drop <- markernames(cross2)[mar_missing > max_mar_miss]
cross2 <- drop.markers(cross = cross2, markers = mar_drop)
cross2

# Filter individuals by missingness
ind_missing <- nmissing(cross = cross2, what = "ind") / sum(nmar(cross2))
ind_keep <- ind_missing <= max_ind_miss
cross2 <- subset(x = cross2, ind = ind_keep)
cross2

# Remove markers with low MAF
geno_mat <- pull.geno(cross2) - 1
# compute MAF
af <- colMeans(geno_mat, na.rm = TRUE) / 2
plot(af)
maf <- pmin(af, 1 - af)
hist(maf)

# Remove markers with < 0.05 MAF
(min_maf * nrow(geno_mat))
mar_drop <- colnames(geno_mat)[maf < min_maf]
cross2 <- drop.markers(cross = cross2, markers = mar_drop)

# Compute marker genotype frequencies
gt_table <- geno.table(cross = cross2)

# calculate tests for seg distortion
cross2_seg_dist <- gt_table %>%
  rownames_to_column("marker") %>%
  left_join(., pmap_df) %>%
  rename(chisq_pvalue = P.value) %>%
  mutate(chr = fct_inseq(chr))

# Plot segregation distortion
chrom_colors <- ggsci::pal_aaas()(2)
chrom_colors <- rep(chrom_colors, length.out = 12)

ybreaks <- pretty(range(-log10(cross2_seg_dist$chisq_pvalue)))

g_seg_dist <- cross2_seg_dist %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(aes(y = -log10(chisq_pvalue), color = chr), size = 1) +
  geom_hline(yintercept = -log10(10^-15)) +
  facet_grid(~ chr, scales = "free_x", switch = "x") +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = ybreaks, limits = range(ybreaks),  guide = guide_axis(cap = TRUE)) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")
g_seg_dist
ggsave(filename = paste0("mxo_rqtl_", prefix, "_segregation_distortion.jpg"), plot = g_seg_dist, path = fig_dir,
       width = 10, height = 5, dpi = 300)


# Plot genotype frequencies
gt_table1 <- gt_table %>%
  mutate(neg_log_pvalue_seg_dist = -log10(P.value)) %>%
  mutate(nGenotypes = MM + OO + MO,
         nAlleles = nGenotypes * 2,
         pStevens = ((2 * MM) + MO) / nAlleles,
         pOxy = ((2 * OO) + MO) / nAlleles,
         pStevensHom = MM / nGenotypes,
         pHet = MO / nGenotypes,
         pOxyHom = OO / nGenotypes) %>%
  mutate(chr = fct_inseq(chr)) %>%
  # Add map information
  rownames_to_column("marker") %>%
  left_join(., pmap_df)

# h lines for expected frequencies
hline_df <- data.frame(type = fct_inorder(c("Expected P(Hom)", "Expected P(Het)")), y = c(0.25, 0.5))

smooth_span <- 0.1
g_geno_freq <- gt_table1 %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_hline(data = hline_df, aes(yintercept = y, lty = type)) +
  geom_smooth(aes(y = pStevensHom, color = "P(MM) - Ben Lear"), method = "loess", span = smooth_span, se = FALSE) +
  geom_smooth(aes(y = pOxyHom, color = "P(OO) - Oxycoccos"), method = "loess", span = smooth_span, se = FALSE) +
  geom_smooth(aes(y = pHet, color = "P(MO) - Heterozygote"), method = "loess", span = smooth_span, se = FALSE) +
  facet_wrap(~ chr, ncol = 4, scales = "free_x") +
  scale_color_manual(values = col, name = NULL) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "Genotype frequency", breaks = pretty, limits = c(0, 1)) +
  scale_linetype_discrete(name = NULL) +
  neyhart::theme_manhattan()
g_geno_freq
ggsave(filename = paste0("mxo_genotype_freq_", prefix, ".jpg"), plot = g_geno_freq, width = 8, height = 6, dpi = 300, path = fig_dir)




## Remove markers that exceed the cutoff for segregation distortion unless they seem biologically relevant
# Try to do this using smoothing
cross2_seg_dist_smooth <- cross2_seg_dist %>%
  split(.$chr) %>%
  map_dfr(~{
    .x$neg_log_p <- -log10(.x$chisq_pvalue)
    fit <- loess(neg_log_p ~ pos, data = .x, span = 0.15)
    pred <- predict(object = fit, newdata = .x, se = TRUE)
    .x$fitted <- pred$fit
    .x$se <- pred$se.fit
    .x
  })

fit_offset <- 1

cross2_seg_dist_smooth %>%
  subset(chr == 12) %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(aes(y = neg_log_p), size = 0.5, color = "grey50") +
  # geom_ribbon(aes(ymin = fitted - (3 * se), ymax = fitted + (10 * se)), fill = alpha("white", 0), color = "slateblue") +
  geom_line(aes(y = fitted + fit_offset), lwd = 1, color = "slateblue") +
  facet_wrap(~ chr, scales = "free_x", nrow = 3) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")



marker_keep <- subset(cross2_seg_dist_smooth, neg_log_p <= fitted + fit_offset, marker, drop = TRUE)
markers_drop_segdist <- setdiff(cross2_seg_dist_smooth$marker, marker_keep)
chisq.pval.thresh <- 10^-max(cross2_seg_dist_smooth$fitted + fit_offset)

# Remove these markers from the cross object
cross3 <- drop.markers(cross = cross2, markers_drop_segdist)
cross3$seg.dist <- markers_drop_segdist


# Identify redundant markers
cross3_dup_mark <- findDupMarkers(cross = cross3, exact.only = FALSE)
# Unique duplicate markers
cross3_uniq_dup_mark <- unlist(cross3_dup_mark, use.names = FALSE)
length(cross3_uniq_dup_mark)
if (is.null(cross3_uniq_dup_mark)) {
  cross3_uniq_dup_mark <- character(0)
}
cross3$duplicate <- cross3_uniq_dup_mark
# Remove the duplicate markers
cross3 <- drop.markers(cross = cross3, markers = cross3_uniq_dup_mark)
cross3

## Remove highly correlated markers
# Pull out the genotype matrix
geno_mat <- pull.geno(cross = cross3) - 1
# Impute with the mean
geno_mat <- apply(X = geno_mat, MARGIN = 2, FUN = function(snp) {
  snp[is.na(snp)] <- mean(snp, na.rm = TRUE)
  snp
})

# Split by chromosome
chrs <- chrnames(cross3)
mar_keep_list <- list()
for (chr in chrs) {
  chr_mars <- markernames(cross = cross3, chr = chr)
  geno_mat_chr <- geno_mat[,chr_mars, drop = FALSE]
  # Compute the correlation
  geno_mat_r2 <- cor(geno_mat_chr)^2
  # Prune based on LD
  geno_mat_pruned <- prune_LD(geno_mat_chr, cor.mat = geno_mat_r2, r2.max = max_snp_r2, check.matrix = FALSE)
  # Store the markers to be kept
  mar_keep_list[[chr]] <- colnames(geno_mat_pruned)

}

# Vector of markers to keep
mar_keep <- unlist(mar_keep_list)
# Vector of markers that were removed
ld_pruned_mars <- setdiff(markernames(cross3), mar_keep)

# Store this information in a new cross
cross4 <- drop.markers(cross = cross3, markers = ld_pruned_mars)
cross4$ld.pruned <- ld_pruned_mars
cross4


# Identify and remove potential duplicate individuals
#
cross4_cg <- comparegeno(cross4)
hist(cross4_cg[lower.tri(cross4_cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cross4_cg[lower.tri(cross4_cg)])


cross4_cg_dup <- which(cross4_cg > min_geno_dup_sim, arr=TRUE)
cross4_cg_dup <- cross4_cg_dup[cross4_cg_dup[,1] < cross4_cg_dup[,2],]
cross4_cg_dup

if (nrow(cross4_cg_dup) > 0) {
  cross4 <- subset(cross4, ind = -cross4_cg_dup[,2])
}

cross4







### Construct the genetic map -----------------------------------------------

# Estimate recombination fraction
cross5 <- est.rf(cross = cross4)

# Check alleles
(allele_check <- checkAlleles(cross = cross5, threshold = 3))

# ## Plot LOD score versus rf
# rf <- pull.rf(cross5)
# lod <- pull.rf(cross5, what="lod")
# plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

# Plot the recombination fraction
jpeg(filename = file.path(fig_dir, paste0("mxo_", prefix, "_recombination_matrix.jpg")), width = 8, height = 8,
     units = "in", res = 300)
plotRF(x = cross5)
dev.off()

# Switch alleles at markers
cross6 <- switchAlleles(cross = cross5, markers = as.character(allele_check$marker))
# Re-run est.RF
cross6 <- est.rf(cross = cross6)

# re-plot the recombination fraction
plotRF(x = cross6)

(allele_check <- checkAlleles(cross = cross6, threshold = 3))


# Look at double XOs
cross6_xo <- countXO(cross = cross6)
summary(cross6_xo)
plot(cross6_xo)

# Remove troublesome individuals and re-estimate RF
cross7 <- subset(cross6, ind = (cross6_xo < 200))
# Re-estimate recombination fractions
cross7 <- est.rf(cross = cross7)
plotRF(cross7)

# Look at overall error rates and log-likelihood
loglik <- err <- seq(0.01, 0.10, by = 0.01)
for (i in seq_along(err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(cross7, error.prob=err[i], map.function = "kosambi", verbose = TRUE, n.cluster = 12)
  loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)

# Plot error rate and loglikelihood
plot(err, lod, xlab="Genotyping error rate",
     ylab=expression(paste(log[10], " likelihood")))
# Which error rate maximizes the MLE?
(error_rate_opt <- err[which.max(lod)])


# Estimate the map based on the physical order
cross7_map <- est.map(cross = cross7, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE, n.cluster = 12)

# New cross object
cross8 <- cross7

plotMap(x = cross7_map)


# Maximum cM length of a map
max_cm <- 150
# Minimum LDiff to drop a marker
min_ldiff <- 15

# Several issues; try running dropone to find problematic markers
# Start with chromosome 8
dropone <- droponemarker(cross = cross8, chr = 8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Chromosome 3
dropone <- droponemarker(cross = cross8, chr = 3, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 3, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Chromosome 2
dropone <- droponemarker(cross = cross8, chr = 2, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff])
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 2, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Chromosome 12
dropone <- droponemarker(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= 5])
cross8 <- drop.markers(cross8, intersect(markernames(cross8), badmar))
new_map_i <- est.map(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Chromosome 12 - run again
dropone <- droponemarker(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
(badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= 5])
cross8 <- drop.markers(cross8, intersect(markernames(cross8), badmar))
new_map_i <- est.map(cross = cross8, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Chromosome 9
dropone <- droponemarker(cross = cross8, chr = 9, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff]
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 9, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Chromosome 1
dropone <- droponemarker(cross = cross8, chr = 1, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff]
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 1, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Chromosome 5
dropone <- droponemarker(cross = cross8, chr = 5, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff]
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 5, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)


# Re-estimate the whole map
new_map <- est.map(cross = cross8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE, n.cluster = 12)
plot(new_map)



# Identify all the markers that were dropped and store them
dropped_mar <- setdiff(markernames(cross7), markernames(cross8))
dropped_mar
cross8$duplicate <- cross5$duplicate
cross8$seg.dist <- cross5$seg.dist
cross8$ld.pruned <- cross5$ld.pruned
cross8$dropped.mar <- dropped_mar


## Add the map to the cross object
cross9 <- replace.map(cross = cross8, map = new_map)
# summarize the map
summaryMap(object = cross9)
# Export the map
map_toprint <- pull.map(cross = cross9, as.table = TRUE) %>%
  rownames_to_column("marker") %>%
  rename(cM = pos) %>%
  left_join(., pmap_df)

# Save
write_csv(x = map_toprint, file = file.path(results_dir, paste0("mxo_rqtl_", prefix, "_base_map.csv")))

# Plot phyical vs genetic
map_toprint %>%
  ggplot(aes(x = pos, y = cM)) +
  geom_point() +
  facet_wrap(~ chr)



### Add markers that were previously removed

# Vector of markers to be added back in to the map
# markers_to_add <- c(cross9$duplicate, cross9$ld.pruned)
markers_to_add <- cross9$duplicate
# markers_to_add <- cross9$ld.pruned
markers_to_add <- sort(unique(markers_to_add))
length(markers_to_add)

# Vector of all markers in the map
all_markers <- sort(unique(c(map_toprint$marker, markers_to_add)))
length(all_markers)

# A data frame of markers to be interpolated
all_mar <- pmap_df %>%
  subset(marker %in% all_markers)

# Run interpolation
cross_interp_map <- interpolate_genetic_pos(target = all_mar, genetic.map = map_toprint)


# Save
write_csv(x = cross_interp_map, file = file.path(results_dir, paste0("mxo_rqtl_", prefix, "_interpolated_map.csv")))

# Plot phyical vs genetic
cross_interp_map %>%
  ggplot(aes(x = pos, y = cM)) +
  geom_point() +
  facet_wrap(~ chr)



### Augment the cross object and export to qtl2 --------------------------------

# Create a new cross object with this new map
cross10_map <- cross_interp_map %>%
  as.data.frame() %>%
  column_to_rownames("marker") %>%
  select(chr, cM) %>%
  table2map()
# Remove markers with duplicated positions
cross10_map1 <- map2table(cross10_map)
cross10_map1 <- table2map(cross10_map1[!duplicated(cross10_map1), ])

cross10 <- drop.markers(cross = cross2, markers = setdiff(markernames(cross2), row.names(map2table(cross10_map1))))
cross10 <- replace.map(cross = cross10, map = cross10_map1)
cross10 <- subset(x = cross10, ind = as.character(cross9$pheno$id))
cross10

# Convert to qtl2
library(qtl2)

cross_qtl <- cross10
cross_qtl2 <- convert2cross2(cross = cross10)
cross_qtl2$pmap <- cross_interp_map %>%
  subset(marker %in% markernames(cross_qtl)) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("marker") %>%
  select(chr, pos) %>%
  table2map()

# Save
save("cross_qtl", "cross_qtl2", file = file.path(data_dir, paste0("mxo_rqtl2_", prefix, "_cross_object.RData")))





