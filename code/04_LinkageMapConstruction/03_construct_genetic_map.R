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
max_snp_r2 <- 0.8
# Percent similarity to call genotypes duplicates
min_geno_dup_sim <- 0.98
min_maf <- 0.05


# List the QTL files
qtl_files <- list.files(path = "data", pattern = "rqtl", full.names = TRUE)



# Family 1 / BenLear genome -----------------------------------------------

prefix <- "F1-CNJ98-325-33_BenLear"


# # Filter the data for R/qtl and mappoly
#
# # First start with the larger F2 family
# cross1_file <- grep(pattern = prefix, x = qtl_files, value = TRUE)
# cross1_geno_file <- grep(pattern = "geno", x = cross1_file, value = TRUE)
#
# # Read the CSV file
# cross1_geno_df <- read.csv(file = cross1_geno_file)
#
# # Get chrom and pos information
# snp_info <- data.frame(marker = names(cross1_geno_df)[-1],
#                        chrom = paste0("chr_", str_pad(string = unlist(cross1_geno_df[1,-1]), width = 2, side = "left", pad = "0")),
#                        position = as.numeric(unlist(cross1_geno_df[2,-1])) * 1e6,
#                        row.names = NULL)
# dim(snp_info)
#
# # Extract the genotype matrix
# geno_mat <- cross1_geno_df[-1:-2, -1]
# # Add row names
# row.names(geno_mat) <- cross1_geno_df$id[-1:-2]
# # Convert to a matrix
# geno_mat <- as.matrix(geno_mat)
# # Convert to numeric
# geno_mat[geno_mat == "B"] <- 0
# geno_mat[geno_mat == "H"] <- 1
# geno_mat[geno_mat == "A"] <- 2
# geno_mat <- apply(geno_mat, 1, as.numeric)
# row.names(geno_mat) <- colnames(cross1_geno_df)[-1]
#
# # Look at the MAF
# af <- rowMeans(geno_mat, na.rm = TRUE) / 2
# plot(af)
# maf <- pmin(af, 1 - af)
# hist(maf)
# plot(maf)
#
# # Remove markers with < 0.05 MAF
# (min_maf * ncol(geno_mat))
# geno_mat1 <- geno_mat[maf >= min_maf, ]
#
# # Remove markers and individuals with > 0.05 missing data
# geno_mat1 <- geno_mat1[rowMeans(is.na(geno_mat1)) <= max_mar_miss, ]
# geno_mat1 <- geno_mat1[, colMeans(is.na(geno_mat1)) <= max_ind_miss]
#
# # Remove highly correlated markers within each chromosome
# mar_list <- split(snp_info$marker, snp_info$chrom) %>%
#   map(~intersect(., row.names(geno_mat1)))
#
# geno_mat1_list <- mar_list %>%
#   map(~{
#     geno_mat_x <- geno_mat1[.x, ] - 1
#     # Impute with the mean
#     geno_mat_x <- apply(X = geno_mat_x, MARGIN = 1, FUN = function(snp) {
#       snp[is.na(snp)] <- mean(snp, na.rm = TRUE)
#       snp
#     })
#     # Compute the correlation
#     geno_mat_x_r2 <- cor(geno_mat_x)^2
#     # Prune and return the new matrix
#     prune_LD(geno_mat_x, cor.mat = geno_mat_x_r2, r2.max = max_snp_r2, check.matrix = FALSE)
#
#   })
#
# # Re-merge
# markers_keep <- colnames(do.call(cbind, geno_mat1_list))
# length(markers_keep)
# # Subset the above geno mat for these markers
# geno_mat3 <- geno_mat1[markers_keep, ]
# idx <- match(x = markers_keep, table = snp_info$marker)
#
# # Format for mappoly
# geno_mappoly <- data.frame(marker = row.names(geno_mat3),
#                            STEVENS = 1, `NJ96-20` = 1,
#                            chrom = parse_number(snp_info$chrom[idx]),
#                            position = snp_info$position[idx],
#                            row.names = NULL)
# geno_mappoly <- cbind(geno_mappoly, geno_mat3)
#
# # Save
# write_csv(x = geno_mappoly, file = file.path(data_dir, "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear.csv"))
#
#
#
# # Format for r/QTL
# geno_mat4 <- t(geno_mat3)
# geno_mat4[geno_mat4 == 0] <- "B"
# geno_mat4[geno_mat4 == 1] <- "H"
# geno_mat4[geno_mat4 == 2] <- "A"
#
# geno_df <- as.data.frame(geno_mat4) %>%
#   rownames_to_column("id")
#
# geno_df <- rbind(c("", as.character(geno_mappoly$chrom)),
#                  c("", geno_mappoly$position / 1e6),
#                  geno_df)
#
# # Read in the phenotype CSV and subset clones
# cross1_pheno_file <- grep(pattern = "pheno", x = cross1_file, value = TRUE)
#
# # Read the CSV file
# pheno_df <- read.csv(file = cross1_pheno_file) %>%
#   subset(id %in% geno_df$id)
#
# # Save all
# write.csv(x = geno_df, file = sub(pattern = ".csv", replacement = "_filtered.csv", x = grep(pattern = "geno", x = cross1_file, value = TRUE)),
#           row.names = FALSE, quote = FALSE, col.names = TRUE)
#
# write.csv(x = pheno_df, file = sub(pattern = ".csv", replacement = "_filtered.csv", x = grep(pattern = "pheno", x = cross1_file, value = TRUE)),
#           row.names = FALSE, quote = FALSE, col.names = TRUE)
#




## Use R/qtl ------------------------------------------

### Load data into qtl ------------------------------------------------------

# First start with the larger F2 family
cross1_files <- basename(grep(pattern = "CNJ98-325-33_BenLear", x = qtl_files, value = TRUE))

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
cross2_seg_dist <- cross2_gt_table %>%
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

# h lines for expected frequencies
hline_df <- data.frame(type = fct_inorder(c("Expected P(Hom)", "Expected P(Het)")), y = c(0.5, 0.25))

smooth_span <- 0.1
g_geno_freq <- gt_table1 %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_hline(data = hline_df, aes(yintercept = y, lty = type)) +
  geom_smooth(aes(y = pStevensHom, color = "P(AA) - Stevens"), method = "loess", span = smooth_span, se = FALSE) +
  geom_smooth(aes(y = pOxyHom, color = "P(BB) - Oxycoccos"), method = "loess", span = smooth_span, se = FALSE) +
  geom_smooth(aes(y = pHet, color = "P(AB) - Heterozygote"), method = "loess", span = smooth_span, se = FALSE) +
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
# Drop the marker
summary(dropone, lodcolumn = 2)
badmar <- rownames(summary(dropone, lod.column=2))
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plot(new_map_i)


# Now look at chromosome 3
dropone <- droponemarker(cross = cross8, chr = 3, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
# Find all potentially bad markers
badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff]
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 3, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Re-estimate the whole map
new_map <- est.map(cross = cross8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plot(new_map)


# Now look at chromosome 5
dropone <- droponemarker(cross = cross8, chr = 5, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
# Find all potentially bad markers
badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff]
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross8, chr = 5, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
plotMap(new_map_i)

# Re-estimate the whole map
new_map <- est.map(cross = cross8, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE, n.cluster = 12)
plot(new_map)


# Now look at chromosome 12, but more cautiously
dropone <- droponemarker(cross = cross9, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
# Plot the difference from dropping markers
plot(dropone, lod = 2)
## It's one marker, so OK dropping it
# Find all potentially bad markers
badmar <- row.names(dropone)[sign(dropone$LOD) == 1 | dropone$Ldiff >= min_ldiff]
cross8 <- drop.markers(cross8, badmar)
new_map_i <- est.map(cross = cross9, chr = 12, error.prob = error_rate_opt, map.function = "kosambi", verbose = TRUE)
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
markers_to_add <- c(cross9$duplicate, cross9$ld.pruned)
markers_to_add <- sort(unique(markers_to_add))
length(markers_to_add)

# Vector of all markers in the map
all_markers <- sort(unique(c(map_toprint$marker, markers_to_add)))
length(all_markers)

# A data frame of markers to be interpolated
all_mar <- pmap_df %>%
  subset(marker %in% all_markers)

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


# Run interpolation
cross_interp_map <- interpolate_genetic_pos(target = all_mar, genetic.map = map_toprint)

# Save
write_csv(x = cross_interp_map, file = file.path(results_dir, paste0("mxo_rqtl_", prefix, "_interpolated_map.csv")))





### Augment the cross object and export to qtl2 --------------------------------

cross10 <- drop.markers(cross = cross2, markers = setdiff(markernames(cross2), cross_interp_map$marker))
cross10 <- subset(x = cross10, ind = as.character(cross9$pheno$id))
cross10

# Convert to qtl2
library(qtl2)

cross_qtl2 <- convert2cross2(cross = cross10)

# Save
save("cross_qtl2", file = file.path(data_dir, paste0("mxo_rqtl2_", prefix, "_cross_object.RData")))



# Use mappoly -------------------------------------------------------------

## Family 1 ----------------------------------------------------------------


### Convert the QTL files for mappoly ---------------------------------------



### Read data into mappoly --------------------------------------------------

dat <- read_geno_csv(file.in = file.path(data_dir, "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear.csv"), ploidy = 2)

print(dat,detailed = T)

plot(dat)

# Filter on missing marker data
dat2 <- filter_missing(dat, type = "marker", filter.thres = max_mar_miss, inter = FALSE)

# Filter on missing individual data
dat3 <- filter_missing(input.data = dat2, type = "individual", filter.thres = max_ind_miss, inter = FALSE)

dat3_seg_dist <- data.frame(marker = names(dat3$chisq.pval), chrom = fct_inseq(dat3$chrom), position = dat3$genome.pos,
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
    fit <- loess(neg_log_p ~ position, data = .x, span = 0.5)
    pred <- predict(object = fit, newdata = .x, se = TRUE)
    .x$fitted <- pred$fit
    .x$se <- pred$se.fit
    .x
  })


dat3_seg_dist_smooth %>%
  # subset(chrom == "chr_12") %>%
  ggplot(aes(x = position / 1e6)) +
  geom_point(aes(y = neg_log_p), size = 0.5, color = "grey50") +
  # geom_ribbon(aes(ymin = fitted - (3 * se), ymax = fitted + (10 * se)), fill = alpha("white", 0), color = "slateblue") +
  geom_line(aes(y = fitted + 12), lwd = 1, color = "slateblue") +
  facet_wrap(~ chrom, scales = "free_x", nrow = 3) +
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
mrks_chi_filt <- list(keep = marker_keep, exclude = markers_drop, chisq.pval.thres = chisq.pval.thresh,
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
length(seq_init$seq.mrk.names)

# Run two-point analysis
n_cores <- 12
tpt <- est_pairwise_rf(seq_init, ncpus = n_cores)

# Convert to a matrix
m <- rf_list_to_matrix(tpt)

# Get the genomic order
id_order <- get_genomic_order(seq_init, verbose = TRUE)
plot(id_order)

# Create a sequence in the genomic order
sgo <- make_seq_mappoly(id_order)

# Plot recombination fraction using genomic order
plot(m, ord = sgo, fact = 5)


# Group markers
grs <- group_mappoly(input.mat = m,
                     expected.groups = 12,
                     comp.mat = TRUE,
                     inter = TRUE)
grs

# Data frame of marker information
data_gz <- data.frame(marker = seq_init$seq.mrk.names, chrom = seq_init$chrom, position = seq_init$genome.pos, row.names = NULL)
dim(data_gz)

## Estimate the map using chromosome and genomic info
MAP <-vector("list", 12)
for(i in 1:12){
  markers_i <- data_gz$marker[data_gz$chrom == i]
  s <- make_seq_mappoly(dat3, markers_i)
  #genome order
  id_order <- get_genomic_order(s)
  # creates a sequence of markers in the genome order
  sgo <- make_seq_mappoly(id_order)
  #Estimating the genetic map for a given order involves the computation of recombination fraction
  map <- est_rf_hmm_sequential(sgo, twopt = tpt, extend.tail = 10, verbose = TRUE)
  # update the recombination fractions by allowing a global error in the HMM recombination fraction re-estimation
  MAP[[i]] <- est_full_hmm_with_global_error(map, error = .05)
}














# ## try mappoly2?
#
# dat <- read_geno_csv(file.in = file.path(data_dir, "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear.csv"), ploidy.p1 = 2, ploidy.p2 = 2,
#                      name.p1 = "STEVENS", name.p2 = "NJ96-20")
#
# print(dat,detailed = T)
#
# plot(dat)
#
# # Rename
# dat1 <- dat
#
# # Quality control - missingness and segregation distortion
# dat1_geno_dose <- dat1$geno.dose
# dat1_geno_dose <- dat1_geno_dose[rowMeans(is.na(dat1_geno_dose)) <= max_mar_miss, ]
# dat1_geno_dose <- dat1_geno_dose[, colMeans(is.na(dat1_geno_dose)) <= max_ind_miss]
#
# dat3_seg_dist <- dat1$QAQC.values$markers %>%
#   rownames_to_column("marker") %>%
#   left_join(., data.frame(marker = names(dat1$chrom), chrom = dat1$chrom, position = dat1$genome.pos, row.names = NULL)) %>%
#   mutate(chrom = fct_inseq(chrom))
#
# # Plot segregation distortion
# chrom_colors <- ggsci::pal_aaas()(2)
# chrom_colors <- rep(chrom_colors, length.out = 12)
#
# g_seg_dist <- dat3_seg_dist %>%
#   ggplot(aes(x = position / 1e6)) +
#   geom_point(aes(y = -log10(chisq.pval), color = chrom), size = 1) +
#   geom_hline(yintercept = -log10(10^-15)) +
#   facet_grid(~ chrom, scales = "free_x", switch = "x") +
#   scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
#   scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
#   scale_color_manual(values = chrom_colors, guide = "none") +
#   theme_classic() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
#         panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")
# g_seg_dist
# ggsave(filename = "mxo_mappoly_geno_F1-CNJ98-325-33_BenLear_segregation_distortion.jpg", plot = g_seg_dist, path = fig_dir,
#        width = 10, height = 5, dpi = 300)
#
#
# ## Remove markers that exceed the cutoff for segregation distortion unless they seem biologically relevant
# # Try to do this using moving average
# dat3_seg_dist_smooth <- dat3_seg_dist %>%
#   split(.$chrom) %>%
#   map_dfr(~{
#     .x$neg_log_p <- -log10(.x$chisq.pval)
#     fit <- loess(neg_log_p ~ position, data = .x, span = 0.5)
#     pred <- predict(object = fit, newdata = .x, se = TRUE)
#     .x$fitted <- pred$fit
#     .x$se <- pred$se.fit
#     .x
#   })
#
#
# dat3_seg_dist_smooth %>%
#   # subset(chrom == "chr_12") %>%
#   ggplot(aes(x = position / 1e6)) +
#   geom_point(aes(y = neg_log_p), size = 0.5, color = "grey50") +
#   # geom_ribbon(aes(ymin = fitted - (3 * se), ymax = fitted + (10 * se)), fill = alpha("white", 0), color = "slateblue") +
#   geom_line(aes(y = fitted + 12), lwd = 1, color = "slateblue") +
#   facet_wrap(~ chrom, scales = "free_x", nrow = 3) +
#   scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
#   scale_y_continuous(name = "-log10(P-value)", breaks = pretty, guide = guide_axis(cap = TRUE)) +
#   scale_color_manual(values = chrom_colors, guide = "none") +
#   theme_classic() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
#         panel.spacing.x = unit(0, "line"), strip.background = element_blank(), strip.placement = "outside")
#
#
# marker_keep <- subset(dat3_seg_dist_smooth, neg_log_p <= fitted + 12, marker, drop = TRUE)
# markers_drop <- setdiff(dat3_seg_dist_smooth$marker, marker_keep)
#
# dat1 <- filter_data(dat1, mrk.thresh = 1, ind.thresh = 1, chisq.pval.thresh = 1)
# dat1$screened.data$thresholds$miss.mrk <- max_mar_miss
# dat1$screened.data$thresholds$miss.ind <- max_ind_miss
# dat1$screened.data$ind.names <- colnames(dat1_geno_dose)
# dat1$screened.data$mrk.names <- intersect(row.names(dat1_geno_dose), marker_keep)
#
# ## Plot
# plot(dat1)
# plot(dat1, type = "density")
#
# # Estimate pairwise recombination
# dat1_twp <- pairwise_rf(x = dat1, mrk.scope = "all", ncpus = 12)
# plot(dat1_twp, type = "rf")
#
# # Filter based on recombination fraction
# dat2_twp <- rf_filter(dat1_twp,
#                       thresh.LOD.ph = 5,
#                       thresh.LOD.rf = 5,
#                       thresh.rf = 0.15,
#                       probs = c(0.01, 0.99))
#
# dat2_twp
#
# plot(dat2_twp)
#
#
# # Group markers
# dat2_grp <- group(x = dat2_twp, expected.groups = 12, comp.mat = TRUE, inter = FALSE)
# print(dat2_grp)
# plot(dat2_grp)
#
# # Create sequences
# dat2_seq <- make_sequence(x = dat2_grp, ch = as.list(1:12))
# print(x = dat2_seq, type = "genome")
#
# # Order markers according to physical genomic position and then linkage information
# dat2_ord <- order_sequence(x = dat2_seq, type = "genome")
# print(dat2_ord, type = "genome")
#
#
# # Phase the parents
# dat2_phase <- pairwise_phasing(x = dat2_ord,
#                                type = "genome",
#                                thresh.LOD.ph = 3,
#                                thresh.LOD.rf = 3,
#                                thresh.rf = 0.5,
#                                max.search.expansion.p1 = 10,
#                                max.search.expansion.p2 = 10)
#
# print(dat2_phase, type = "genome")
#
#
# ## Estimate maps
# ##
# ## Parent 1 + Parent 2
# dat2_map <- mapping(dat2_phase, type = "genome", parent = "p1", ncpus = 12,  )
# print(dat2_map, type = "genome")
#
# # Calculate haplotypes initially
# dat2_hap_init <- calc_haplotypes(x = dat2_map, type = "genome", ncpus = 12)






