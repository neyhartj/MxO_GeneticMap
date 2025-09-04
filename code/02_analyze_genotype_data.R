## MxO Genetic map
##
## Analyze the genotypic data
##
## This script will analyze genotype frequencies and segregation patterns
##

# Load packages
library(tidyverse)
library(qtl)
library(qtl2)
library(ggpubr)

pop_metadata <- read_csv("data/population_metadata.csv")

# Colors
col <- c("slateblue", "violetred", "green3")

# Load data into qtl ------------------------------------------------------

# Start with the larger family
qtl_files <- list.files(path = "data", pattern = "CNJ98-325-33", full.names = FALSE)

# Read into qtl
cross <- read.cross(format = "csvs", dir = "data", genfile = qtl_files[1], phefile = qtl_files[2], crosstype = "f2")
# Count missing marker proportion
# Remove markers with > 10% missing
marker_missing <- nmissing(cross = cross, what = "mar") / nind(cross)
hist(marker_missing)
markers_missing_remove <- names(marker_missing)[marker_missing > 0.10]

# Identify duplicate markers
dup_markers <- findDupMarkers(cross = cross)
dup_markers_remove <- unlist(dup_markers)

cross_1 <- drop.markers(cross = cross, markers = union(dup_markers_remove, markers_missing_remove))
c(cross1 = sum(nmar(cross_1)), cross = sum(nmar(cross)))

# Remove high missing individuals
ind_missing <- nmissing(cross = cross, what = "ind") / sum(nmar(cross))
hist(ind_missing)
ind_missing_remove <- which(ind_missing > 0.10)
cross_1 <- cross_1[,-ind_missing_remove]


# Get the map
cross_map_df <- map2table(pull.map(cross_1)) %>%
  rownames_to_column("marker") %>%
  mutate(chrom = paste0("chr", str_pad(chr, 2, "left", "0")))

# Convert to qtl2
cross2 <- convert2cross2(cross = cross_1)
# Extract the map
cross_map <- pull.map(cross = cross_1)
# Calculate genotype probabilities
cross2_genoprob <- calc_genoprob(cross = cross2, map = cross_map)


# Calculate genotype frequencies
cross_geno_tabl <- geno.table(cross = cross_1) %>%
  rownames_to_column("marker") %>%
  left_join(., cross_map_df)

cross_geno_tabl1 <- cross_geno_tabl %>%
  mutate(allele_count = 2* (AA + AB + BB),
         pMacro = ((2 * AA) + AB) / allele_count,
         pOxy = ((2 * BB) + AB)  / allele_count,
         pAA = AA / (allele_count/2),
         pAB = AB / (allele_count/2),
         pBB = BB / (allele_count/2),
         pMissing = (missing * 2) / allele_count)



# Plot
gg1 <- cross_geno_tabl1 %>%
  ggplot(aes(x = pos, y = pMacro)) +
  geom_point(color = "grey85", size = 0.5) +
  geom_smooth(method = "loess") +
  geom_hline(yintercept = c(0.42, 0.58), linetype = 2) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_continuous(name = "Position (Mb)", breaks = pretty) +
  scale_y_continuous(name = "V. macrocarpon (Stevens) allele frequency", breaks = pretty) +
  # coord_cartesian(ylim = c(0.2)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.spacing.x = unit(0.2, "lines"))
gg1

# Save
ggsave(filename = "mxo_Vmac_allele_frequency.jpg", plot = gg1, path = "output/", width = 6, height = 4, dpi = 500)

# Plot
gg2 <- cross_geno_tabl1 %>%
  ggplot(aes(x = pos, y = -log10(P.value))) +
  geom_point(color = "grey85", size = 0.5) +
  geom_smooth(method = "loess") +
  facet_grid(~ chrom, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_continuous(name = "Position (Mb)", breaks = pretty) +
  scale_y_continuous(name = "-log10(p-value) for segregation test", breaks = pretty) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.spacing.x = unit(0.2, "lines"))
gg2


# Plot genotype probabilities ----------------------------------------------------------

max_genoprob <- maxmarg(probs = cross2_genoprob, return_char = TRUE)

qtl_geno_df <- do.call("cbind", max_genoprob) %>%
  as.data.frame() %>%
  rownames_to_column("geno_id") %>%
  gather(marker, genotype, -geno_id) %>%
  left_join(., cross_map_df) %>%
  mutate(chrom = fct_inorder(chrom)) %>%
  arrange(chrom, pos)

g_geno <- qtl_geno_df %>%
  mutate(genotype = as.factor(genotype),
         marker = fct_inorder(marker)) %>%
  ggplot(aes(x = marker, y = geno_id)) +
  geom_tile(aes(fill = genotype)) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(name = NULL) +
  scale_y_discrete(name = "Clone") +
  scale_fill_manual(values = c("AA" = "slateblue", "AB" = "violetred", "BB" = "green3", "NA" = "black"), labels = c("Vmac", "Het", "Voxy", "NA")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing.x = unit(0.1, "lines"), strip.background = element_blank())
ggsave(filename = "mxo_genotypes_family_CNJ98-325-33.jpg", plot = g_geno, path = "output/", width = 10, height = 6, dpi = 500)


g_geno1 <- qtl_geno_df %>%
  subset(chrom %in% c("chr01", "chr02", "chr12")) %>%
  mutate(genotype = as.factor(genotype),
         marker = fct_inorder(marker)) %>%
  ggplot(aes(x = marker, y = geno_id)) +
  geom_tile(aes(fill = genotype)) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(name = NULL) +
  scale_y_discrete(name = "Clone") +
  scale_fill_manual(name = "Genotype", values = c("AA" = "slateblue", "AB" = "violetred", "BB" = "green3", "NA" = "black"), labels = c("AA/Vmac", "AB/Het", "BB/Voxy", "NA")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing.x = unit(0.1, "lines"), strip.background = element_blank())
ggsave(filename = "mxo_genotypes_family_CNJ98-325-33_subset.jpg", plot = g_geno1, path = "output/", width = 5, height = 5, dpi = 1000)





# Add phenotypic data to the cross object ---------------------------------

# traits to study
traits <- c("DaysToFlw", "AvgFirstFrtWt", "PollenViability")
# Trait labels
trait_labs <- setNames(c("Days to flowering", "Avg. fruit weight", "Pollen viability"), traits)

pheno_df <- read_csv(file = "data/mxo_all_family_pheno.csv") %>%
  as.data.frame() %>%
  column_to_rownames("id")

# Plot phenotypic data with parents
pheno_df1 <- pheno_df %>%
  rownames_to_column("id") %>%
  gather(trait, BLUE, -id) %>%
  subset(trait %in% traits)

# Parents
parents <- c("STEVENS", "NJ96-20")
f1 <- na.omit(unique(pop_metadata$S0_name))
ids <- subset(pop_metadata, !is.na(family_name), individual, drop = TRUE)

gg_pheno <- pheno_df1 %>%
  subset(id %in% ids) %>%
  ggplot(aes(x = BLUE)) +
  geom_histogram(fill = "grey85", col = "black", lwd = 0.25) +
  geom_vline(data = subset(pheno_df1, id %in% c(parents, f1)), aes(xintercept = BLUE, color = id)) +
  scale_x_continuous(name = "Phenotypic value", breaks = pretty) +
  scale_y_continuous(name = NULL, breaks = pretty) +
  scale_color_manual(values = setNames(c(col[1], col[2], col[2], col[3]), c(parents[1], f1, parents[2])),
                     labels = c("Hybrid", "Voxy", "Vmac"), name = "Parent/F1", guide = guide_legend(override.aes = list(linewidth = 3))) +
  facet_wrap(~ trait, scales = "free_x", switch = "x", labeller = labeller(trait = trait_labs)) +
  theme_minimal() +
  theme(strip.placement = "outside", axis.text.y = element_blank())

ggsave(filename = "mxo_phenotypic_dist.jpg", plot = gg_pheno, path = "output", width = 8, height = 3, dpi = 1000)



pheno_mat <- as.matrix(pheno_df[ind_ids(cross2),])

cross2$covar <- NULL
cross2$pheno <- pheno_mat

summary(cross2)


# Map! --------------------------------------------------------------------

# Calculate kinship
kin <- calc_kinship(probs = cross2_genoprob, type = "loco")

# Run genomewide scan using the mixed model
lmm_scan_out <- scan1(genoprobs = cross2_genoprob, pheno = cross2$pheno, kinship = kin, reml = TRUE)
# Run a permutation test
lmm_scan_perm_out <- scan1perm(genoprobs = cross2_genoprob, pheno = cross2$pheno, kinship = kin, reml = TRUE, n_perm = 1000, cores = 8)
(lmm_scan_thresholds <- summary(lmm_scan_perm_out, alpha = c(0.01, 0.05, 0.10, 0.20)))

# Save
save("lmm_scan_out", "lmm_scan_thresholds", file = "output/lmm_scan_results.RData")
load("output/lmm_scan_results.RData")

# Plot
ymx <- maxlod(lmm_scan_out)
i = 1
plot(lmm_scan_out, map = cross_map, lodcolumn = i, col = "slateblue", main = colnames(lmm_scan_out)[i], ylim = c(0, ymx+1))

lmm_scan_out_df <- lmm_scan_out %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  select(marker, all_of(traits)) %>%
  left_join(., cross_map_df) %>%
  gather(trait, score, -marker, -chrom, -chr, -pos)

g_lmm_scan <- lmm_scan_out_df %>%
  ggplot(aes(x = pos, y = score, color = chrom)) +
  geom_line() +
  geom_hline(data = subset(gather(rownames_to_column(as.data.frame(lmm_scan_thresholds), "signif"), trait, thresh, -signif), trait %in% traits & signif == "0.05"),
             aes(yintercept = thresh), lty = 2, lwd = 0.5) +
  facet_grid(trait ~ chr, switch = "both", labeller = labeller(trait = trait_labs), scales = "free_x", space = "free_x") +
  scale_color_manual(guide = "none", values = rep(c("slateblue", "grey50"), length.out = 12)) +
  scale_y_continuous(name = "LOD score", breaks = pretty, limits = c(0, ymx)) +
  scale_x_continuous(name = NULL, labels = NULL) +
  theme_minimal() +
  theme(strip.placement = "outside", panel.spacing.x = unit(0, "line"), panel.spacing.y = unit(1.5, "line"))
ggsave(filename = "lmm_scan_all_traits.jpg", plot = g_lmm_scan, path = "output", width = 8, height = 5, dpi = 1000)


# Plot just pollen viability
g_lmm_scan <- lmm_scan_out_df %>%
  subset(trait == "PollenViability" & chr %in% c(1, 2, 12)) %>%
  ggplot(aes(x = pos, y = score, color = chrom)) +
  geom_line() +
  facet_grid(trait ~ chr, switch = "both", labeller = labeller(trait = trait_labs), scales = "free_x", space = "free_x") +
  scale_color_manual(guide = "none", values = rep(c("slateblue", "grey50"), length.out = 12)) +
  scale_y_continuous(name = "LOD score", breaks = pretty, limits = c(0, ymx)) +
  scale_x_continuous(name = NULL, labels = NULL) +
  theme_minimal() +
  theme(strip.placement = "outside", panel.spacing.x = unit(0, "line"), panel.spacing.y = unit(1.5, "line"))
ggsave(filename = "lmm_scan_pollen_viability_subset.jpg", plot = g_lmm_scan, path = "output", width = 4, height = 3, dpi = 1000)



# Plot three traits
for (trt in traits) {
  i <- which(colnames(lmm_scan_out) == trt)
  jpeg(filename = file.path("output", paste0("lmm_scan1_", trt, "_out.jpg")), height = 5, width = 8, units = "in", res = 1000)
  plot(lmm_scan_out, map = cross_map, lodcolumn = i, col = "slateblue", main = colnames(lmm_scan_out)[i], ylim = c(0, ymx+1))
  abline(h = lmm_scan_thresholds["0.05",trt], lty = 2)
  dev.off()
}


# Plot pollen viability results
plot(lmm_scan_out, map = cross_map, lodcolumn = 6, col = "slateblue", main = colnames(lmm_scan_out)[i], ylim = c(0, ymx+1), chr = c(1,2,12))


# Get QTL
peaks <- find_peaks(scan1_output = lmm_scan_out, map = cross_map, threshold = lmm_scan_thresholds["0.1",], drop = 1.5)
peaks
# bayes_interval <- bayes_int(scan1_output = lmm_scan_out, map = cross_map)

# Fit QTL effects
i = 1
trait <- colnames(cross2$pheno)[i]

genoprobs_i <- pull_genoprobpos(genoprobs = cross2_genoprob, map = cross_map, chr = peaks$chr[peaks$lodcolumn == trait], peaks$pos[peaks$lodcolumn == trait])
qtl_fit <- fit1(genoprobs = genoprobs_i, pheno = cross2$pheno[,i], kinship = kin[["2"]], )

# Colors
col <- c("slateblue", "violetred", "green3")

# Estimate QTL effects
j = 2
g <- maxmarg(probs = cross2_genoprob, map = cross_map, chr = peaks$chr[j], pos = peaks$pos[j])
qtl_eff <- scan1coef(genoprobs = cross2_genoprob[,peaks$chr[j]], pheno = cross2$pheno[,peaks$lodcolumn[j]], kinship = kin[[peaks$chr[j]]],
                     se = TRUE, reml = TRUE, zerosum = TRUE)
plot(qtl_eff, map = cross_map, columns = 1:3, col = col)
last_coef <- unclass(qtl_eff)[nrow(qtl_eff),]
for(i in seq(along=last_coef)) axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis = col[i])
(qtl_eff_j <- qtl_eff[paste0("chr", peaks$chr[j], "_", peaks$pos[j]),])
plot_pxg(geno = g, pheno = cross2$pheno[,peaks$lodcolumn[j]] - qtl_eff_j["intercept"], sort = TRUE, SEmult = 2, swap_axes = TRUE, jitter = 0.1,
         xlab = peaks$lodcolumn[j])

j = 2
g <- maxmarg(probs = cross2_genoprob, map = cross_map, chr = peaks$chr[j], pos = peaks$pos[j])
qtl_eff <- scan1coef(genoprobs = cross2_genoprob[,peaks$chr[j]], pheno = cross2$pheno[,peaks$lodcolumn[j]], kinship = kin[[peaks$chr[j]]],
                     se = TRUE, reml = TRUE, zerosum = TRUE)
plot(qtl_eff, map = cross_map, columns = 1:3, col = col)
last_coef <- unclass(qtl_eff)[nrow(qtl_eff),]
for(i in seq(along=last_coef)) axis(side=4, at=last_coef[i], labels = names(last_coef)[i], tick=FALSE, col.axis = col[i])
# Plot with lod score
plot(qtl_eff, map = cross_map, columns = 1:3, col = col, scan1_output = lmm_scan_out, legend = "topright")
(qtl_eff_j <- qtl_eff[paste0("chr", peaks$chr[j], "_", peaks$pos[j]),])
plot_pxg(geno = g, pheno = cross2$pheno[,peaks$lodcolumn[j]] - qtl_eff_j["intercept"], sort = TRUE, SEmult = 2, swap_axes = TRUE, jitter = 0.1,
         xlab = peaks$lodcolumn[j])

j = 6
g <- maxmarg(probs = cross2_genoprob, map = cross_map, chr = peaks$chr[j], pos = peaks$pos[j])
qtl_eff <- scan1coef(genoprobs = cross2_genoprob[,peaks$chr[j]], pheno = cross2$pheno[,peaks$lodcolumn[j]], kinship = kin[[peaks$chr[j]]],
                     se = TRUE, reml = TRUE, zerosum = TRUE)
plot(qtl_eff, map = cross_map, columns = 1:3, col = col)
last_coef <- unclass(qtl_eff)[nrow(qtl_eff),]
for(i in seq(along=last_coef)) axis(side=4, at=last_coef[i], labels = names(last_coef)[i], tick=FALSE, col.axis = col[i])
# Plot with lod score
plot(qtl_eff, map = cross_map, columns = 1:3, col = col, scan1_output = lmm_scan_out, legend = "topright")
(qtl_eff_j <- qtl_eff[paste0("chr", peaks$chr[j], "_", peaks$pos[j]),])
plot_pxg(geno = g, pheno = cross2$pheno[,peaks$lodcolumn[j]] - qtl_eff_j["intercept"], sort = FALSE, SEmult = 2, swap_axes = TRUE, jitter = 0.1,
         xlab = peaks$lodcolumn[j])


