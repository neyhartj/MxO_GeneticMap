## MxO Genetic map
##
## Analyze the genotypic data
##
## This script will analyze genotype frequencies and segregation patterns
##

# Load packages
library(tidyverse)
library(qtl2)
library(ggpubr)

pop_metadata <- read_csv("data/population_metadata.csv")

# Colors
col <- c("slateblue", "violetred", "green3")



# Load data ------------------------------------------------------

# Read in the interpolated map
gmap_tab <- read_csv("results/mxo_rqtl_F1-CNJ98-325-33_BenLear_interpolated_map.csv")

# Load the qtl2 cross object
data_files_load <- list.files(data_dir, pattern = "rqtl2", full.names = TRUE)
for (file in data_files_load) {
  load(file)
}

# Add the map to the cross object
cross_qtl2$pmap <- sapply(cross_qtl2$gmap, function(x) x * 1e6)
gmap <- gmap_tab %>%
  as.data.frame() %>%
  select(marker, chr, pos = cM) %>%
  column_to_rownames("marker") %>%
  qtl::table2map()

cross_qtl2$gmap <- gmap


# Load phenotypes
pheno <- read_csv(file = list.files(path = "data_raw/", pattern = "pheno.csv", full.names = TRUE))

# Add the phenotypes to the cross object
pheno_add <- cross_qtl2$pheno %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  select(-fake) %>%
  left_join(., pheno) %>%
  column_to_rownames("id") %>%
  as.matrix()

# Recalculate flowering interval
pheno_add[,"FlwInterval"] <- pheno_add[,"DaysToSecFlw"] - pheno_add[,"DaysToFlw"]

# Add reflowering as a binary trait
pheno_add <- cbind(pheno_add, ReflowerBinary = as.numeric(!is.na(pheno_add[,"FlwInterval"])))
# Add NA if first flower is NA
pheno_add[is.na(pheno_add[,"DaysToFlw"]), "ReflowerBinary"] <- NA

cross_qtl2_map <- cross_qtl2
cross_qtl2_map$pheno <- pheno_add
# Remove covariate
cross_qtl2_map$covar <- NULL

cross_qtl2_map


##

# Run genome scans -------------------------------------------------------------

# Trait names
traits <- pheno_names(cross_qtl2_map)
trait_labs <- neyhart::str_add_space(traits)
# Binary traits
binary_traits <- "ReflowerBinary"
ols_traits <- setdiff(traits, binary_traits)

# Get the map
gmap <- insert_pseudomarkers(map = cross_qtl2_map$gmap)
pmap <- cross_qtl2_map$pmap

# Calculate genotype probabilities
genoprob <- calc_genoprob(cross = cross_qtl2_map, map = gmap, error_prob = 0.06)

# Count crossovers
genos <- maxmarg(probs = genoprob)
crossovers <- count_xo(geno = genos, cores = 8)
crossovers_sum <- rowSums(x = crossovers)

# Choose genotypes without excessive crossover counts
geno_select <- names(subset(crossovers_sum, crossovers_sum <= 50))
# geno_select <- names(crossovers_sum)

cross_qtl2_map2 <- subset(cross_qtl2_map, ind = geno_select)
genoprob <- calc_genoprob(cross = cross_qtl2_map2, map = gmap, error_prob = 0.06)
genos <- maxmarg(probs = genoprob)
crossovers <- count_xo(geno = genos, cores = 8)

cross_qtl2_map2$pheno <- cbind(cross_qtl2_map2$pheno, Crossovers = rowSums(crossovers))


# Calculate kinship
kin <- calc_kinship(probs = genoprob, type = "loco")

# Phenotypes for OLS
phenos_ols <- cross_qtl2_map2$pheno
# Run genomewide scan using the mixed model
lmm_scan_out_ols <- scan1(genoprobs = genoprob, pheno = phenos_ols, kinship = kin, cores = 8)

# # Compare with model without kinship
# scan_out_ols <- scan1(genoprobs = genoprob, pheno = phenos_ols, cores = 8)
# plot(lmm_scan_out_ols, map = gmap, lodcolumn = 1, col = "firebrick")
# plot(scan_out_ols, map = gmap, lodcolumn = 1, add = TRUE, col = "slateblue")


# # Phenotypes for binary
# phenos_bin <- cross_qtl2_map2$pheno[, binary_traits, drop = FALSE]
# # Run genomewide scan using the mixed model
# lmm_scan_out_bin <- scan1(genoprobs = genoprob, pheno = phenos_bin, cores = 8, model = "binary")
#
# # Compare results for this trait using the normal model versus binary
# ## Pretty much the same answer - just use the LMM model
# plot(lmm_scan_out_ols, map = gmap, lodcolumn = binary_traits, col = "firebrick")
# plot(lmm_scan_out_bin, map = gmap, lodcolumn = binary_traits, add = TRUE, col = "slateblue")
#
# # Merge
# scan1_ans <- cbind(lmm_scan_out_ols, lmm_scan_out_bin)

scan1_ans <- lmm_scan_out_ols

# Run a permutation test
scan1_ans_perm <- scan1perm(genoprobs = genoprob, pheno = phenos_ols, kinship = kin, reml = TRUE, n_perm = 10000, cores = 12)
(scan1_ans_perm_thresholds <- summary(scan1_ans_perm, alpha = c(0.01, 0.05, 0.10, 0.20)))


# Save
save_file <- file.path(results_dir, "qtl_mapping/scan_results.RData")
save("scan1_ans", "scan1_ans_perm", "scan1_ans_perm_thresholds", file = save_file)
load(save_file)


# Plot
(trait_names <- colnames(scan1_ans))
i = 6
trt <- trait_names[i]
plot(scan1_ans, map = pmap, lodcolumn = i, main = trt)
add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.05",])

# Zoom in
plot(scan1_ans, map = gmap, lodcolumn = i, main = trt, chr = 2)




# Estimate additive and dominance effects
# chromosome 2 reflowering
trt <- trait_names[2]
chr <- "9"
eff <- scan1coef(genoprob[,chr], phenos_ols[,trt], kinship = kin[[chr]], reml = FALSE)
plot(eff, map = gmap, columns = 1:3)
last_coef <- unclass(eff)[nrow(eff),] # pull out last coefficients
for(i in seq(along=last_coef)) axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

eff <- scan1coef(genoprob[,chr], phenos_ols[,trt], kinship = kin[[chr]], contrasts=cbind(mu=c(1,1,1), a=c(2, 1, 0), d=c(0, 1, 0)))
plot(eff, map = gmap, columns = 2:3)
last_coef <- unclass(eff)[nrow(eff),] # pull out last coefficients
for(i in seq(along=last_coef)) axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


# chromosome 1 pollen viability
trt <- trait_names[7]
chr <- "4"
eff <- scan1coef(genoprob[,chr], phenos_ols[,trt], kinship = kin[[chr]], reml = FALSE)
plot(eff, map = pmap, columns = 1:3, legend = "topright", scan1_output = scan1_ans[,trt, drop = FALSE], main = trt)

eff <- scan1coef(genoprob[,chr], phenos_ols[,trt], kinship = kin[[chr]], contrasts=cbind(mu=c(1,1,1), a=c(2, 1, 0), d=c(0, 1, 0)))
plot(eff, map = pmap, columns = 2:3, legend = "topright", scan1_output = scan1_ans[,trt, drop = FALSE], main = trt)


## Create plot of the whole-genome scan results
ymx <- maxlod(scan1_ans)

# DF of the map
cross_map_df <- gmap_tab %>%
  rename(chrom = chr)

# df of threshold information
scan1_ans_perm_thresholds_df <- scan1_ans_perm_thresholds %>%
  as.data.frame() %>%
  rownames_to_column("signif") %>%
  gather(trait, thresh, -signif) %>%
  subset(trait %in% traits & signif == "0.05")

scan1_ans_df <- scan1_ans %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  select(marker, all_of(traits)) %>%
  left_join(., cross_map_df) %>%
  gather(trait, score, -marker, -chrom, -pos, -cM) %>%
  mutate(chrom = fct_inorder(as.character(chrom)))

# Iterate over traits and create a plot
for (trt in traits) {
  dat <- subset(scan1_ans_df, trait == trt)
  thresh <- subset(scan1_ans_perm_thresholds_df, trait == trt, thresh, drop = TRUE)

  g_scan_trt <- dat %>%
    ggplot(aes(x = pos, y = score, color = chrom)) +
    geom_line() +
    geom_hline(yintercept = thresh, lty = 2, lwd = 0.5) +
    facet_grid(. ~ chrom, switch = "both", labeller = labeller(trait = trait_labs), scales = "free_x", space = "free_x") +
    scale_color_manual(guide = "none", values = rep(c("slateblue", "grey50"), length.out = 12)) +
    scale_y_continuous(name = "LOD score", breaks = pretty, limits = c(0, ymx), guide = guide_axis(cap = TRUE)) +
    scale_x_continuous(name = "Physical position", labels = NULL) +
    labs(subtitle = trt) +
    theme_classic() +
    theme(strip.placement = "outside", panel.spacing.x = unit(0, "line"), panel.spacing.y = unit(1.5, "line"),
          strip.background.x = element_blank(), axis.ticks.x = element_blank())

  # Save
  filename <- paste0("scan1_plot_", trt, ".jpg")
  ggsave(filename = filename, plot = g_scan_trt, path = fig_dir, width = 8, height = 4, dpi = 500)

}




# Summarize QTL peaks -----------------------------------------------------

# Helpful plotting to look at lod scores
plot(scan1_ans, map = pmap, lodcolumn = 1, chr = 2); add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.05",])

peaks05 <- find_peaks(scan1_output = scan1_ans, map = pmap, threshold = scan1_ans_perm_thresholds["0.05",], drop = 1.5, peakdrop = 2)
peaks10 <- find_peaks(scan1_output = scan1_ans, map = pmap, threshold = scan1_ans_perm_thresholds["0.1",], drop = 1.5, peakdrop = 2)
peaks05
peaks10

# Estimate kinship across genome (no LOCO)
kin_full <- calc_kinship(probs = genoprob, type = "overall")

## Start trait loop here

(trt <- traits[6])
peaks_trt <- subset(peaks05, lodcolumn == trt)
# Build the additive design matrix for each marker
X <- NULL
for (j in seq_len(nrow(peaks_trt))) {
  # Pull the genotype probabilieis
  prob_j <- pull_genoprobpos(genoprobs = genoprob, map = pmap, chr = peaks_trt[j,"chr"], pos = peaks_trt[j,"pos"])
  mar_j <- find_marker(map = pmap, chr = peaks_trt[j,"chr"], pos = peaks_trt[j,"pos"])
  add_dos_j <- prob_j %*% matrix(c(1, 0, -1)) # Positive = Stevens
  dom_dos_j <- prob_j %*% matrix(c(0, 1, 0))
  # Create a matrix
  X_j <- cbind(add_dos_j, dom_dos_j)
  colnames(X_j) <- c(mar_j, paste0(mar_j, ".dom"))
  X <- cbind(X, X_j)
}

# Fit the QTL model
dat_use <- as.data.frame(X)
dat_use$id <- row.names(phenos_ols)
dat_use[[trt]] <- as.numeric(phenos_ols[,trt])
# Fit the null model
fixed0 <- reformulate(termlabels = "1", response = trt)
fit0 <- mmes(fixed = fixed0, random = ~ vsm(ism(id), Gu = kin_full), data = dat_use, verbose = FALSE)

mar_names <-  colnames(X)
fixed_full <- reformulate(termlabels = mar_names, response = trt)
fit_full <- mmes(fixed = fixed_full, random = ~ vsm(ism(id), Gu = kin_full), data = dat_use, verbose = FALSE)

remove_markers <- TRUE
while(remove_markers) {
  # Backwards elimination of interactions
  elim <- data.frame(marker = mar_names, p_value = as.numeric(0))
  for (k in seq_len(nrow(elim))) {
    fixed_red <- formula(drop.terms(termobj = terms(fixed_full), dropx = k, keep.response = TRUE))
    fit_red <- mmes(fixed = fixed_red, random = ~ vsm(ism(id), Gu = kin_full), data = dat_use, verbose = FALSE)
    stdout <- capture.output(lrt <- anova(object = fit_full, object2 = fit_red))
    elim$p_value[k] <- pchisq(q = as.numeric(lrt$Chisq[-1]), df = as.numeric(lrt$ChiDf[-1]), lower.tail = FALSE)
  }
  # Remove the marker with the highest p-value above 0.05
  elim <- subset(elim, p_value > 0.05)
  remove_markers <- nrow(elim) > 0
  mar_elim <- elim$marker[which.max(elim$p_value)]
  mar_names <- setdiff(mar_names, mar_elim)
  fixed_full <- reformulate(termlabels = mar_names, response = trt)
  fit_full <- mmes(fixed = fixed_full, random = ~ vsm(ism(id), Gu = kin_full), data = dat_use, verbose = FALSE)
}

# Add dominance terms
fixed_full <- reformulate(termlabels = c(mar_names, paste0(mar_names, ".dom")), response = trt)
fit_full <- mmes(fixed = fixed_full, random = ~ vsm(ism(id), Gu = kin_full), data = dat_use, verbose = FALSE)




# Map crossover count -----------------------------------------------------






covar <- maxmarg(probs = genoprob, map = pmap, chr = peaks05_trt[,"chr"], pos = peaks05_trt[,"pos"])
est_herit(pheno = phenos_ols[,trt], kinship = kin_full, addcovar = covar)

eff <- scan1coef(genoprob[,chr], phenos_ols[,trt], kinship = kin[[chr]], reml = FALSE)

# fit a qtl model
fit_model <- fit1(genoprobs = genoprob, pheno = phenos_ols[,trt], kinship = kin[[peaks05_trt[1,"chr"]]])
summary(fit_model)





g_scan <- scan1_ans_df %>%
  ggplot(aes(x = pos, y = score, color = chrom)) +
  geom_line() +
  geom_hline(data = scan1_ans_perm_thresholds_df,  aes(yintercept = thresh), lty = 2, lwd = 0.5) +
  facet_grid(trait ~ chrom, switch = "both", labeller = labeller(trait = trait_labs), scales = "free_x", space = "free_x") +
  scale_color_manual(guide = "none", values = rep(c("slateblue", "grey50"), length.out = 12)) +
  scale_y_continuous(name = "LOD score", breaks = pretty, limits = c(0, ymx)) +
  scale_x_continuous(name = NULL, labels = NULL) +
  theme_minimal() +
  theme(strip.placement = "outside", panel.spacing.x = unit(0, "line"), panel.spacing.y = unit(1.5, "line"))
g_scan
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

















# Other -------------------------------------------------------------------






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

