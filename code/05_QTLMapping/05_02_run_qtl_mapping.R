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
library(rrBLUP)

pop_metadata <- read_csv("data/population_metadata.csv")

# Colors
col <- c("slateblue", "violetred", "green3")



# Load data ------------------------------------------------------

# Load phenotypes
pheno_blues <- read_csv(file = file.path(data_dir, "mxo_all_family_pheno_blues.csv"))
pheno_cat <- read_csv(file = file.path(data_dir, "mxo_all_family_pheno_categorical.csv"))

# Load the cross objects
data_files_load <- list.files(data_dir, pattern = "rqtl2", full.names = TRUE)
for (file in data_files_load) {
  load(file)
  fam <- str_extract(basename(file), "CNJ[0-9]{2}-[0-9]{3}-[0-9]{2}")

  # Add the phenotypes to the cross object
  pheno_add <- cross_qtl$pheno %>%
    as.data.frame() %>%
    select(-fake) %>%
    left_join(., pheno_blues, by = c("id" = "geno_id")) %>%
    column_to_rownames("id")
  pheno_add2 <- pheno_add %>%
    as.matrix()

  cross_qtl$pheno <- pheno_add
  cross_qtl2$pheno <- pheno_add2
  # Remove covariate
  cross_qtl2$covar <- NULL

  # Rename
  cross <- ifelse(fam == "CNJ98-309-19", "BLxOxy", "STxOxy")
  assign(x = paste0("cross_qtl2_", cross), value = cross_qtl2)
  assign(x = paste0("cross_qtl_", cross), value = cross_qtl)


}

rm(cross_qtl)
rm(cross_qtl2)




# Use R/qtl ---------------------------------------------------------------

# Find duplicate individuals
cg <- comparegeno(cross = cross_qtl_STxOxy)
cg_dup <- which(cg > 0.99, arr=TRUE)
cg_dup <- cg_dup[cg_dup[,1] < cg_dup[,2],]
cg_dup

cross_qtl_STxOxy <- subset(x = cross_qtl_STxOxy, ind = -cg_dup[,2])

# Calcuate genotype probabilities
cross_qtl_STxOxy <- calc.genoprob(cross = cross_qtl_STxOxy, error.prob = 0.03, map.function = "kosambi")
# Add crossovers
cross_qtl_STxOxy$pheno <- cbind(cross_qtl_STxOxy$pheno, Crossovers = countXO(cross_qtl_STxOxy))

probs <- calc_genoprob(cross = cross_qtl2_STxOxy)
kin <- calc_kinship(probs = probs, type = "loco")


pheno <- pull.pheno(cross_qtl_STxOxy)
pheno_add <- select(pheno_cat, id = geno_id, year, ReFlwBin) %>%
  group_by(id, year) %>%
  slice(1) %>%
  ungroup() %>%
  spread(year, ReFlwBin) %>%
  rename_all(~c("id", "ReFlwBin2021", "ReFlwBin2022"))
pheno1 <- pheno %>%
  rownames_to_column("id") %>%
  left_join(., pheno_add) %>%
  column_to_rownames("id")
cross_qtl_STxOxy$pheno <- pheno1

pheno <- pull.pheno(cross_qtl_STxOxy)
(traits <- colnames(pheno))
i <- 9
(trt <- traits[i])

# Run scan1
scan_ans <- scanone(cross = cross_qtl_STxOxy, pheno.col = i, method = "hk")
plot(scan_ans, main = trt)

binary <- all(na.omit(unique(pull.pheno(cross = cross_qtl_STxOxy, pheno.col = trt))) %in% c(0, 1))
if (binary) {
  cim_scan_ans <- scanone(cross = cross_qtl_STxOxy, pheno.col = i, method = "hk", model = "binary")
} else {
  cim_scan_ans <- cim(cross = cross_qtl_STxOxy, n.marcovar = 5, pheno.col = trt, method = "hk", window = 10)
}
plot(cim_scan_ans,  main = trt)
# plot(cim_scan_ans,  main = trt, chr = 9)


# Run lmm scan1
lmm_scan_ans <- scan1(genoprobs = probs, kinship = kin, pheno = pheno[,trt, drop = FALSE])
plot(lmm_scan_ans, map = cross_qtl2_STxOxy$pmap, col = col[1], ylim = c(0, 14), main = trt)


## Plot all
scan_ans1 <- qtl2convert::scan_qtl_to_qtl2(scanone_output = scan_ans)
cim_scan_ans1 <- qtl2convert::scan_qtl_to_qtl2(scanone_output = cim_scan_ans)


plot(lmm_scan_ans, map = cross_qtl2_STxOxy$pmap, col = col[1], ylim = c(0, 14), main = trt)
plot(scan_ans1$scan1, map = cross_qtl2_STxOxy$pmap, add = TRUE, col = col[2])
plot(cim_scan_ans1$scan1, map = cross_qtl2_STxOxy$pmap, add = TRUE, col = col[3])
legend("topleft", col = col, legend = c("LMM", "HK", "CIM"), lwd = 2)

# Get peaks from CIM
peaks_gmap <- find_peaks(scan1_output = cim_scan_ans1$scan1, map = cross_qtl2_STxOxy$gmap, threshold = 4, peakdrop = 1.5, drop = 1.5)
peaks_pmap <- find_peaks(scan1_output = cim_scan_ans1$scan1, map = cross_qtl2_STxOxy$pmap, threshold = 4, peakdrop = 1.5, drop = 1.5)
peaks_pmap

# Fit a model
qtl <- makeqtl(cross = cross_qtl_STxOxy, chr = peaks_gmap$chr, pos = peaks_gmap$pos, what = "prob")
formula <- reformulate(qtl$altname, response = "y")

# stepwise selection
step_qtl <- stepwiseqtl(cross = cross_qtl_STxOxy, pheno.col = trt, max.qtl = 4, method = "hk", qtl = qtl,
                        refine.locations = TRUE)
step_qtl
formula <- formula(step_qtl)

# Fit the QTL model
qtl_model <- fitqtl(cross = cross_qtl_STxOxy, pheno.col = trt, qtl = step_qtl, formula = y ~ Q1 + Q2 + Q3 + Q2:Q3,
                    method = "hk", get.ests = T)
summary(qtl_model)


find_markerpos(cross_qtl2_STxOxy, find.marker(cross = cross_qtl_STxOxy, chr = step_qtl$chr, pos = step_qtl$pos))

# plot
mar_name <- find.marker(cross = cross_qtl_STxOxy, chr = step_qtl$chr, pos = step_qtl$pos)
cross_qtl_STxOxy_sim <- sim.geno(cross = cross_qtl_STxOxy, error.prob = 0.03)
effectplot(cross = cross_qtl_STxOxy_sim, pheno.col = trt, mname1 = mar_name[1])
eff_scan <- effectscan(cross_qtl_STxOxy_sim, pheno.col = trt, chr = 3, get.se = TRUE)


plotPXG(x = cross_qtl_STxOxy, marker = mar_name[1], pheno.col = trt)



 # Run genome scans -------------------------------------------------------------

# Trait names
traits <- pheno_names(cross_qtl2_STxOxy)
trait_labs <- neyhart::str_add_space(traits)


## Ben Lear X Oxy ----------------------------------------------------------

# Get the map
gmap <- insert_pseudomarkers(map = cross_qtl2_BLxOxy$gmap)
pmap <- cross_qtl2_BLxOxy$pmap

# Calculate genotype probabilities
genoprob <- calc_genoprob(cross = cross_qtl2_BLxOxy, map = gmap, error_prob = 0.07)

# Count crossovers
genos <- maxmarg(probs = genoprob)
crossovers <- count_xo(geno = genos, cores = 8)
crossovers_sum <- rowSums(x = crossovers)
plot(crossovers_sum)

# Choose genotypes without excessive crossover counts
geno_select <- names(subset(crossovers_sum, crossovers_sum <= 50))
# geno_select <- names(crossovers_sum)

cross_qtl2_BLxOxy <- subset(cross_qtl2_BLxOxy, ind = geno_select)
genoprob <- calc_genoprob(cross = cross_qtl2_BLxOxy, map = gmap, error_prob = 0.07)
genos <- maxmarg(probs = genoprob)
crossovers <- count_xo(geno = genos, cores = 8)

cross_qtl2_BLxOxy$pheno <- cbind(cross_qtl2_BLxOxy$pheno, Crossovers = rowSums(crossovers))



### Run genomewide scan -----------------------------------------------------



# Calculate kinship
kin <- calc_kinship(probs = genoprob, type = "loco")

# Phenotypes for LMM
phenos_lmm <- cross_qtl2_BLxOxy$pheno
phenos_lmm <- phenos_lmm[,colSums(!is.na(phenos_lmm)) >= 15]
# Run genomewide scan using the mixed model
lmm_scan_out <- scan1(genoprobs = genoprob, pheno = phenos_lmm, kinship = kin, cores = 8)

# Prelim plotting
ymx <- max(pretty(range(lmm_scan_out)))
i <- 1
plot(lmm_scan_out, map = pmap, lodcolumn = i, main = colnames(phenos_lmm)[i], ylim = c(0, ymx))

# Compare with model without kinship
scan_out_ols <- scan1(genoprobs = genoprob, pheno = phenos_lmm, cores = 8)
i = 5
plot(lmm_scan_out, map = gmap, lodcolumn = i, col = "firebrick", ylim = c(0, ymx))
plot(scan_out_ols, map = gmap, lodcolumn = i, add = TRUE, col = "slateblue")
legend("topleft", legend = c("LMM", "HK"), col = c("firebrick", "slateblue"), lty = 1, lwd = 2)

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

# Run a permutation test
scan1_ans_perm <- scan1perm(genoprobs = genoprob, pheno = phenos_lmm, kinship = kin, reml = TRUE, n_perm = 1000, cores = 12)
(scan1_ans_perm_thresholds <- summary(scan1_ans_perm, alpha = c(0.01, 0.05, 0.10, 0.20)))


## Run chisq tests for categorical traits
cat_traits <- names(pheno_cat)[-1:-2]
cat_scan_out <- list()
for (trt in cat_traits) {
  dat <- pheno_cat[c("geno_id", "year", trt)]
  dat <- dat[!is.na(dat[[trt]]), ]
  dat <- subset(dat, geno_id %in% row.names(genoprob[[1]]))
  dat <- as.data.frame(dat)
  if (length(unique(dat$year)) == 1) {
    dat <- dat[c("geno_id", trt)]
  }

  # Index of rows in genoprobs to merge with dat
  idx <- match(x = dat$geno_id, table = row.names(genoprob[[1]]))
  dat <- dat[-1]

  form <- reformulate(termlabels = ".", response = trt)

  # fit the base model
  fit0 <- nnet::multinom(formula = form, data = dat, maxit = 500)

  scores_out <- NULL

  # Iterate over chromosome
  for (i in seq_along(genoprob)) {
    genoprob_i <- genoprob[[i]]
    scores_out_i <- matrix(0, nrow = dim(genoprob_i)[[3]], dimnames = list(dimnames(genoprob_i)[[3]], trt))
    for (j in seq_len(dim(genoprob_i)[[3]])) {
      # print(j)
      probs <- genoprob_i[idx,,j]
      dat_j <- cbind(dat, probs)
      out <- capture.output(fit <- nnet::multinom(formula = form, data = dat_j))
      # LRT
      lrt <- anova(fit,  fit0)
      # record the score
      scores_out_i[j,] <- -log10(lrt$`Pr(Chi)`[2])
    }
    scores_out <- rbind(scores_out, scores_out_i)
  }

  # Add to the list
  cat_scan_out_list[[trt]] <- list(score = scores_out, sample.size = length(unique(idx)))

}

# Bind columns
cat_scan_ans <- do.call(cbind, lapply(cat_scan_out_list, "[[", "score"))
attr(cat_scan_ans, "hsq") <- matrix(NA, nrow = length(genoprob), ncol = length(cat_scan_out_list), dimnames = list(seq_along(genoprob), names(cat_scan_out_list)))
attr(cat_scan_ans, "sample_size") <- do.call(c, lapply(cat_scan_out_list, "[[", "sample.size"))
class(cat_scan_ans) <- c("scan1", "matrix")

# Add to the scan output

scan1_ans <- lmm_scan_out_ols
scan1_ans <- cbind(scan1_ans, cat_scan_ans)
colnames(scan1_ans)


# Save
save_file <- file.path(results_dir, "qtl_mapping/scan_results.RData")
save("scan1_ans", "scan1_ans_perm", "scan1_ans_perm_thresholds", file = save_file)
load(save_file)


# Plot
(trait_names <- colnames(scan1_ans))
i = 10
trt <- trait_names[i]
plot(scan1_ans, map = pmap, lodcolumn = i, main = trt, ylim = c(0, ymx))
add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.05",], col = col[1])
add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.1",], col = col[2])


# Zoom in
plot(scan1_ans, map = gmap, lodcolumn = i, main = trt, chr = 2)




# Estimate additive and dominance effects
# chromosome 2 reflowering
i = 7
trt <- trait_names[i]
chr <- "5"
eff <- scan1coef(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]], reml = FALSE)
plot(eff, map = gmap, columns = 1:3, scan1_output = scan1_ans[,trt, drop = FALSE], main = trt)

eff <- scan1coef(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]], contrasts=cbind(mu=c(1,1,1), a=c(2, 1, 0), d=c(0, 1, 0)))
plot(eff, map = gmap, columns = 2:3)
last_coef <- unclass(eff)[nrow(eff),] # pull out last coefficients
for(j in seq(along=last_coef)) axis(side=4, at=last_coef[j], names(last_coef)[j], tick=FALSE, col.axis=col[j])


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




### Summarize QTL peaks -----------------------------------------------------

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








## Stevens X Oxy ----------------------------------------------------------

# Change alleles
cross_qtl2_STxOxy$alleles <- c("M", "O")

# Get the map
gmap <- insert_pseudomarkers(map = cross_qtl2_STxOxy$gmap)
pmap <- cross_qtl2_STxOxy$pmap

# Calculate genotype probabilities
genoprob <- calc_genoprob(cross = cross_qtl2_STxOxy, map = gmap, error_prob = 0.06)

# Count crossovers
genos <- maxmarg(probs = genoprob)
crossovers <- count_xo(geno = genos, cores = 8)
crossovers_sum <- rowSums(x = crossovers)
plot(crossovers_sum)

# Choose genotypes without excessive crossover counts
geno_select <- names(subset(crossovers_sum, crossovers_sum <= 50))
# geno_select <- names(crossovers_sum)

cross_qtl2_STxOxy <- subset(cross_qtl2_STxOxy, ind = geno_select)
genoprob <- calc_genoprob(cross = cross_qtl2_STxOxy, map = gmap, error_prob = 0.06)
genos <- maxmarg(probs = genoprob)
crossovers <- count_xo(geno = genos, cores = 8)

cross_qtl2_STxOxy$pheno <- cbind(cross_qtl2_STxOxy$pheno, Crossovers = rowSums(crossovers))




### Run genomewide scan -----------------------------------------------------

# Calculate kinship
kin <- calc_kinship(probs = genoprob, type = "loco")

# Phenotypes for LMM
phenos_lmm <- cross_qtl2_STxOxy$pheno
phenos_lmm <- phenos_lmm[,colSums(!is.na(phenos_lmm)) >= 15]

fit <- lm(ReFlwBin ~ geno_id + year, data = pheno_cat)
emm <- emmeans::emmeans(object = fit, specs = "geno_id")

# Get binary reflowering traits for each year
phenos_reflw <- pheno_cat %>%
  select(geno_id, year, ReFlwBin) %>%
  left_join(data.frame(geno_id = ind_ids(cross_qtl2_STxOxy)), .) %>%
  mutate(year = paste0("ReFlwBin_", year)) %>%
  spread(year, ReFlwBin) %>%
  left_join(., as.data.frame(emm)[,c("geno_id", "emmean")]) %>%
  as.data.frame() %>%
  column_to_rownames("geno_id") %>%
  rename(ReFlwBinBLUE = emmean) %>%
  as.matrix()

# Add to the cross object
phenos_lmm <- cbind(phenos_lmm, phenos_reflw)


# Run genomewide scan using the mixed model
lmm_scan_out <- scan1(genoprobs = genoprob, pheno = phenos_lmm, kinship = kin, cores = 8)

# Run a permutation test
scan1_ans_perm <- scan1perm(genoprobs = genoprob, pheno = phenos_lmm, kinship = kin, n_perm = 1000, cores = 12)
(scan1_ans_perm_thresholds <- summary(scan1_ans_perm, alpha = c(0.01, 0.05, 0.10, 0.20)))

scan1_ans <- lmm_scan_out


# Prelim plotting
ymx <- max(pretty(range(scan1_ans)))
(trait_names <- colnames(scan1_ans))
# Plot
i = 19
trt <- trait_names[i]
plot(scan1_ans, map = pmap, lodcolumn = i, main = trt, ylim = c(0, ymx))
add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.05",], col = col[1])
add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.1",], col = col[2])


# Plot LOD score plots for each trait
for (i in seq_len(ncol(scan1_ans))) {
  trt <- colnames(scan1_ans)[i]
  filename <- file.path(fig_dir, paste0("scan1_plot_", trt, ".jpg"))

  jpeg(filename = filename, width = 8, height = 5, units = "in", res = 300)
  plot(scan1_ans, map = pmap, lodcolumn = i, main = trt, ylim = c(0, ymx))
  add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.05", i], col = col[1])
  add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.1", i], col = col[2])
  legend(x = "topright", legend = c("0.05", "0.10"), title = "P value threshold", col = col[1:2], lwd = 3)
  dev.off()

}


# Save
save_file <- file.path(results_dir, "qtl_mapping/scan_results.RData")
save("scan1_ans", "scan1_ans_perm", "scan1_ans_perm_thresholds", file = save_file)
load(save_file)




### Fit QTL models ----------------------------------------------------------

# Find peaks
peaks <- find_peaks(scan1_output = scan1_ans, map = pmap, threshold = scan1_ans_perm_thresholds["0.1", ], peakdrop = 1.5,
                    prob = 0.90)

# Prep the cross object for modeling
ped <- data.frame(id = ind_ids(cross_qtl2_STxOxy), parent1 = "STEVENS", parent2 = "NJ96-20")
cross_qtl2_STxOxy_prep <- prep_cross2(cross = cross_qtl2_STxOxy, ped = ped, ploidy = 2, error.prob = 0.06,
                                      cores = 12)


# Unique traits to model
trt_model <- unique(peaks$lodcolumn)
# List to store output
qtl_model_list <- list()

# Iterate over traits
for (trt in trt_model) {

  # Get the markers to model
  peaks_trt <- peaks[peaks$lodcolumn == trt, ]
  chr <- peaks_trt$chr
  pos <- peaks_trt$pos

  markers_model <- find_marker(map = cross_qtl2_STxOxy$pmap, chr = chr, pos = pos)

  qtl <- data.frame(marker = markers_model, dominance = 2)
  epistasis <- as.data.frame(t(combn(x = markers_model, 2)))
  names(epistasis) <- c("marker1", "marker2")

  # Fit a QTL model and perform backwards elimination
  qtl_model_fit <- fit_qtl_lmm(data = cross_qtl2_STxOxy_prep, pheno = cross_qtl2_STxOxy$pheno[,trt, drop = FALSE],
                               trait = trt, ploidy = 2, qtl = qtl, epistasis = epistasis, polygenic = TRUE, max.cor.sing = 0.93)

  qtl_model_list[[trt]] <- qtl_model_fit

}




# Find peaks with a liberal threshold


# For each trait-chr combination with peaks, plot the effect plot and LOD score
peaks_unique <- unique(peaks[, c("lodindex", "lodcolumn", "chr")])
for (i in seq_len(nrow(peaks_unique))) {

  trt <- peaks_unique$lodcolumn[i]
  chr <- peaks_unique$chr[i]
  eff <- scan1blup(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]], cores = 8)
  # eff <- scan1coef(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]])

  # Where to place the legend
  idx <- which(eff == max(eff[,1:3]), arr.ind = TRUE)
  if (idx[,1] >= nrow(eff) / 2) {
    legend <- "topleft"
  } else {
    legend <- "topright"
  }


  filename <- file.path(fig_dir, paste0("scan1_blup_effect_plot_", trt, "_chr", str_pad(chr, 2, "left", "0"), ".jpg"))
  jpeg(filename = filename, width = 8, height = 6, units = "in", res = 300)
  plot(eff, map = pmap, columns = 1:3, scan1_output = scan1_ans[,trt, drop = FALSE], main = trt, legend = legend)
  # add_threshold(map = pmap, thresholdA = scan1_ans_perm_thresholds["0.05", trt])
  dev.off()

}


# # Extra peaks of interest
# i <- 18
# chr <- "6"
# trt <- colnames(scan1_ans)[i]
# eff <- scan1blup(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]], cores = 8)
#
# # Where to place the legend
# idx <- which(eff == max(eff[,1:3]), arr.ind = TRUE)
# if (idx[,1] >= nrow(eff) / 2) {
#   legend <- "topleft"
# } else {
#   legend <- "topright"
# }
#
# plot(eff, map = pmap, columns = 1:3, scan1_output = scan1_ans[,trt, drop = FALSE], main = trt, legend = legend)
#
# j <- 18
# chr <- as.character(peaks[18,"chr"])
# pos <- peaks[18, "pos"]
# trt <- peaks[18, "lodcolumn"]
# geno <- maxmarg(probs = genoprob, map = pmap, chr = chr, pos = pos, return_char = TRUE)
# plot_pxg(geno = geno, pheno = phenos_lmm[,trt], sort = T, SEmult = TRUE, omit_points = T)

peaks_mqm <- find_peaks(scan1_output = scan1_ans, map = pmap, threshold = scan1_ans_perm_thresholds["0.05", , drop = FALSE],
                        peakdrop = 2, drop = 2)

# Calculate overall kinship
kin_overall <- calc_kinship(probs = genoprob, type = "overall")

# Run QTL scans while holding major QTL constant
mqm_scan_ans <- list()
peak_traits <- unique(peaks_mqm$lodcolumn)

for (trt in peak_traits) {
  trt_peaks <- peaks_mqm[peaks_mqm$lodcolumn == trt, , drop = FALSE]
  covar_trt <- NULL

  while (nrow(trt_peaks) > 0) {
    which_lod_max <- order(trt_peaks$lod, decreasing = TRUE)
    peak_mar <- find_marker(map = pmap, chr = trt_peaks$chr, pos = trt_peaks$pos)
    if (is.null(covar_trt)) {
      idx <- which_lod_max[1]
    } else {
      r2 <- numeric(length = length(peak_mar))
      for (k in seq_along(peak_mar)) {
        geno <- maxmarg(probs = genoprob, map = pmap, chr = trt_peaks$chr[k], pos = trt_peaks$pos[k])
        cor_geno <- cor(geno, covar_trt, use = "pairwise.complete.obs")^2
        r2[k] <- max(cor_geno)
      }
      which_lod_max <- which_lod_max[r2 <= 0.80]
      idx <- which_lod_max[1]
    }

    if (is.na(idx)) break

    covar_mar <- find_marker(map = pmap, chr = trt_peaks$chr[idx], pos = trt_peaks$pos[idx])
    geno <- maxmarg(probs = genoprob, map = pmap, chr = trt_peaks$chr[idx], pos = trt_peaks$pos[idx])
    geno <- as.matrix(geno) - 2
    colnames(geno) <- covar_mar
    covar_trt <- cbind(covar_trt, geno)

    genoprob_use <- calc_genoprob(cross = drop_markers(cross = cross_qtl2_STxOxy, markers = colnames(covar_trt)), map = gmap, error_prob = 0.06,
                                  cores = 8)

    scan1_covar_out <- scan1(genoprobs = genoprob_use, pheno = phenos_lmm[,trt, drop = FALSE], kinship = kin,
                             addcovar = covar_trt, cores = 8)

    plot_scan1(x = scan1_covar_out, map = pmap, ylim = c(0, ymx))

    trt_peaks <- find_peaks(scan1_output = scan1_covar_out, map = pmap, threshold = scan1_ans_perm_thresholds["0.05", trt, drop = FALSE],
                            peakdrop = 2, drop = 2)
  }

  # return the list of markers selected in the multiple qtl mapping algorithm
  mqm_scan_ans[[trt]] <- data.frame(marker = colnames(covar_trt), row.names = NULL)


}





#













## Run chisq tests for categorical traits
cat_traits <- names(pheno_cat)[-1:-2]
cat_scan_out <- list()
for (trt in cat_traits) {
  dat <- pheno_cat[c("geno_id", "year", trt)]
  dat <- dat[!is.na(dat[[trt]]), ]
  dat <- subset(dat, geno_id %in% row.names(genoprob[[1]]))
  dat <- as.data.frame(dat)
  if (length(unique(dat$year)) == 1) {
    dat <- dat[c("geno_id", trt)]
  }

  # Index of rows in genoprobs to merge with dat
  idx <- match(x = dat$geno_id, table = row.names(genoprob[[1]]))
  dat <- dat[-1]

  form <- reformulate(termlabels = ".", response = trt)

  # fit the base model
  fit0 <- nnet::multinom(formula = form, data = dat, maxit = 500)

  scores_out <- NULL

  # Iterate over chromosome
  for (i in seq_along(genoprob)) {
    genoprob_i <- genoprob[[i]]
    scores_out_i <- matrix(0, nrow = dim(genoprob_i)[[3]], dimnames = list(dimnames(genoprob_i)[[3]], trt))
    for (j in seq_len(dim(genoprob_i)[[3]])) {
      # print(j)
      probs <- genoprob_i[idx,,j]
      dat_j <- cbind(dat, probs)
      out <- capture.output(fit <- nnet::multinom(formula = form, data = dat_j))
      # LRT
      lrt <- anova(fit,  fit0)
      # record the score
      scores_out_i[j,] <- -log10(lrt$`Pr(Chi)`[2])
    }
    scores_out <- rbind(scores_out, scores_out_i)
  }

  # Add to the list
  cat_scan_out_list[[trt]] <- list(score = scores_out, sample.size = length(unique(idx)))

}

# Bind columns
cat_scan_ans <- do.call(cbind, lapply(cat_scan_out_list, "[[", "score"))
attr(cat_scan_ans, "hsq") <- matrix(NA, nrow = length(genoprob), ncol = length(cat_scan_out_list), dimnames = list(seq_along(genoprob), names(cat_scan_out_list)))
attr(cat_scan_ans, "sample_size") <- do.call(c, lapply(cat_scan_out_list, "[[", "sample.size"))
class(cat_scan_ans) <- c("scan1", "matrix")

# Add to the scan output

scan1_ans <- lmm_scan_out_ols
scan1_ans <- cbind(scan1_ans, cat_scan_ans)
colnames(scan1_ans)






# Zoom in
plot(scan1_ans, map = gmap, lodcolumn = i, main = trt, chr = 2)




# Estimate additive and dominance effects
# chromosome 2 reflowering
i = 7
trt <- trait_names[i]
chr <- "5"
eff <- scan1coef(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]], reml = FALSE)
plot(eff, map = gmap, columns = 1:3, scan1_output = scan1_ans[,trt, drop = FALSE], main = trt)

eff <- scan1coef(genoprob[,chr], phenos_lmm[,trt], kinship = kin[[chr]], contrasts=cbind(mu=c(1,1,1), a=c(2, 1, 0), d=c(0, 1, 0)))
plot(eff, map = gmap, columns = 2:3)
last_coef <- unclass(eff)[nrow(eff),] # pull out last coefficients
for(j in seq(along=last_coef)) axis(side=4, at=last_coef[j], names(last_coef)[j], tick=FALSE, col.axis=col[j])


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




### Summarize QTL peaks -----------------------------------------------------

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

