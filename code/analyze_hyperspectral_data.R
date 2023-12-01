## MXO S1 Genetic mapping
##
## Examine hyperspectral data for the MxO S1 population
##

# Load packages, startup
library(tidyverse)
library(waves)
library(lme4)

source("startup.R")

# Cran directory
cran_dir <- proj_dir %>%
  str_split("/") %>%
  first() %>%
  {.[seq_len(str_which(., "Cranberry"))]} %>%
  paste0(., collapse = "/")
# Breeding directory
breeding_dir <- file.path(cran_dir, "Breeding")

# Read in the field book for this trial
trial_metadata <- read_csv(file.path(breeding_dir, "2021/Greenhouse/mxof2_2021_greenhouse_fieldbook_withParents.csv"))


# Read in the hyperspectral data ------------------------------------------

load(file.path(breeding_dir, "2021/Hyperspec/MXOF2_LeafImages/processed_scan_data",
               "mxo_2021_gh_hyperspecLeafScan_tidy.RData"))

# Remove data points intended only for measure leaf sizes
scan_data_leaf_df2 <- scan_data_leaf_df1 %>%
  split(.$unique_id) %>%
  modify_if(.x = ., .p = sapply(X = ., FUN = function(x) n_distinct(x$target) > 1), .f = ~filter(.x, target != "leafSize")) %>%
  bind_rows() %>%
  select(unique_id, scan_date, leaf_id = object_id, target, scan_data) %>%
  unnest(scan_data) %>%
  # Rename
  rename_at(vars(contains("nm")), ~paste0("X", str_remove(., "nm")))

# Remove previous data
rm(scan_data_leaf_df1)


# Summarize by leaf within sample
scan_data_leaf_df4 <- scan_data_leaf_df2 %>%
  group_by(unique_id, scan_date, leaf_id, target) %>%
  summarize_at(vars(starts_with("X", ignore.case = FALSE)), mean) %>%
  ungroup() %>%
  select(-target)

# Use the waves package to pre-process spectra
# Range of wavelengths
wl_range <- 593:1707
scan_data_preprocessed <- DoPreprocessing(df = scan_data_leaf_df4, preprocessing.method = 1, wavelengths = wl_range)


# Gather
scan_data_preprocessed_long <- scan_data_preprocessed %>%
  gather(wavelength, value, starts_with("X", ignore.case = FALSE)) %>%
  inner_join(select(trial_metadata, -note), .) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(vars(scan_date), as.factor)

# Iterate over wavelengths
hyperspec_blues <- scan_data_preprocessed_long %>%
  filter(wavelength %in% sample(levels(wavelength), 8)) %>%
  group_by(wavelength) %>%
  do({
    dat1 <- .

    # Model?
    fit <- lmer(value ~ -1 + individual + (1|scan_date), data = dat1)
    # fit1 <- lmer(value ~ individual + (1|individual:leaf_id) + (1|scan_date), data = dat1,
    #              control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))

    # Get blues and return
    summary(fit)$coef %>%
      as.data.frame() %>%
      rownames_to_column("individual") %>%
      mutate(individual = str_remove(individual, "individual")) %>%
      select(individual, estimate = Estimate)

  }) %>% ungroup()


# Create a matrix
hyperspec_blues_mat <- hyperspec_blues %>%
  spread(wavelength, estimate) %>%
  as.data.frame() %>%
  column_to_rownames("individual") %>%
  as.matrix() %>%
  scale()

# Create relationship matrix
H_mat <- tcrossprod(hyperspec_blues_mat) / ncol(hyperspec_blues_mat)

# PCA
H_pca <- prcomp(H_mat)
stdev <- H_pca$sdev / sum(H_pca$sdev)

H_pca_plot <- H_pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("individual") %>%
  as_tibble() %>%
  mutate(family = map_chr(str_split(individual, "-"), ~paste0(.x[1:2], collapse = "-"))) %>%
  filter(! individual %in% c("CNJ98-12-3", "CNJ98-317-48", "CNJ98-320-4", "CNJ98-7-5"))

ggplot(data = NULL, aes(x = PC1, y = PC2)) +
  geom_point(data = subset(H_pca_plot, family %in% names(subset(table(H_pca_plot$family), table(H_pca_plot$family) > 1))),
             aes(color = family), size = 2) +
  geom_point(data = subset(H_pca_plot, ! family %in% names(subset(table(H_pca_plot$family), table(H_pca_plot$family) > 1))),
             aes(shape = individual), size = 2)










