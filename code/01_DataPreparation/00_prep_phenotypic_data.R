## MxO Genetic map
##
## Compile phenotypic data
##

# Load packages
library(tidyverse)
library(readxl)
library(lme4)


# Read in the phenotypic data ---------------------------------------------

pheno <- read_excel(path = "data_raw/mxof2_raw_phenotypic_data_compiled.xlsx", sheet = "data_tidy") %>%
  # Edit columns
  mutate_at(vars(DaysToSecFlw, PowderyMildewSeverity, TotalSecondFrtWt), parse_number) %>%
  mutate_at(vars(year, bench), as.factor) %>%
  mutate_at(vars(row, col), list(id = as.factor))


# Visualize
pheno %>%
  subset(!is.na(bench)) %>%
  ggplot(aes(x = row, y = col, fill = DaysToFlw)) +
  geom_tile() +
  facet_grid(year ~ bench)

pheno %>%
  subset(!is.na(bench)) %>%
  ggplot(aes(x = row, y = col, fill = PowderyMildewSeverity)) +
  geom_tile() +
  facet_grid(year ~ bench)

pheno %>%
  subset(!is.na(bench)) %>%
  ggplot(aes(x = row, y = col, fill = AvgFirstFrtWt)) +
  geom_tile() +
  facet_grid(year ~ bench)



# Model -------------------------------------------------------------------

# Traits to model
traits_model <- c("DaysToFlw", "DaysToSecFlw", "FlwInterval", "AvgFirstFrtWt")

# Calculate BLUEs
pheno_blues_list <- list()
for (trt in traits_model) {
  formula <- reformulate(termlabels = c("individual", "(1|year)"), response = trt, intercept = FALSE)
  fit <- lmer(formula = formula, data = pheno)
  blues <- fixef(fit)
  blues_df <- data.frame(geno_id = sub(pattern = "individual", replacement = "", x = names(blues)), BLUE = blues, row.names = NULL)
  names(blues_df) <- c("geno_id", trt)
  pheno_blues_list[[trt]] <- blues_df

}

pheno_blues_df <- reduce(pheno_blues_list, full_join) %>%
  mutate(FlwIntervalBLUE = DaysToSecFlw - DaysToFlw)


# Summarize other traits
traits_summ <- c("PollenViability", "PowderyMildewSeverity", "AverageFlowersPerUpright", "NumberFloweringUprights", "Prolificacy")

pheno_summ_list <- list()
for (trt in traits_summ) {
  formula <- reformulate(termlabels = c("individual"), response = trt)
  blues_df <- aggregate(formula, pheno, FUN = mean, na.rm = TRUE)
  names(blues_df) <- c("geno_id", trt)
  pheno_summ_list[[trt]] <- blues_df

}

pheno_summ_df <- reduce(pheno_summ_list, full_join)

pheno_df <- full_join(pheno_blues_df, pheno_summ_df) %>%
  rename(id = geno_id)

# Save
write_csv(x = pheno_df, file = "data/mxo_all_family_pheno.csv")












