## MxO Genetic map
##
## Analyze phenotypic data for input to QTL mapping
##

# Load packages
library(tidyverse)
library(readxl)
library(lme4)
library(bestNormalize)
library(emmeans)
library(lmerTest)

# Get the population metadata
pop_metadata <- read_csv("data/population_metadata.csv")

col <- c("slateblue", "violetred", "green3")


# Read in the phenotypic data ---------------------------------------------

pheno <- read_excel(path = "data_raw/mxof2_raw_phenotypic_data_compiled.xlsx", sheet = "data_tidy", na = c("", "NA")) %>%
  # Edit columns
  rename(geno_id = individual) %>%
  mutate_at(vars(year, bench, geno_id), as.factor) %>%
  mutate_at(vars(row, col), list(id = as.factor))


# Adjust some traits
pheno1 <- pheno %>%
  rename(DaysToFirstFlw = DaysToFlw, DaysToSecondFlw = DaysToSecFlw) %>%
  mutate(FlwInterval = DaysToSecondFlw - DaysToFirstFlw,
         ReFlw = case_when(
           !is.na(DaysToFirstFlw) & !is.na(DaysToSecondFlw) ~ TRUE,
           !is.na(DaysToFirstFlw) & is.na(DaysToSecondFlw) ~ FALSE,
           is.na(DaysToFirstFlw) & !is.na(DaysToSecondFlw) ~ FALSE,
           .default = as.logical(NA)
         ),
         ReFlwBin = as.numeric(ReFlw))

# # Data checking
# # If you flowered but didn't have fruit, set fruit count to 0
# pheno2 <- pheno1 %>%
#   mutate(FirstFrtCt = ifelse(!is.na(DaysToFirstFlw) & is.na(FirstFrtCt), 0, FirstFrtCt),
#          SecondFrtCt = ifelse(!is.na(DaysToSecondFlw) & is.na(SecondFrtCt), 0, SecondFrtCt))
pheno2 <- pheno1

# Create alternative flowering traits
# 1. DaysToFirstFlw_Overall - days to first flower, regardless of round
# 2. DaysToFlw_Once - days to first flower, if flowered one time
pheno3 <- pheno2 %>%
  mutate(DaysToFirstFlw_Overall = pmin(DaysToFirstFlw, DaysToSecondFlw, na.rm = TRUE),
         DaysToFlw_Once = ifelse(ReFlw, as.numeric(NA), DaysToFirstFlw_Overall),
         FloweringType = case_when(
           ReFlw ~ "reflowering",
           !is.na(DaysToFirstFlw) & is.na(DaysToSecondFlw) ~ "early_flowering",
           is.na(DaysToFirstFlw) & !is.na(DaysToSecondFlw) ~ "late_flowering"
         ))





# Trait names
all_traits <- grep(pattern = "^[A-Z]", x = names(pheno3), ignore.case = FALSE, value = TRUE)
all_traits <- setdiff(all_traits, "ReFlw")

# Classify traits
all_traits_type <- sapply(pheno3[all_traits], class)
all_traits_type["ReFlwBin"] <- "binary"

# Remove extra comments and such
# Compute the number of non-NA values for each trait per year
pheno4 <- pheno3 %>%
  select(-note, -pollen_staining_notes)

pheno_obs_summ <- sapply(all_traits, FUN = function(trt) tapply(X = pheno4[[trt]], INDEX = pheno4$year, FUN = function(x) sum(!is.na(x))) )
pheno_obs_summ


# Traits to be analyzed using BLUEs
model_traits_blues <- colMeans(pheno_obs_summ == 0) == 0
model_traits_blues <- names(model_traits_blues)[model_traits_blues]
model_traits_blues <- intersect(model_traits_blues, names(all_traits_type)[all_traits_type == "numeric"])

# Traits to analyze by simple means
model_traits_avg <- colMeans(pheno_obs_summ == 0) > 0
model_traits_avg <- names(model_traits_avg)[model_traits_avg]
model_traits_avg <- intersect(model_traits_avg, names(all_traits_type)[all_traits_type == "numeric"])

# Traits to be included as per-year because of character or binary nature
traits_per_year <- names(all_traits_type)[all_traits_type != "numeric"]


# Vizualize phenotypic values
pheno4_numeric_tidy <- pheno4 %>%
  select(geno_id, year, any_of(c(model_traits_blues, model_traits_avg))) %>%
  gather(trait, value, -geno_id, -year)

g_pheno_dist <- pheno4_numeric_tidy %>%
  ggplot(aes(x = value, fill = year)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ trait, ncol = 4, scales = "free")
g_pheno_dist
ggsave(filename = "phenotypic_values_dist.jpg", plot = g_pheno_dist, path = fig_dir, height = 8, width = 10, dpi = 500)


# Visualize counts of categorical data
pheno4_cat_tidy <- pheno4 %>%
  select(geno_id, year, any_of(traits_per_year)) %>%
  gather(trait, category, -geno_id, -year) %>%
  group_by(year, trait, category) %>%
  summarize(n = n()) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  subset(p < 1)

g_pheno_cat <- pheno4_cat_tidy %>%
  ggplot(aes(x = category, y = n, fill = year)) +
  geom_col(position = position_dodge(0.9)) +
  facet_wrap(~ trait, ncol = 2, scales = "free")
g_pheno_cat
ggsave(filename = "phenotypic_categorical_count.jpg", plot = g_pheno_cat, path = fig_dir, height = 6, width = 8, dpi = 500)





# Calculate genotype means  -------------------------------------------------------------------

## Compute BLUEs

# Calculate BLUEs
pheno_blues_list <- list()

for (trt in c(model_traits_blues)) {
  pheno4_tmp <- pheno4[!is.na(pheno4[[trt]]), ]
  # Normalize the data using quantiles
  norm_out <- bestNormalize(x = pheno4_tmp[[trt]], standardize = FALSE)
  trt_norm <- paste0(trt, ".norm")
  pheno4_tmp[[trt_norm]] <- norm_out$x.t

  fixed <- c("geno_id", "year")
  rand <- c("(1|bench:year)", "(1|row_id:bench:year)", "(1|col_id:bench:year)")

  formula <- reformulate(termlabels = c(fixed, rand), response = trt)
  fit <- lmer(formula = formula, data = pheno4_tmp, REML = FALSE)

  # Reduce variables
  fit_red <- ranova(fit)
  reduce <- any(fit_red$`Pr(>Chisq)` > 0.05)
  while (reduce) {
    rand <- rand[-(which.max(fit_red$`Pr(>Chisq)`) - 1)]
    if (length(rand) == 0) break
    formula <- reformulate(termlabels = c(fixed, rand), response = trt)
    fit <- lmer(formula = formula, data = pheno4_tmp, REML = FALSE, na.action = na.pass)
    # Reduce variables
    fit_red <- ranova(fit)
    reduce <- any(fit_red$`Pr(>Chisq)` > 0.05, na.rm = TRUE)
  }

  if (length(rand) == 0) {
    fit <- lm(formula = reformulate(termlabels = fixed, response = trt), data = pheno4_tmp)
  } else {
    fit <- lmer(formula = reformulate(termlabels = c(fixed, rand), response = trt, intercept = FALSE), data = pheno4_tmp)
  }

  blues <- emmeans(fit, specs = "geno_id")
  blues_df <- as.data.frame(blues)[,c("geno_id", "emmean")]
  names(blues_df)[2] <- trt


  # # Get genotype blues
  # blues_df <- data.frame(geno_id = row.names(blues), BLUE = blues[[1]], row.names = NULL)
  # blues_df <- subset(blues_df, grepl(pattern = "^geno_id", x = blues_df$geno_id))
  # blues_df$geno_id <- sub(pattern = "geno_id", replacement = "", x = blues_df$geno_id)
  # names(blues_df)[2] <- trt

  pheno_blues_list[[trt]] <- blues_df

}

pheno_blues_df <- reduce(pheno_blues_list, full_join)





## Compute means for single-year data
pheno_avg_list <- list()

for (trt in model_traits_avg) {
  pheno4_tmp <- pheno4[!is.na(pheno4[[trt]]), ]
  formula <- reformulate(termlabels = c("geno_id"), response = trt)
  blues_df <- aggregate(formula, pheno, FUN = mean, na.rm = TRUE)
  names(blues_df) <- c("geno_id", trt)
  pheno_avg_list[[trt]] <- blues_df
}

pheno_avg_df <- reduce(pheno_avg_list, full_join)

# Combine all numeric traits
pheno_geno_blues_df <- full_join(pheno_blues_df, pheno_avg_df)


# Organize categorical traits
pheno_cat_df <- pheno4[c("geno_id", "year", traits_per_year)]

# Output data to CSVs
write_csv(x = pheno_geno_blues_df, file = file.path(data_dir, "mxo_all_family_pheno_blues.csv"))
write_csv(x = pheno_cat_df, file = file.path(data_dir, "mxo_all_family_pheno_categorical.csv"))




