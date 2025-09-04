## MxO Genetic Map
##
## Analyze genome synteny
##

library(tidyverse)

# Read the coords file (skip header lines if present)
coords <- read_table("analysis/GenomeAlignment/genome_alignment_filtered.coords", col_names = c("S1", "E1", "S2", "E2", "LENR", "LENQ", "IDY", "READFRAMER", "READFRAMEQ", "REFID", "QUERYID"))

# Optional: Filter for long, high-identity matches
coords_filtered <- coords %>%
  filter(LENR > 1000, IDY > 90) %>%
  mutate(orientation = ifelse((S1 - E1) * (S2 - E2) > 0, "forward", "inversion")) %>%
  filter(grepl(pattern = "chr", x = QUERYID))


ggplot(coords_filtered) +
  geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = orientation)) +
  scale_color_manual(values = c("forward" = "blue", "inversion" = "red")) +
  theme_minimal() +
  labs(
    title = "Synteny Plot",
    x = "Reference genome (bp)",
    y = "Query genome (bp)",
    color = "Orientation"
  )

ggplot(coords_filtered) +
  geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = orientation)) +
  facet_grid(. ~ QUERYID, scales = "free", space = "free") +
  theme_bw()
