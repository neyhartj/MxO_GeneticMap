# MxO Genetic Map
#
# Visualize genome coverage
#

library(tidyverse)
library(Gviz)



# List coverage files
cov_dir <- file.path(results_dir, "variant_calling/genome_coverage/")
cov_files <- list.files(path = cov_dir, pattern = "coverage", full.names = TRUE)

track_list <- map(cov_files, ~{
  # Load bedGraph data
  track_data <- read_table(file = .x, col_names =  c("chr", "start", "end", "depth"))

  track_data <- track_data %>%
    mutate(chr = str_extract(string = track_data$chr, pattern = "chr[0-9]{1,}"),
           chr = paste0("chr", str_pad(str_extract(chr, "[0-9]{1,2}"), width = 2, side = "left", pad = "0"))) %>%
    filter(chr != "chrNA")

  track_data_gr <- GRanges(track_data)

  genome <- sub(pattern = "_alignment", replacement = "", x = str_extract(basename(.x), "[A-Za-z]{1,}_alignment"))
  DataTrack(track_data_gr, genome = genome, name = basename(.x))

})


plotTracks(track_list, chromosome = "chr12", col = "firebrick", )

# Plot as a line or area plot
track_data %>%
  filter(chr == "chr12") %>%
  ggplot(aes(x=start, y=depth)) +
  geom_area(fill="steelblue", alpha=1) +
  theme_classic() +
  labs(x="Genomic Position", y="Depth")

