#!/usr/bin/env Rscript
#

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage:\n")
  cat("Rscript merge_blast_map.R <blast_file> <genetic_map_file> <output_file>\n")
  quit(status = 1)
}

blast_file <- args[1]
genetic_map_file <- args[2]
output_file <- args[3]

# Read in the blast output
blast <- read.table(blast_file, header = FALSE, sep = "\t",
                    col.names = c("qseqid", "sseqid", "sstart", "send", "pident", "length", "evalue", "bitscore"))

# Get the best hit per marker
blast_best <- blast %>%
  group_by(qseqid) %>%
  top_n(n = 1, wt = bitscore) %>%
  ungroup()


# Only use markers that align uniquely
blast_best_unique <- blast_best %>%
  group_by(qseqid) %>%
  filter(n() == 1) %>%
  ungroup()

# Calculate the midpoint of each contextual sequence
blast_best_unique1 <- blast_best_unique %>%
  mutate(position = floor((sstart + send) / 2)) %>%
  # Convert the qseqid to data in the bed file
  separate(qseqid, c("scaffold", "range"), sep = ":", remove = FALSE) %>%
  separate(range, c("start", "end"), sep = "-") %>%
  mutate_at(vars(start, end), as.numeric) %>%
  # Create the marker name
  mutate(marker.name = paste0(scaffold, "_", ((start + end) / 2) + 1))

# Read in the genetic map
genetic_map <- read.csv(genetic_map_file, col.names = c("chrom", "marker.name", "locus.name", "type", "position"))

# Merge with the blast results
genetic_map_merged <- genetic_map %>%
  inner_join(., blast_best_unique1, by = "marker.name") %>%
  # Get the columns for allmaps
  select(scaffold_id = sseqid, scaffold_position = position.y, chrom, genetic_position = position.x) %>%
  mutate(chrom = paste0("chr", sub("Cranberry-Composite_map-F1.LG", "", chrom)),
         genetic_position = as.numeric(sub(" cM", "", genetic_position)))


# Save the merged file
write.csv(genetic_map_merged, output_file, row.names = FALSE, quote = FALSE)









