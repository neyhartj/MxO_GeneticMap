<!-- README.md is generated from README.Rmd. Please edit that file -->

# MxO_GeneticMap

## Variant Calling Pipeline

Steps for completing the variant calling pipeline

-   [x]  Create indexed references for Stevens and NJ96-20

-   [ ] Align reads to both genomes using BWA

-   [ ] Remove duplicate alignments (important because this is targeted sequencing) and merge .bam files

-   [ ] Use freebayes to call SNPs separately based on alignments

-   [ ] Filter SNPs for quality

-   [ ] Construct genetic maps using SNPs from either alignment
