# In this script I obtain the hit information from the exon GFF file
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
exon_gff_file = args[1]
out_tsv_file = args[2]

# Read in exon GFF file
gff_file <- read_tsv(exon_gff_file, col_names = FALSE,
    show_col_types = FALSE)

# Group exons by hit name and obtain hit information
if (nrow(gff_file) > 0){
    tsv_out <- gff_file %>%
    mutate(X11 = str_remove(X10, "_exon[0-9]+")) %>%
    group_by(X11) %>%
    summarize(
    hit_length = sum(X5 - X4),
    hit_strand = unique(X7),
    hit_chrom = unique(X1),
    hit_start = min(X4),
    hit_end = max(X5),
    hit_exon_count = n())

    write_tsv(tsv_out, out_tsv_file, col_names = FALSE)
}