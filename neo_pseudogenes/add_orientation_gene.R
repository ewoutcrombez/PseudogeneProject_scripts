# add orientation (strand information) to gene TSV file
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
region_txt_file = args[1]
bed_file = args[2]

# Read in files
region <- read_tsv(region_txt_file, col_names = FALSE,
    show_col_types = FALSE)
bed <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "name"),
    show_col_types = FALSE)

# Add orientation to BED file from region txt file
inner_join(bed, region, by = c("name" = "X1")) %>%
    select(chrom, start, end, name, X4) %>%
    write_tsv(bed_file, col_names = FALSE)
