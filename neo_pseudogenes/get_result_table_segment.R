
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
start_region = args[1]
location = args[2]
scores = args[3]
out_result = args[4]
gff_file = args[5]
output = args[6]
multiplicon = args[7]


# get filtered result tables that contains perc id
scores_result <- read_tsv(scores, 
                   col_names = c("hit", "gene", "hit_aln_start", "hit_aln_end", 
                                 "hit_aln_length", "gene_length", "strand",
                                 "hit_strand", "perc_id"),
                   show_col_types = FALSE)

# get coordinates of start regions
start_region <- read_tsv(start_region, 
                          col_names = 
                            c("start"),
                          show_col_types = FALSE)

# get coordinates of hits and change local coordinates to global coordinates
location_result <- read_tsv(location, 
                             col_names = 
                               c("hit", "hit_length", "gene_strand", "hit_chrom",
                                 "hit_start", "hit_end", "hit_exon_count"),
                             show_col_types = FALSE) 

if (exists("location_result") && nrow(location_result) > 0) {
  location_result <- location_result %>%
  mutate(hit_start = hit_start + start_region$start,
         hit_end = hit_end + start_region$start)
}

# write gff files and update coordinates from local to global
gff_table <- read_tsv(gff_file, 
                      col_names = 
                        c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "hit_name"),
                      show_col_types = FALSE) %>%
                        mutate(start = start + start_region$start,
                                end = end + start_region$start) %>%
                        mutate(hit = str_extract(hit_name, ".*(?=_exon[0-9]+)")) %>%
                        group_by(hit) %>%
                        arrange(start, .by_group = TRUE)
write_tsv(gff_table, gff_file, col_names = FALSE)

# create final result table
result <- left_join(scores_result, location_result, by = "hit") %>%
  select(hit, hit_chrom, gene, perc_id, hit_start, hit_end, hit_length, hit_strand,
         hit_aln_start, hit_aln_end, hit_aln_length, gene_length)
write_tsv(result, out_result, append = TRUE)