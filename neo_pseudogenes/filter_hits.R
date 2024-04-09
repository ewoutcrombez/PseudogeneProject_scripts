# In this script I filter hits
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
tsv_file = args[1]
gene_tsv_file = args[2]
exon_gff_file = args[3]
out_tsv_file = args[4]

# Read in TSV file
res_tsv <- read_tsv(tsv_file, col_names =
  c("query_id", "target_id", "q_strand", "t_strand",
    "q_length", "t_length", "t_aln_start", "t_aln_end", "t_aln_length",
    "q_aln_start", "q_aln_end", "q_aln_length", "perc_id",
    "raw_score", "total_equivalenced", "identical",
    "similarity", "mismatches"),
    show_col_types = FALSE)

# Filter based on query alignment length and percent identity
head(res_tsv)
res_tsv <- res_tsv %>%
  filter(q_length > 100 | q_length > (0.5 * t_length))
  
if (nrow(res_tsv) > 0){
  # only select best hit (with highest percent identity) for further analysis
  res_tsv <- res_tsv %>%
    group_by(query_id) %>%
    summarize(
      target_id = unique(target_id),
      hit_aln_start = min(q_aln_start),
      hit_aln_end = max(q_aln_end),
      hit_aln_length = max(q_aln_end) - min(q_aln_start),
      t_length = unique(t_length),
      q_strand = unique(q_strand),
      t_strand = unique(t_strand),
      avg_perc_id = (sum(
        perc_id * q_aln_length)) / (sum(q_aln_length))) %>%
    filter(hit_aln_length > 0.2 * t_length) %>%
    filter(avg_perc_id > 20) %>%
    ungroup()}
  
if (nrow(res_tsv) > 0){
  res_tsv <- res_tsv %>%
    group_by(target_id) %>%
    filter(hit_aln_length == max(hit_aln_length))
  write_tsv(res_tsv,
            out_tsv_file, col_names = FALSE)
  # remove target_id from query_id for further visualization
  res_tsv_only_hitnumber <- res_tsv %>%
    mutate(query_id = str_remove(query_id, target_id))
  write_tsv(res_tsv_only_hitnumber,
    str_replace(out_tsv_file, ".tsv", "_hitname.tsv"),
    col_names = FALSE)
}

# Filter gene TSV file
res_gene_bed <- read_tsv(gene_tsv_file, col_names = FALSE,
  show_col_types = FALSE)
bed_gene <- res_gene_bed %>%
  filter(X1 %in% res_tsv$query_id)
write_tsv(bed_gene,
  gene_tsv_file, col_names = FALSE)
# remove target_id from query_id for further visualization
bed_gene_only_hitnumber <- bed_gene %>%
  mutate(X1 = str_extract(X1, "_hit.*"))
write_tsv(bed_gene_only_hitnumber,
  str_replace(gene_tsv_file, ".tsv", "_hitname.tsv"),
  col_names = FALSE)

# Filter exon GFF file
res_exon_gff <- read_tsv(exon_gff_file, col_names = FALSE,
  show_col_types = FALSE)
gff_exons <- res_exon_gff %>%
  filter(str_detect(X10, str_c(res_tsv$query_id, collapse = "|")))
write_tsv(gff_exons,
  exon_gff_file, col_names = FALSE)