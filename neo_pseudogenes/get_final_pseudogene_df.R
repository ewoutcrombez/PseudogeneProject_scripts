#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
working_dir <- args[1]

# load in data
## Exonerate results
exonerate_res <- 
  read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/result_table.tsv"),
           col_names = c("pseudogene", "chrom", "funct_paralogue", "perc_id",
                         "start", "end", "pseudogene_length", "strand", "aln_start",
                         "aln_end", "aln_length", "gene_length", "multiplicon")) %>%
  # only keep one entry of a pseudogene in a multiplicon (i.e. remove duplicate entries)
  distinct(pseudogene, multiplicon, .keep_all = TRUE) %>%
  # only keep one entry of exactly the same pseudogene which matches the same gene
  distinct(chrom, start, end, funct_paralogue, .keep_all = TRUE) %>%
  # for each functional paralogue, only keep the longest hit on a specified chromosome
  group_by(funct_paralogue, chrom) %>%
  slice_max(pseudogene_length, n = 1) %>%
  # if it is the same, keep the one with highest perc id
  slice_max(perc_id, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # if a hit matches multiple functional paralogues, only keep the longest one
  # and if same length, the one with highest perc id
  mutate(hit_name = str_extract(pseudogene, "hit(merged)?[0-9]+")) %>%
  mutate(multiplicon = as.character(multiplicon)) %>%
  mutate(hit_name = paste(hit_name, multiplicon, sep = "-")) %>%
  group_by(hit_name) %>%
  slice_max(pseudogene_length, n = 1) %>%
  slice_max(perc_id, n = 1, with_ties = FALSE) %>%
  ungroup()
  
## MACSE results
MACSE_res <-
  read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/result_table_MACSE.tsv"),
           col_names = c("pseudogene",  "pseudogene_length", "funct_paralogue",
                         "gene_length", "perc_id_cds", "fraction_cds", "length_aln_cds",
                         "perc_id_prot", "fraction_prot", "length_aln_prot",
                         "num_fs", "num_stop", "start_codon", "stop_codon",
                         "multiplicon")) %>%
  # only keep one entry of a pseudogene in a multiplicon (i.e. remove duplicates)
  distinct(pseudogene, multiplicon, .keep_all = TRUE) %>%
  mutate(multiplicon = as.character(multiplicon))

# Combine exonerate and MACSE results
res <- left_join(exonerate_res, MACSE_res, by = c("pseudogene", "multiplicon",
                                                  "funct_paralogue")) %>%
  mutate(psg_subg_vs_g_subg = paste(str_sub(chrom, 1, 1),
                                    str_sub(funct_paralogue, 1, 1),
                                    sep = "-")) %>%
  mutate(psg_subg = str_sub(chrom, 1,1)) %>%
  mutate(psg_subg_vs_g_subg = as_factor(psg_subg_vs_g_subg)) %>%
  mutate(subg_vs_subg = paste(str_sub(pmin(chrom, funct_paralogue), 1, 1),
                              str_sub(pmax(chrom, funct_paralogue), 1, 1),
                              sep = "-")) %>%
  mutate(subg_vs_subg = as.factor(subg_vs_subg))

# Filter out pseudogenes that are unannotated genes (no detrimental mutations)
res_unannot <- res %>%
  filter(num_fs == 0 & num_stop == 1 & start_codon != "absent" &
           stop_codon != "absent" & fraction_prot >= 0.95 &
           (perc_id_prot >= 0.95 | perc_id >= 95))

res_filt <- anti_join(res, res_unannot)

# Write out results
write_tsv(res_filt, paste0(working_dir, "/results/i-ADHoRe-run/results_table_final.tsv"))
# Write out unannotated genes
write_tsv(res_unannot, paste0(working_dir, "/results/i-ADHoRe-run/results_table_unannotated_genes.tsv"))