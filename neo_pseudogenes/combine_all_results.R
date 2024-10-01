#!/usr/bin/env Rscript

# Load library
library(tidyverse)

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
working_dir <- args[1]

# Read in data
ap_groups <- read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/anchorpair_groups_all.tsv"))
pseudogenes <- read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/results_table_final.tsv"))
translocated <- read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/translocated_genes.tsv"), col_names = c("group_id", "level", "transloc_genes", "description"))
unannotated <- read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/results_table_unannotated_genes.tsv"))

# Combine anchor pair groups with translocated genes
df <- left_join(ap_groups, translocated, by = "group_id") %>%
  # Count number of translocated genes
  mutate(num_translocated = (str_count(transloc_genes, ",") + 1))

# Convert group table to long format (per gene) to be able to join with pseudogenes
gene_to_group <- ap_groups %>% separate_rows(genes, sep =",") %>%
  select(group_id, genes)
# Combine pseudogenes with gene_to_group
pseudogenes <- left_join(pseudogenes, gene_to_group, by = c("funct_paralogue" = "genes"))
# If hits correspond to functional paralogue that are of the same collinear group
# filter them out if they are on the same chromosome
pseudogenes <- pseudogenes %>% group_by(group_id, chrom) %>%
  slice_head(n = 1)

# Combine unannotated genes with gene_to_group
unannotated <- left_join(unannotated, gene_to_group, by = c("funct_paralogue" = "genes"))
# If hits correspond to functional paralogue that are of the same collinear group
# filter them out if they are on the same chromosome
unannotated <- unannotated %>% group_by(group_id, chrom) %>%
  slice_head(n = 1)

# Convert back to wide format (per group)
pseudogene_hits_per_group <- pseudogenes %>%
  mutate(group_id = as.factor(group_id)) %>%
  group_by(group_id) %>%
  summarise(pseudogenes = paste(pseudogene, collapse = ","),
            chrom_psgs = paste(chrom, collapse = ","),
            subgenomes_psgs = paste(unique(psg_subg), collapse = ",")) %>%
  ungroup()
# Convert back to wide format (per group)
unannotated_hits_per_group <- unannotated %>%
  mutate(group_id = as.factor(group_id)) %>%
  group_by(group_id) %>%
  summarise(unannotated = paste(pseudogene, collapse = ","),
            chrom_unann = paste(chrom, collapse = ","),
            subgenomes_unann = paste(unique(psg_subg), collapse = ",")) %>%
  ungroup()

# Combine all results
df <- left_join(df, pseudogene_hits_per_group, by = "group_id") %>%
  mutate(num_pseudogenes = (str_count(pseudogenes, ",") + 1)) %>%
  mutate(num_subg_pseudogenes = (str_count(subgenomes_psgs, ",") + 1))

df <- left_join(df, unannotated_hits_per_group, by = "group_id") %>%
  mutate(num_unannotated = (str_count(unannotated, ",") + 1)) %>%
  mutate(num_subg_unannotated = (str_count(subgenomes_unann, ",") + 1))

df <- df %>% select(group_id, num_aps = number_of_haplotypes,
              num_psgs = num_subg_pseudogenes,
              num_trans = num_translocated,
              num_unan = num_subg_unannotated,
              ap_subgenomes = haplotypes,
              psgs_subgenomes = subgenomes_psgs,
              unan_subgenomes = subgenomes_unann,
              genes, transloc_genes, pseudogenes, 
                     unannotated) %>%
  mutate(across(c(num_psgs, num_trans, num_unan), ~ replace_na(., 0)))

# write to file
write_tsv(df, paste0(working_dir, "/results/i-ADHoRe-run/table_level_all.tsv"))