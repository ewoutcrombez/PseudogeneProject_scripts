#!/usr/bin/env Rscript
library(tidyverse)
library(igraph)
library(ggthemes)
library(ggtext)

args = commandArgs(trailingOnly=TRUE)
working_dir <- args[1]

# dataframe containing all homologous gene pairs in collinear segments
df <- 
  read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/hom_pairs_all.tsv"),
           col_names = c("hom1", "hom2", "multiplicon"))

all_genes <- 
  read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/all_genes.txt"),
           col_names = "gene")

# Function to sort first two columns alphabetically per row
df_sorted <- df %>%
  rowwise() %>%
  mutate(
    sorted_cols = list(sort(c(hom1, hom2)))
  ) %>%
  mutate(
    hom1 = sorted_cols[[1]],
    hom2 = sorted_cols[[2]]
  ) %>%
  select(-sorted_cols) %>%  # Remove helper column
  ungroup() %>%
  distinct(hom1, hom2, .keep_all = TRUE)

# Function to get anchor pair groups
get_ap_groups <- function(){
  ap_groups <- df_sorted %>%
    select(hom1, hom2) %>%
    graph_from_data_frame(directed = FALSE) %>%
    components()
  ap_groups_df <- tibble(
    group_id = ap_groups$membership,
    gene = names(ap_groups$membership)
  ) %>%
    mutate(haplotype = str_sub(gene, 1, 1)) %>%
    group_by(group_id) %>%
    summarize(genes = paste(gene, collapse = ","),
              haplotypes = paste(sort(unique(haplotype)), collapse = ","),
              number_of_haplotypes = length(unique(haplotype)))
}

# Combine all anchor pairs that have the same genes
ap_groups <- get_ap_groups()

# Separate genes in each group
levels <- ap_groups %>% 
  separate_rows(genes, sep = ",")
# Get genes that are not part of any anchor pair and are thus collinear level 1 (only in one sub-genome)
level1 <- all_genes %>% filter(!gene %in% levels$genes)
level1_df <- tibble(group_id = paste0("l", seq(1, level1$gene %>% length())), 
       genes = level1$gene,
       haplotypes = level1$gene %>% str_sub(start = 1, end = 1),
       number_of_haplotypes = rep("1", level1$gene %>% length()))
levels <- rbind(levels, level1_df)

# Get the count of genes per collinearity level
levels_count <- levels %>% group_by(number_of_haplotypes) %>%
  count()

# Plot the gene retention levels
plot <- ggplot(levels, aes(x = number_of_haplotypes)) +
  geom_bar() +
  theme_clean(base_size = 24) +
  xlab("Gene retention level on sub-genomes") +
  ylab("") +
  geom_richtext(data = levels_count, aes(y = n + 500, label = n),
                        size = 6, fill = NA, label.colour = NA)
ggsave(paste0(working_dir, "/results/i-ADHoRe-run/gene_level_collinearity.svg"), 
       plot,
       width = 17.5, height = 22.25, units = "cm")

# Write the different gene retention levels to file (including the level 1 genes)
result_groups <- rbind(ap_groups, level1_df)
write_tsv(result_groups, paste0(working_dir, "/results/i-ADHoRe-run/anchorpair_groups_all.tsv"))