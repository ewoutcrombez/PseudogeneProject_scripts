# In this script lonely genes that are actually tandem repeats of a paired gene
# will be removed from lonely gene list.
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
working_dir = args[1]
multiplicon = args[2]
lonely_file = args[3]
chrom1 = args[4]
chrom2 = args[5]
region1 = args[6]
region2 =args[7]
#region2_rev = args[8]
species1 = args[9]
species2 = args[10]

setwd(working_dir)
getwd()

genes <- read_tsv("results/i-ADHoRe-run2/genes.txt",
                  show_col_types = FALSE)

lonely_g <- read_tsv(lonely_file, col_names = "id",
                     show_col_types = FALSE)

paired_g <- read_tsv(str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/paired_genes.txt"), col_names = "id",
                     show_col_types = FALSE)

multiplicon_t <- read_tsv(str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/multiplicon_{multiplicon}.txt"), 
                          col_names = FALSE,
                          show_col_types = FALSE)

# Inner join of lonely and paired genes should be empty
print("If the dataframe below is not empty, the lonely genes are not (all) lonely:")
inner_join(lonely_g, paired_g) # OK!

# FILTER OUT TANDEM REPEATS
# Get lonely genes that are tandem and have a tandem_representative that is
# paired
tandem_rep <- right_join(genes, lonely_g, by = "id") %>% 
  filter(is_tandem == -1) %>%
  select(id, tandem_representative) %>% unique() %>%
  mutate(paired_tandem_rep = tandem_representative %in% paired_g$id) %>%
  filter(paired_tandem_rep == TRUE)

write_tsv(tandem_rep %>% select(id), str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/tandem_genes.txt"), col_names = FALSE)

# Get the homologue of the tandem gene that is on the other segment by 
# selecting the gene that is homologous to the tandem representative
# may be in column 4 or 5 of the multiplicon file
tandem_with_homologue <- tandem_rep %>%
  left_join(multiplicon_t, by = c("tandem_representative" = "X5")) %>%
  select(-X1, -X2, -X3, -X6, -X7, -X8) %>%
  left_join(multiplicon_t, by = c("tandem_representative" = "X4")) %>%
  select(-X1, -X2, -X3, -X6, -X7, -X8) %>%
  mutate(matched_tandem = ifelse(is.na(X4), X5, X4)) %>%
  select(-X4, -X5) %>% select(1,4)
write_tsv(tandem_with_homologue, 
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/tandem_withhom.txt"), 
          col_names = FALSE)

# Filter out genes that are tandem repeats of a paired gene
lonely_g <- lonely_g %>% filter(!id %in% tandem_rep$id)

# FILTER OUT LONELY GENES THAT PAIR WITH EACH OTHER
# Get blast results of lonely genes (not the ones assigned as tandem)
blastp <- read_tsv(str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/tmp_blast"),
                   col_names = FALSE, show_col_types = FALSE) %>% 
  filter(!X1 %in% tandem_rep$id | !X2 %in% tandem_rep$id)
head(blastp)

# Select pairs of genes that are both lonely genes and on different segments
blastp_f <- blastp %>% filter(X1 %in% lonely_g$id & X2 %in% lonely_g$id) %>%
  filter(X1 != X2) %>% filter(str_extract(X1, "^[^-]+") != str_extract(X2, "^[^-]+")) %>%
  filter(X3 >= 60 & X4 >= 0.7*X13) # additional filter %identity >= 60 and alignment length >= 0.7xquery_length
write_tsv(blastp_f, 
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_pairs.txt"),
          col_names = FALSE)
head(blastp_f)

# Select the genes that are part of such pairs
matched_genes <- lonely_g %>% filter(id %in% blastp_f$X1 | id %in% blastp_f$X2)
write_tsv(matched_genes, str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_genes.txt"), 
          col_names = FALSE)

# Filter out the genes that are part of "lonely gene pair"
lonely_g <- lonely_g %>% filter(!id %in% blastp_f$X1 & !id %in% blastp_f$X2)
head(lonely_g)
# FILTER OUT LONELY GENES THAT PAIR WITH OTHER ON SAME CHROMOSOME
# Get gff information
gff <- read_tsv(str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/tmp_gff"), 
                col_names = c("chrom", "source", "type", "start", "end", "x",
                              "or", "y", "gene"),
                show_col_types = FALSE)

# Get the lonely genes that are part of the haplotypes (species)
blastp_m <- blastp %>% filter(X1 %in% lonely_g$id | X2 %in% lonely_g$id) %>%
  filter(X1 != X2) %>%
  filter((startsWith(X2, species2) & startsWith(X1, species1)) | 
           (startsWith(X2, species1) & startsWith(X1, species2)))
# Get the ones that are on the same chromosomes
blastp_m <- inner_join(blastp_m, gff, by = c("X1" = "gene")) %>% select(1:15) %>%
  rename(chrom_X1 = chrom) %>% inner_join(gff, by = c("X2" = "gene")) %>%
  select(1:16) %>%
  rename(chrom_X2 = chrom) %>%
  filter((chrom_X1 == chrom1 | chrom_X1 == chrom2) & 
           (chrom_X2 == chrom2 | chrom_X2 == chrom1)) %>%
  filter(X3 >= 60 & X4 >= 0.7*X13) # additional filter %identity >= 60 and alignment length >= 0.7xquery_length
write_tsv(blastp_m, str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_chrom_table.txt"),
          col_names = FALSE)

# Check which genes are part of multiplicon
match_in_multiplicon <-
  blastp_m %>% filter(X1 %in% c(multiplicon_t$X4, multiplicon_t$X5) |
                            X2 %in% c(multiplicon_t$X4, multiplicon_t$X5))
write_tsv(match_in_multiplicon %>% select(X1, X2), 
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_in_mult_pairs.txt"),
          col_names = FALSE)
match_g_in_mult <- c(match_in_multiplicon$X1, match_in_multiplicon$X2) %>%
  as_tibble() %>% filter(value %in% lonely_g$id)
lonely_g <- lonely_g %>% filter(!id %in% match_g_in_mult$value)

# Get genes that have match on chromosome but not in multiplicon
match_in_chrom <- c(blastp_m$X1, blastp_m$X2) %>% unique() %>% as_tibble() %>%
  filter(value %in% lonely_g$id) %>%
  filter(!value %in% c(match_in_multiplicon$X1, match_in_multiplicon$X2))
write_tsv(match_in_chrom, str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_chrom_genes.txt"),
          col_names = FALSE)
lonely_g <- lonely_g %>% filter(!id %in% match_in_chrom$value)

# FILTER OUT LONELY GENES THAT PAIR WITH OTHER ON DIFFERENT CHROMOSOME
# (BUT BETWEEN SAME HAPLOTYPES)
blastp_l <- blastp %>% filter(X1 %in% lonely_g$id) %>% 
  filter(X1 != X2) %>%
  filter((startsWith(X2, species2) & startsWith(X1, species1)) | 
           (startsWith(X2, species1) & startsWith(X1, species2)))
head(blastp_l)
nrow(blastp_l)
if (nrow(blastp_l) > 0) {
blastp_l_max <- blastp_l %>%
  group_by(X1) %>% filter(X11 == min(X11)) %>%
  filter(X12 == max(X12)) %>%
  filter(X3 >= 60 & X4 >= 0.7*X13) %>% # additional filter %identity >= 60 and alignment length >= 0.7xquery_length
  rename("lonely_gene" = "X1", "id_match" = "X2") %>%
  inner_join(gff, by = c("id_match" = "gene"))
}

if (nrow(blastp_l) > 0) {
write_tsv(blastp_l_max, str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_scat_table.txt"),
          col_names = FALSE)
match_scattered <- blastp_l_max %>% select(lonely_gene)
write_tsv(match_scattered,
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/matched_scat_genes.txt"),
          col_names = FALSE)
# Save remaining "lonely" genes
lonely_g <- lonely_g %>% filter(!id %in% match_scattered$lonely_gene)
# Write remaining "lonely genes" to file
write_tsv(lonely_g, lonely_file, col_names = FALSE)
} else {
  match_scattered <- tibble(lonely_gene = NA)
}

# Get lonely genes per segment
region_1 <- read_tsv(str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/region_{species1}_{chrom1}.txt"),
              col_names = FALSE,
              show_col_types = FALSE)
lonely_g_r1 <- lonely_g %>% filter(id %in% region_1$X1)
write_tsv(lonely_g_r1, 
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/lonely_g_{species1}_{chrom1}.txt"),
          col_names = FALSE)
region_2 <- read_tsv(str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/region_{species2}_{chrom2}.txt"), 
                     col_names = FALSE, show_col_types = FALSE)
lonely_g_r2 <- lonely_g %>% filter(id %in% region_2$X1)
write_tsv(lonely_g_r2, 
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/lonely_g_{species2}_{chrom2}.txt"),
          col_names = FALSE)

# Get best hits for lonely genes
if (nrow(blastp_l %>% filter(X1 %in% lonely_g$id)) > 0) {
blastp_lonely_g <- blastp_l %>% filter(X1 %in% lonely_g$id) %>%
  group_by(X1) %>% filter(X11 == min(X11)) %>% filter(X12 == max(X12)) %>%
  arrange(desc(X12))
write_tsv(blastp_lonely_g, 
          str_glue("results/i-ADHoRe-run2/processing/{multiplicon}/lonely_genes_blastp.txt"))
}

# Color genes in region files that will be used for visualization according to
# whether they are tandem repeats, matches from each other or have a match
# in the same chromosomes
color_genes <- function(file, lonely_g, tandem_rep, matched_genes, match_in_chrom){
  region <- read_tsv(file, 
                     col_names = c("name", "start", "end", 
                                   "strand", "col", "fill"),
                     show_col_types = FALSE)
  region <- region %>% 
    # color black if it is "lonely" gene
    mutate(col = ifelse(name %in% lonely_g$id, "black", col), 
           fill = ifelse(name %in% lonely_g$id, "black", fill)) %>%
    # color blue if tandem repeat
    mutate(col = ifelse(name %in% tandem_rep$id, "blue", col), 
           fill = ifelse(name %in% tandem_rep$id, "blue", fill)) %>%
    # color green if match with other "lonely" gene
    mutate(col = ifelse(name %in% matched_genes$id, "green", col), 
           fill = ifelse(name %in% matched_genes$id, "green", fill)) %>%
    # color purple if match with gene between same chromosomes
    mutate(col = ifelse(name %in% match_in_chrom$value, "purple", col), 
           fill = ifelse(name %in% match_in_chrom$value, "purple", fill)) %>%
    # color yellow if match within multiplicon but not between lonely genes
    mutate(col = ifelse(name %in% match_g_in_mult$value, "yellow", col), 
           fill = ifelse(name %in% match_g_in_mult$value, "yellow", fill)) %>%
    # color pink if match with gene between different chromosomes
    mutate(col = ifelse(name %in% match_scattered$lonely_gene, "pink", col), 
           fill = ifelse(name %in% match_scattered$lonely_gene, "pink", fill))
}

# Color all genes according to type
r1 <- color_genes(region1, lonely_g, tandem_rep, matched_genes, match_in_chrom)
write_tsv(r1, region1, col_names = FALSE)
r2 <- color_genes(region2, lonely_g, tandem_rep, matched_genes, match_in_chrom)
write_tsv(r2, region2, col_names = FALSE)
#r2_rev <- color_genes(region2_rev, lonely_g, tandem_rep, matched_genes, match_in_chrom)
#write_tsv(r2_rev, region2_rev, col_names = FALSE)