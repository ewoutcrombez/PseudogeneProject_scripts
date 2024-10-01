# In this script I will extract the regions that are part of a collinear region 
# and mask the genic regions in this regions
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
working_dir = args[1]
multiplicon = args[2]
haplotype_1 = args[3]
chrom_1 = args[4]
gff_1 = args[5]
region_1 = args[6]
haplotype_2 = args[7]
chrom_2 = args[8]
gff_2 = args[9]
region_2 = args[10]

setwd(working_dir)

# read in data of first region
gff_t1 <- read_tsv(str_glue("results/i-ADHoRe-run/processing/{multiplicon}/tmp_gff_{haplotype_1}"), 
                   col_names = 
                     c("chrom", "source", "type", "start", "end", "x",
                                 "or", "y", "gene"),
                   show_col_types = FALSE)
reg_t1 <- read_tsv(str_glue("results/i-ADHoRe-run/processing/{multiplicon}/region_{haplotype_1}_{chrom_1}.txt"), 
                        col_names = FALSE, show_col_types = FALSE)
# read in data of second region
gff_t2 <- read_tsv(str_glue("results/i-ADHoRe-run/processing/{multiplicon}/tmp_gff_{haplotype_2}"), 
                   col_names = 
                     c("chrom", "source", "type", "start", "end", "x",
                       "or", "y", "gene"),
                   show_col_types = FALSE)
reg_t2 <- read_tsv(str_glue("results/i-ADHoRe-run/processing/{multiplicon}/region_{haplotype_2}_{chrom_2}.txt"), 
                   col_names = FALSE, show_col_types = FALSE)

# combine gff and region data
coords_r1 <- inner_join(reg_t1, gff_t1, by = c("X1" = "gene"))
coords_r2 <- inner_join(reg_t2, gff_t2, by = c("X1" = "gene"))

# get bed files for these regions
bed_r1 <- tibble(chrom = coords_r1$chrom[1],
              start = max(min(c(coords_r1$start, coords_r1$end)) - 1001, 0), # cannot go below 0
              end =max(c(coords_r1$start, coords_r1$end)) + 1000,
              name = str_glue("{haplotype_1}_{chrom_1}"))
write_tsv(bed_r1, str_glue("results/i-ADHoRe-run/processing/{multiplicon}/region_{haplotype_1}_{chrom_1}.bed"),
          col_names = FALSE)
bed_r2 <- tibble(chrom = coords_r2$chrom[1],
                 start = max(min(c(coords_r2$start, coords_r2$end)) - 1001, 0), # cannot go below 
                 end = max(c(coords_r2$start, coords_r2$end)) + 1000,
                 name = str_glue("{haplotype_2}_{chrom_2}"))
write_tsv(bed_r2, str_glue("results/i-ADHoRe-run/processing/{multiplicon}/region_{haplotype_2}_{chrom_2}.bed"),
          col_names = FALSE)

# change global coordinates of genes to the local region coordinates
coords_r1 %>% mutate(start = start - bed_r1$start - 1,
                     end = end - bed_r1$start,
                     name = str_glue("{haplotype_1}_{chrom_1}")) %>%
  select(name, start, end, X1) %>%
  write_tsv(str_glue("results/i-ADHoRe-run/processing/{multiplicon}/genes_{haplotype_1}_{chrom_1}.bed"),
            col_names = FALSE)
coords_r2 %>% mutate(start = start - bed_r2$start - 1,
                     end = end - bed_r2$start,
                     name = str_glue("{haplotype_2}_{chrom_2}")) %>%
  select(name, start, end, X1) %>%
  write_tsv(str_glue("results/i-ADHoRe-run/processing/{multiplicon}/genes_{haplotype_2}_{chrom_2}.bed"),
            col_names = FALSE)

# save region start coordinates
write_tsv(bed_r1 %>% select(start), 
          str_glue("results/i-ADHoRe-run/processing/{multiplicon}/region_{haplotype_1}_{chrom_1}_start.txt"),
          col_names = FALSE)
write_tsv(bed_r2 %>% select(start), 
          str_glue("results/i-ADHoRe-run/processing/{multiplicon}/region_{haplotype_2}_{chrom_2}_start.txt"),
          col_names = FALSE)
