library(tidyverse)
setwd("/scratch/recent_wgds/paleo_polyploids")

args <- commandArgs(trailingOnly = TRUE)

species <- args[1]

blastp <- 
  read_tsv(stringr::str_glue("allvsall_blastp/{species}_allvsallblastp_filtered_wlen.tsv"),
           col_names = FALSE)

# remove self-hits
blastp <- blastp %>% filter(X1 != X2)

# select first two columns to create blast_table as input for i-ADHoRe
# and only keep distinct pairs
blastp <- blastp %>% select(X1, X2) %>%
  distinct()

psg <- read_tsv(stringr::str_glue("/scratch/recent_wgds/paleo_polyploids/collinearity/filtered_no_intron/{species}_filtered_no_intronic.txt"),
                col_names = FALSE) %>%
  select(X5, X15) %>%
  rename(X1 = X5, X2 = X15)

table <- rbind(blastp, psg)


write_tsv(table, stringr::str_glue("collinearity/i-adhore/no_introns2/{species}/table.txt"),
          col_names = FALSE)
