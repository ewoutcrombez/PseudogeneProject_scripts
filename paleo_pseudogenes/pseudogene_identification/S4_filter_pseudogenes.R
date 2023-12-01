library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
species = args[1]
plaza_name = args[2]

# load in data
psg <- read_tsv(str_glue("../PseudoPipe_results/{species}_pgenes.txt"),
                    col_names =
                      c("chr", "start", "end", "strand", "query", "frac", "ins",
                        "del", "shift", "stop", "expect",
                        "ident", "polya", "type", "psg_name")) %>%
                        mutate(length = end - start + 1)
filt_psg <- psg %>%
	filter(expect <= 1e-5) %>% 
	filter(ident >= 0.2) %>% 
	filter(frac >= 0.05) %>%
  filter(length >= 90)

gff <- read_tsv(str_glue("/scratch/recent_wgds/data/{species}/gff/annotation.selected_transcript.all_features.{plaza_name}.gff3"),
         col_names = FALSE) %>%
  filter(X3 == "mRNA") %>%
  mutate(X9 = str_extract(X9, "ID=(.*?);") %>%
  str_replace("ID=", "") %>%
  str_replace(";", "")) %>%
  select(X1, X9)

unproc_psg <- left_join(filt_psg, gff, by = c("query" = "X9")) %>%
  # This should not be done for Amborella which consists of scaffolds
  filter(!str_detect(X1, "scaff")) %>%
  filter(!str_detect(X1, "scaf")) %>%
  filter(!str_detect(X1, "random")) %>%
  filter(!str_detect(X1, "super")) %>%
  filter(chr != "ChrM" & !startsWith(query, "ATM"))

unproc_psg <- unproc_psg %>%
  select(-X1)

PSSD <- filt_psg %>%
	filter(type == "PSSD") %>%
	nrow() %>%
	as.character()
FRAG <- filt_psg %>%
	filter(type == "FRAG") %>%
	nrow() %>%
	as.character()

DUP <- filt_psg %>%
	filter(type == "DUP") %>%
	nrow() %>%
	as.character()

# print number of pseudogenes in each category (PSSD, FRAG, DUP)
print(str_glue("number of PSSD pseudogenes: {PSSD}"))
print(str_glue("number of FRAG pseudogenes: {FRAG}"))
print(str_glue("number of DUP pseudogenes: {DUP}"))
# Total number of pseudogenes
print(str_glue("Total number of identified pseudogenes: {nrow(filt_psg)}"))

# get unprocessed pseudogenes
unproc_psg <- unproc_psg %>%
  filter(type != "PSSD")

## Save pseudogenes
write_tsv(unproc_psg, str_glue("{species}_filtered_low_pgenes_unprocessed.txt"),
          col_names = FALSE)

write_tsv(filt_psg, str_glue("{species}_filtered_low_pgenes_all.txt"),
          col_names = FALSE)