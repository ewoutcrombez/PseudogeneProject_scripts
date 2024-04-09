library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
tsv_file = args[1]
gff_file = args[2]
gff_out_file = args[3]
tsv_out_file = args[4]

# Read in files
## Here we assume that hits are in same order in both files
## (which should always be the case from exonerate output)
tsv <- read_tsv(tsv_file, col_names =
    c("query_id", "target_id", "q_strand", "t_strand",
    "q_length", "t_aln_start", "t_aln_end", "t_aln_length",
    "q_aln_start", "q_aln_end", "q_aln_length", "perc_id",
    "raw_score", "total_equivalenced", "identical",
    "similarity", "mismatches"), show_col_types = FALSE) %>%
    mutate(hit_num = paste0("hit_", row_number())) # add hit number
head(tsv)
write_tsv(tsv, tsv_file, col_names = FALSE)
gff <- read_tsv(gff_file, col_names =
    c("chrom", "source", "type", "start", "end", "score",
    "strand", "phase", "attributes"),
    show_col_types = FALSE) %>%
    filter(type == "gene") %>% # only keep gene rows
    mutate(hit_num = paste0("hit_", row_number())) # add hit number
head(gff)
write_tsv(gff, gff_out_file, col_names = FALSE)

if (nrow(tsv) > 0 & nrow(gff) > 0){
  # Pick the hit with the longest length if multiple hits match to same
  # region of the gene
  tsv_filt <- tsv %>% group_by(query_id) %>%
      arrange(q_aln_start) %>%
      mutate(index = cumsum(cummax(lag(q_aln_end, default = first(q_aln_end))) < q_aln_start)) %>% # create index
      group_by(index, query_id) %>%
      filter(q_aln_length == max(q_aln_length)) %>%
      distinct() %>%
      ungroup()
  # Next pick the one (if still multiple because same max length), pick the one
  # with the highest percent identity
  tsv_filt <- tsv_filt %>%
    group_by(index, query_id) %>%
    filter(perc_id == max(perc_id)) %>%
    ungroup() %>%
    select(-index)

  # Select the GFF lines that correspond to filtered hits
  gff_filt <- gff %>%
    filter(gff$hit_num %in% tsv_filt$hit_num)

  # Write out filtered files
  tsv_filt %>%
    write_tsv(tsv_out_file, col_names = FALSE)
  gff_filt %>%
    write_tsv(gff_out_file, col_names = FALSE)
}
