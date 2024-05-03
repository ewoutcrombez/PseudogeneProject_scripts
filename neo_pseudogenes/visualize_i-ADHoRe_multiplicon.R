# In this script I make alignment plot of a multiplicon
# using the genoPlotR library
## Load in library
library(genoPlotR)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
multiplicon = args[1]
species1 = args[2]
chromosome1 = args[3]
species2 = args[4]
chromosome2 = args[5]
output = args[6]
with_pseudogenes = args[7]

## Get table that annotates pseudogenes
if (with_pseudogenes == "YES") {
  ann_table <- read_tsv(sprintf("%s/%s/result_table_%s_MACSE.tsv",
    output, multiplicon, multiplicon), col_names = TRUE,
    show_col_types = FALSE)
}

## Make lists of all genes in multiplicon
if (with_pseudogenes == "NO") {
  df1 <- read_tsv(sprintf("%s/%s/region_%s_%s.txt", output, multiplicon,
                          species1, chromosome1), col_names = FALSE,
                          show_col_types = FALSE)
  df2 <- read_tsv(sprintf("%s/%s/region_%s_%s.txt", output, multiplicon,
                          species2, chromosome2), col_names = FALSE,
                          show_col_types = FALSE)
  #df2_reversed <- read.table(sprintf("%s/%s/region_%s_%s_reversed.txt",
  #  output, multiplicon, species2, chromosome2), header = FALSE)
  colnames(df1) <- c("name", "start", "end", "strand", "col", "fill")
  colnames(df2) <- c("name", "start", "end", "strand", "col", "fill")
  #colnames(df2_reversed) <- c("name", "start", "end", "strand", "col", "fill")
  fh <- read_tsv(sprintf("%s/%s/multiplicon_%s_pairs.txt",
    output, multiplicon, multiplicon), col_names = FALSE,
    show_col_types = FALSE)
} else {
  df1 <- read_tsv(sprintf("%s/%s/%s.tab.txt", output, multiplicon,
                          species1), col_names = TRUE,
                          show_col_types = FALSE)
  df2 <- read_tsv(sprintf("%s/%s/%s.tab.txt", output, multiplicon,
                          species2), col_names = TRUE,
                          show_col_types = FALSE)
  fh <- read_tsv(sprintf("%s/%s/multiplicon_pairs_with_pseudogenes.txt",
    output, multiplicon), col_names = FALSE,
    show_col_types = FALSE)
}

## Make coordinates for comparisons
fh <- fh %>%
  mutate(
    X4 = if_else(str_starts(X1, species1), X1, X2),
    X5 = if_else(str_starts(X2, species2), X2, X1)
  ) %>%
  select(-X1, -X2) %>%
  rename(X1 = X4, X2 = X5)

# unannot_genes <- ann_table %>%
#                     filter(annotation == "intact_unannotated_gene" |
#                            annotation == "likely_unannotated_gene")

# dubious_psg <- ann_table %>%
#                     filter(annotation == "dubious_pseudogene")

# if (with_pseudogenes == "YES") {
#   fh <- fh %>%
#     mutate(
#       X3 = ifelse(X1 %in% unannot_genes$hit_synonym | 
#                   X2 %in% unannot_genes$hit_synonym, "unann_gene", X3)) %>%
#     mutate(
#       X3 = ifelse(X1 %in% dubious_psg$hit_synonym |
#                   X2 %in% dubious_psg$hit_synonym, "dubious_pseudogene", X3))
# }

## Make unique list of all genes
genes <- union(df1$name, df2$name)
genes <- unique(genes)

## Create genoPlotR dna_segs
dna_seg1 <- dna_seg(df1)
dna_seg2 <- dna_seg(df2)
dna_segs <- list(dna_seg1, dna_seg2)

#dna_seg2_reversed <- dna_seg(df2_reversed)
#dna_segs_reversed <- list(dna_seg1, dna_seg2_reversed)

## Obtain coordinates from dna_seg's
start1 <- sapply(fh$X1, function(x) dna_seg1[dna_seg1$name == x, ]$start)
end1 <- sapply(fh$X1, function(x) dna_seg1[dna_seg1$name == x, ]$end)
start2 <- sapply(fh$X2, function(x) dna_seg2[dna_seg2$name == x, ]$start)
end2 <- sapply(fh$X2, function(x) dna_seg2[dna_seg2$name == x, ]$end)


# col <- ifelse(fh$X3 == "anchorpair", "grey",
#               ifelse(fh$X3 == "matched_lonely_pairs", "grey",
#                      ifelse(fh$X3 == "tandem_pairs", "grey",
#                             ifelse(fh$X3 == "multiplicon_pairs", "grey",
#                               ifelse(fh$X3 ==
#                                 "gene_pseudogene_pairs",
#                                      "blue", 
#                                 ifelse(fh$X3 ==
#                                          "unann_gene", "grey", 
#                                          ifelse(fh$X3 ==
#                                             "dubious_pseudogene", "yellow", "black")))))))

col <- ifelse(fh$X3 == "anchorpair", "grey",
  ifelse(fh$X3 == "matched_lonely_pairs", "grey",
        ifelse(fh$X3 == "tandem_pairs", "grey",
                ifelse(fh$X3 == "multiplicon_pairs", "grey",
                  ifelse(fh$X3 ==
                    "gene_pseudogene_pairs",
                        "blue", "black")))))

#start2_reversed <- sapply(
#  fh$X2, function(x) dna_seg2_reversed[dna_seg2_reversed$name == x, ]$start)
#end2_reversed <- sapply(
#  fh$X2, function(x) dna_seg2_reversed[dna_seg2_reversed$name == x, ]$end)

## Add annotation (if short enough multiplicon region)
annot1 <- annotation(x1 = dna_seg1$start, text = dna_seg1$name, rot = 60)
annot2 <- annotation(x1 = dna_seg2$start, text = dna_seg2$name, rot = 60)
annots <- list(annot1, annot2)
## Keeping only one every 4 annots
annots <- lapply(annots, function(x) x[(1:nrow(x))%%3 == 0,])

## Create genoPlotR comparisons
comparison1 <- as.comparison(data.frame(
  start1 = start1, end1 = end1, start2 = start2, end2 = end2, col = col))
comparisons <- list(comparison1)

#comparison1_reversed <- as.comparison(data.frame(
#  start1 = start1, end1 = end1,
#  start2 = start2_reversed, end2 = end2_reversed, col = col))
#comparisons_reversed <- list(comparison1_reversed)

## Adding names to plot
names <- c(sprintf("%s|%s", species1, chromosome1),
  sprintf("%s|%s", species2, chromosome2))
names(dna_segs) <- names
#names(dna_segs_reversed) <- names

# Make figure
jpeg(sprintf("%s/%s/multiplicon_%s_%s_%s.jpeg", output, multiplicon,
             multiplicon, species1, species2), width = 5000, height = 800)
plot_gene_map(dna_segs = dna_segs, comparisons = comparisons, scale = FALSE,
              dna_seg_scale = FALSE, dna_seg_label_cex = 1, annotation_cex = 1,
              seg_plot_height = 1)
dev.off()

#jpeg(sprintf("%s/%s/multiplicon_%s_%s_%s_reversed.jpeg", output, multiplicon,
#             multiplicon, species1, species2), width = 5000, height = 800)
#plot_gene_map(dna_segs = dna_segs_reversed, comparisons = comparisons_reversed,
#              scale = FALSE, dna_seg_scale = FALSE,
#              dna_seg_label_cex = 1, annotation_cex = 1)
#dev.off()
