library(tidyverse)
library(UpSetR)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
working_dir <- args[1]

# Read in the data
level_all <- read.table(paste0(working_dir, "/results/i-ADHoRe-run/table_level_all_upd.tsv"),
                        header = TRUE) %>% as_tibble()

missing_coll <- read_tsv(paste0(working_dir, "/results/i-ADHoRe-run/collinear_missing.txt"), 
                        col_names = c("group_id", "subgenome", "coll"))

# Summarize the data
collinear_missing <- missing_coll %>% 
  filter(coll == TRUE) %>%
  group_by(group_id) %>%
  summarise(subgenome = paste(subgenome, collapse = ",")) %>%
  ungroup()

## Add the sub-genomes where a DNA deletion was detected
df <- left_join(level_all, collinear_missing, by = "group_id") %>%
  rename(del_subgenomes = subgenome) %>%
  mutate(num_del = str_count(del_subgenomes, ",") + 1,
         num_anchored = num_aps + num_unan) %>% # unannotated genes are counted as collinear genes
  mutate(num_del = replace_na(num_del, 0)) %>%
  rowwise() %>%
  mutate(num_noncoll = 
           4 - min(num_anchored + num_del + num_psgs + num_trans, 4)) %>% # these are the ones if non other are detected
  ungroup() %>%
  # If there is nothing homologous to a gene 
  # (no translocated, collinear homologous gene, pseudogene, deleted gene)
  # then we can define it as a singleton
  mutate(num_singleton = ifelse(num_anchored == 1 & num_del == 0 & 
                                  num_trans == 0 & num_psgs == 0,
                                1, 0)) %>%
  mutate(num_anchored = ifelse(num_singleton == 1, 0, num_anchored),
         num_noncoll = ifelse(num_singleton == 1, 0, num_noncoll)) %>%
  select(group_id, num_aps, num_unan, num_anchored, 
         num_trans, num_psgs, num_del, num_noncoll, num_singleton,
         ap_subgenomes,unan_subgenomes, psgs_subgenomes, del_subgenomes, 
         everything())

counts <- df %>%
  rename("Collinear genes" = num_anchored, "Translocated genes" = num_trans,
         "Pseudogenes" = num_psgs, "Deleted genes" = num_del,
         "No collinearity" = num_noncoll,
         "Singleton genes" = num_singleton) %>%
  filter(redundant == "False") %>% # remove redundant groups
  summarise(across(c(`Collinear genes`:`Singleton genes`), sum)) %>%
  pivot_longer(cols = everything(), names_to = "categories")

# Write the summarized data
write_tsv(df, paste0(working_dir, "/results/results_table.tsv"))

# Plot the results
presence <- df %>% 
  mutate(across(num_aps:num_singleton, ~ ifelse(. != 0, 1, 0))) %>% # change the counts to presence/absence
  filter(redundant == "False") %>% # remove redundant groups
  data.frame() %>%
  rename("With collinear genes" = num_anchored, "With translocated genes" = num_trans,
         "With pseudogenes" = num_psgs, "With deleted genes" = num_del,
         "With no collinearity" = num_noncoll,
         "Singleton gene" = num_singleton) 

## UpSet plot
plot <- upset(presence, sets = c("With collinear genes",
                                 "With translocated genes",
                                 "With pseudogenes",
                                 "With deleted genes", 
                                 "With no collinearity",
                                 "Singleton gene"),
              # set.metadata = list(data = counts %>%
              #                       rename("Absolute_Number" = value), 
              #                     plots = list(list(type = "hist", 
              #                                       column = "Absolute_Number",
              #                                       assign = 20,
              #                                       colors = "grey"))),
              order.by = "freq", 
              mainbar.y.label = "Number of collinear gene groups",
              sets.x.label = "Number of collinear gene groups",
              set_size.show = TRUE,
              set_size.numbers_size = 8,
              set_size.scale_max = 45000,
              text.scale = 1.5)

svg(paste0(working_dir, "/results/upsetplot.svg"), width = 14, height = 7)
plot
dev.off()

## Bar/dot plot
### In the UpSet plot, the set size side plot only shows the presence/absence of a category in a group
### it is probably more interesting to have the total number for each category
### e.g. if there are more than one pseudogenes in a group, count them also more than once
dotplot <- ggplot(counts, 
       aes(y = reorder(categories, as.numeric(value), decreasing = TRUE), x = "",
           size = value)) +
  geom_point(stat = "identity") + 
  geom_text(aes(label = value),
                      size = 6, hjust = -0.25) +
  ylab("") + xlab("") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank(),
        legend.position = "none")

svg(paste0(working_dir, "/results/dotplot.svg"), width = 5, height = 5)
dotplot
dev.off()

# barplot <- ggplot(counts %>% filter(!categories %in% c("num_aps", "num_unan")), 
#                   aes(y = reorder(categories, value, decreasing = TRUE), x = value)) +
#   geom_bar(stat = "identity") + 
#   geom_text(aes(label = value),
#             size = 6, hjust = -0.25) +
#   ylab("") + xlab("") +
#   ggthemes::theme_base() +
#   theme(text = element_text(size = 20))
# 
# svg(paste0(working_dir, "/results/barplot.svg"), width = 14, height = 7)
# barplot
# dev.off()

## UpSet plot per level
upsetplot_level <- function(df, level){
  presence <- df %>% filter(num_anchored == level,
                            redundant == "False") %>%
    mutate(across(num_aps:num_noncoll, ~ ifelse(. != 0, 1, 0))) %>%
    data.frame() %>%
    rename("With collinear genes" = num_anchored, "With translocated genes" = num_trans,
         "With pseudogenes" = num_psgs, "With deleted genes" = num_del,
         "With no collinearity" = num_noncoll, "Singleton gene" = num_singleton)
  
  counts <- df %>% filter(num_anchored == level) %>%
    rename("With collinear genes" = num_anchored, "With translocated genes" = num_trans,
           "With pseudogenes" = num_psgs, "With deleted genes" = num_del,
           "With no collinearity" = num_noncoll,
           "Singleton gene" = num_singleton) %>%
    summarise(across(c(`With collinear genes`:`Singleton gene`), sum)) %>%
    pivot_longer(cols = everything(), names_to = "categories")
  
  plot <- upset(presence, sets = c("With collinear genes",
                                   "With translocated genes",
                                   "With pseudogenes",
                                   "With deleted genes",
                                   "With no collinearity",
                                   "Singleton gene"),
              # set.metadata = list(data = counts %>%
              #                       rename("Absolute_Number" = value),
              #                     plots = list(list(type = "hist",
              #                                       column = "Absolute_Number",
              #                                       assign = 20,
              #                                       colors = "grey"))),
              order.by = "freq", mainbar.y.label = "Number of collinear gene groups",
              set_size.show = TRUE,
              set_size.numbers_size = 8,
              text.scale = 1.5)
}

svg(paste0(working_dir, "/results/upsetplot_level_4.svg"), width = 14, height = 7)
plot <- upsetplot_level(df, 4)
plot
dev.off()

svg(paste0(working_dir, "/results/upsetplot_level_3.svg"), width = 14, height = 7)
plot <- upsetplot_level(df, 3)
plot
dev.off()

svg(paste0(working_dir, "/results/upsetplot_level_2.svg"), width = 14, height = 7)
plot <- upsetplot_level(df, 2)
plot
dev.off()

svg(paste0(working_dir, "/results/upsetplot_level_1.svg"), width = 14, height = 7)
plot <- upsetplot_level(df, 1)
plot
dev.off()

## Dot plot per level
dotplot_level <- function(df, level){
  counts <- df %>%
    filter(num_anchored == level,
           redundant == "False") %>%
    rename("Collinear genes" = num_anchored, "Translocated genes" = num_trans,
           "Pseudogenes" = num_psgs, "Deleted genes" = num_del,
           "No collinearity" = num_noncoll,
           "Singleton genes" = num_singleton) %>%
    summarise(across(c(`Collinear genes`:`Singleton genes`), sum)) %>%
    pivot_longer(cols = everything(), names_to = "categories")
  
  dotplot <- ggplot(counts, 
                    aes(y = reorder(categories, as.numeric(value), 
                                    decreasing = TRUE), x = "",
                        size = value)) +
    geom_point(stat = "identity") + 
    geom_text(aes(label = value),
              size = 6, hjust = -0.25) +
    ylab("") + xlab("") +
    theme_minimal() +
    theme(text = element_text(size = 20),
          panel.grid = element_blank(),
          legend.position = "none")
}

# barplot_level <- function(df, level){
#   counts <- df %>%
#     filter(num_anchored == level) %>%
#     rename("With collinear genes" = num_anchored, "With translocated genes" = num_trans,
#             "With pseudogenes" = num_psgs, "With deleted genes" = num_del,
#             "With no collinearity" = num_noncoll,
#            "Singleton gene" = num_singleton) %>%
#     summarise(across(c(`With collinear genes`:`Singleton gene`), sum)) %>%
#     pivot_longer(cols = everything(), names_to = "categories")
# 
#     barplot <- ggplot(counts %>% filter(!categories %in% c("num_aps", "num_unan")), 
#         aes(y = reorder(categories, value, decreasing = TRUE), x = value)) +
#     geom_bar(stat = "identity") + 
#     geom_text(aes(label = value), size = 6, hjust = -0.25) +
#     ylab("") + xlab("") +
#     ggthemes::theme_base() +
#     theme(text = element_text(size = 20))
# }

svg(paste0(working_dir, "/results/dotplot_level_4.svg"), width = 5, height = 5)
plot <- dotplot_level(df, 4)
plot
dev.off()

svg(paste0(working_dir, "/results/dotplot_level_3.svg"), width = 5, height = 5)
plot <- dotplot_level(df, 3)
plot
dev.off()

svg(paste0(working_dir, "/results/dotplot_level_2.svg"), width = 5, height = 5)
plot <- dotplot_level(df, 2)
plot
dev.off()

svg(paste0(working_dir, "/results/dotplot_level_1.svg"), width = 5, height = 5)
plot <- dotplot_level(df, 1)
plot
dev.off()

## Collinearity gene retention level bar plot
level_gene <- df %>% select(num_anchored, genes) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(num_anchored = ifelse(num_anchored == 0, 1, num_anchored))
level_count <- level_gene %>% group_by(num_anchored) %>%
  count()

plot <- ggplot(level_count, aes(x = num_anchored, y = n)) +
  geom_bar(stat = "identity") +
  xlab("Gene retention level on sub-genomes") +
  ylab("Number of genes") +
  geom_text(aes(label = n),
            size = 6, vjust = -0.5) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))

svg(paste0(working_dir, "/results/coll_level_barplot.svg"), width = 6, height = 7)
plot
dev.off()