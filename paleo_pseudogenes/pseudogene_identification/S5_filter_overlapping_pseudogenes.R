library(tidyverse)
library(GenomicRanges)

# load in data
pseudogenes <- read_tsv(
  "pgenes_with_overlapping.txt",
  col_names = c("chrom", "start", "end", "orientation",
                "parent_gene", "fraction_parent_gene", "insertion",
                "deletion", "shift", "stop", "e-val", "ident",
                "polya", "type", "pseudogene_id", "length_pseudogene", 
                "species"))

# Create GRanges object
gr <- with(pseudogenes, 
           GRanges(seqnames = chrom, ranges = IRanges(start = start,
                                                      end = end),
                   strand = orientation, mcols = data.frame(name = pseudogene_id,
                                                                type = type,
                                                                eval = `e-val`,
                                                                ident = ident,
                                                                fraction_parent_gene =
                                                                  fraction_parent_gene
                       )))

# Split GRanges by species
grl <- split(gr, f = pseudogenes$species)

# Get overlaps of pseudogenes (> 30bp)
overlaps <- lapply(grl, function(x) findOverlaps(x, type = "any",
                                                 minoverlap = 30,
                                                 ignore.strand = TRUE))
# Remove self-hits
overlaps <- lapply(overlaps, function(x) x[queryHits(x) != subjectHits(x)])

# Get the pseudogenes that overlap
grl_overlaps_query <- mapply(function(x, y) x[queryHits(y)],
                             x = grl, y = overlaps)
grl_overlaps_subject <- mapply(function(x, y) x[subjectHits(y)],
                               x = grl, y = overlaps)
overlaps_query <- lapply(grl_overlaps_query, as.data.frame)
overlaps_subject <- lapply(grl_overlaps_subject, as.data.frame)

# Combine overlapping pseudogenes into one dataframe
combined_df <- bind_rows(Map(cbind, species = names(overlaps_query), 
                             overlaps_query))
combined_df2 <- bind_rows(Map(cbind, species = names(overlaps_subject), 
                              overlaps_subject))
combined <- cbind(combined_df, combined_df2)
names(combined) <- c("species1", "chrom1", "start1", "end1", "width1",
                     "strand1", "pseudogene_id1", "type1", "eval1",
                     "ident1", "fraction1", "species2", "chrom2", "start2", 
                     "end2", "width2",
                     "strand2", "pseudogene_id2", "type2", "eval2",
                     "ident2", "fraction2")

# Check which hit matches best and select the other to be removed
combined <- combined %>%
  mutate(removed = ifelse(eval1 < eval2,
                          pseudogene_id2,
                          ifelse(eval1 > eval2,
                                 pseudogene_id1, ifelse(
                                   fraction1 > fraction2,
                                   pseudogene_id1,
                                   ifelse(fraction2 > fraction1,
                                          pseudogene_id2, 
                                          ifelse(ident1 > ident2,
                                                 pseudogene_id1,
                                                 pseudogene_id2))
                                 ))))

# Pseudogenes to be removed
remove <- combined %>% select(species1, removed) %>% distinct() %>%
  mutate(spec_pseudogene = paste(species1, removed, sep = "|"))

# Number per species that will be removed
remove %>%
  group_by(species1) %>%
  count()

# Filter out pseudogenes that should be removed because they overlap with a pseudogene
# that is a better hit
pseudogenes <- 
  pseudogenes %>% mutate(spec_pseudogene = paste(species, pseudogene_id, sep = "|"))

filtered <- 
  pseudogenes %>% filter(!spec_pseudogene %in% remove$spec_pseudogene) %>%
  select(-spec_pseudogene)

# Write out filtered pseudogenes
write_tsv(filtered, 
          "pgenes.txt",
          col_names = FALSE)
# Group by species and write each group to a separate file
filtered %>%
  group_split(species) %>%
  walk2(
    ., paste0(unique(filtered$species), "_filtered_nonoverlap.txt"),
    ~ write_tsv(.x, .y, col_names = FALSE)
  )

species <- filtered %>% select(species) %>% distinct()

for (spec in species$species){
  spec_df_filt <- read_tsv(glue::glue("{spec}_filtered_low_pgenes_all.txt"),
                col_names = c("chrom", "start", "end", "orientation",
                              "parent_gene", "fraction_parent_gene", "insertion",
                              "deletion", "shift", "stop", "e-val", "ident",
                              "polya", "type", "pseudogene_id", "length_pseudogene", 
                              "species"))
  spec_df_filt_unproc <- read_tsv(glue::glue("{spec}_filtered_low_pgenes_unprocessed.txt"),
                                col_names = c("chrom", "start", "end", "orientation",
                                              "parent_gene", "fraction_parent_gene", "insertion",
                                              "deletion", "shift", "stop", "e-val", "ident",
                                              "polya", "type", "pseudogene_id", "length_pseudogene", 
                                              "species"))
  spec_df_noverlap <- read_tsv(glue::glue("{spec}_filtered_nonoverlap.txt"),
                                col_names = c("chrom", "start", "end", "orientation",
                                              "parent_gene", "fraction_parent_gene", "insertion",
                                              "deletion", "shift", "stop", "e-val", "ident",
                                              "polya", "type", "pseudogene_id", "length_pseudogene", 
                                              "species"))
  spec_df <- spec_df_filt %>% filter(pseudogene_id %in% spec_df_noverlap$pseudogene_id)
  spec_df_unproc <- spec_df_filt_unproc %>% filter(pseudogene_id %in% spec_df_noverlap$pseudogene_id)
  write_tsv(spec_df,
            glue::glue("{spec}_filtered_full_all.txt"),
                 col_names = FALSE)
  write_tsv(spec_df_unproc,
            glue::glue("{spec}_filtered_full_unprocessed.txt"),
                 col_names = FALSE)
}

# combine all "{spec}_filtered_full_all.txt" files in one using cat > pgenes_full_filtered.txt