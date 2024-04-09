#!/bin/bash

species=$1
plaza_name=$2

module load bedtools

echo "working on ${species}"

# Obtain mRNA GFF file
awk -F "\t" '$3 == "mRNA"' /scratch/recent_wgds/data/$species/gff/annotation.selected_transcript.all_features.${plaza_name}.gff3 > /scratch/recent_wgds/data/$species/gff/annotation.selected_transcript.mRNA.${plaza_name}.gff3

# Check which pseudogenes overlap with an mRNA (so not only exons)
# The ones that overlap with more than 30 bp are considered to be intron pseudogenes
bedtools intersect -a /scratch/recent_wgds/paleo_polyploids/filtering/${species}_filtered_full_unprocessed.txt -b /scratch/recent_wgds/data/${species}/gff/annotation.selected_transcript.mRNA.${plaza_name}.gff3 -wo | awk '$26 > 30' | cut -f15 | sort | uniq > intron_pseudogenes_${plaza_name}.txt

# Create regular expression so that we can remove intron pseudogenes
sed -r 's/(.*)/\1\t\[0-9\]\+/' intron_pseudogenes_${plaza_name}.txt > tmp && mv tmp intron_pseudogenes_${plaza_name}.txt

# Remove intron pseudogenes from list
grep -E -v -f intron_pseudogenes_${plaza_name}.txt /scratch/recent_wgds/paleo_polyploids/filtering/${species}_filtered_full_unprocessed.txt > ${species}_filtered_no_intronic.txt

# Remove PSSD pseudogenes from list
awk -F "\t" '$14 != "PSSD"' ${species}_filtered_no_intronic.txt > tmp && mv tmp ${species}_filtered_no_intronic.txt