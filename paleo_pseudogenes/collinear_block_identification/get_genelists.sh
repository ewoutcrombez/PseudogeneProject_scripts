#!/bin/bash
# In this script i-ADHoRe is run with gene pairs and pseudogene-parent gene pairs as input

## Input in command line
species=$1 # e.g. Arabidopsis_thaliana
plaza_name=$2 # species abbreviation used in PLAZA, e.g. ath

## Directories to use
wd="/scratch/recent_wgds/paleo_polyploids/collinearity"
gff="/scratch/recent_wgds/data/$species/gff/annotation.selected_transcript.all_features.${plaza_name}.gff3"
# after low filtering
#pseudogenes="/scratch/recent_wgds/paleo_polyploids/filtering/${species}_filtered_low_pgenes_unprocessed.txt"
# after filtering out pseudogenes that reside in introns
pseudogenes="/scratch/recent_wgds/paleo_polyploids/collinearity/filtered_no_intron/${species}_filtered_no_intronic.txt"

analysis="i-adhore/no_introns2/${species}"

## Create genome lists
mkdir $analysis/genome_lists
cd $analysis/genome_lists
### Process gff file
awk -F "\t" '$3 == "mRNA"' $gff | sed -r 's/ID=([^;]*);.*/\1/' | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' > genes.txt
### Process pseudogene txt file
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$15}' $pseudogenes > pgenes.txt

### Combine processed files and create a file per chromosome/scaffold/...
cat genes.txt pgenes.txt | sort -k1,1 -k2,2n | awk '{print $1"\t"$5$4}' | awk -F "\t" '{print $2>$1}'
rm genes.txt pgenes.txt
### If chromosome level assembly remove scaffolds, except for Amborella:
if [[ $species != "Amborella_trichopoda" ]]
then
    rm *scaf*
    rm *random*
    rm *super*
    rm chrUn
    rm ChrM
    rm ChrC
fi