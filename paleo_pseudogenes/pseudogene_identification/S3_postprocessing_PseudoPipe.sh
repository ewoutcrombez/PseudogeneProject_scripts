#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 15G # memory pool for all cores
#SBATCH -o slurm.%j_Vitis_post.out # STDOUT
#SBATCH -e slurm.%j_Vitis_post.err # STDERR
# In this script I do some post-processing of the PseudoPipe results file

## Required modules
module load python
module load bedtools

## Input in command line
species=$1 # e.g. Arabidopsis_thaliana
plaza_name=$2 # species abbreviation used in PLAZA, e.g. ath

echo "Started post-processing for species $species!"
## Create result directory
mkdir /scratch/recent_wgds/results/pseudogenes/Pseudopipe/${species}/
cd /scratch/recent_wgds/results/pseudogenes/Pseudopipe/${species}/

## Collect raw results of each chromosome from PseudoPipe
mkdir 0_raw_result
cp /scratch/recent_wgds/analysis/pgenes/ppipe_output/${species}PerChr/*/pgenes/*_pgenes.* 0_raw_result
find . -type f -empty -delete
gunzip 0_raw_result/*.gz
## Concatenate raw results
mkdir 1_concatenated
cat 0_raw_result/*_pgenes.txt | sed '/#chr/d' | sort -k 1,1n -k 2,2n > 1_concatenated/${species}_pgenes.txt
cat 0_raw_result/*_pgenes.align > 1_concatenated/${species}_pgenes.align
## Remove raw results
rm -rf 0_raw_result

## Give each pseudogene a pseudogene number (e.g. "Pseudogene1")
awk '{printf "%s\tPseudogene%s\n",$0, NR}' 1_concatenated/${species}_pgenes.txt > tmp && mv tmp 1_concatenated/${species}_pgenes.txt

## Remove pseudogenes with overlap higher than 30 with exons (quality check)
### Convert pseudogenes txt file and exon annotation file to BED format
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$15}' 1_concatenated/${species}_pgenes.txt > tmp_pseudogenes.bed
awk -F "\t" '{print $1"\t"$4"\t"$5}' /scratch/recent_wgds/data/${species}/gff/annotation.selected_transcript.exon_features.${plaza_name}.gff3 > tmp_exons.bed
### Look for pseudogenes that have overlap (> 30 nt) with exons
bedtools intersect -a tmp_pseudogenes.bed -b tmp_exons.bed -wo | awk '{arr[$4]+=$8} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1n | awk '$2 > 30 {print $1}' > 1_concatenated/highOverlapPseudogenes.txt
sed -r 's/(.*)/\1\$/' 1_concatenated/highOverlapPseudogenes.txt > tmp && mv tmp 1_concatenated/highOverlapPseudogenes.txt
rm tmp_exons.bed tmp_pseudogenes.bed
numHighOverlap=$(cat 1_concatenated/highOverlapPseudogenes.txt | wc -l)
### Remove pseudogenes with overlap higher than 30 nt
grep -vif 1_concatenated/highOverlapPseudogenes.txt 1_concatenated/${species}_pgenes.txt > tmp && mv tmp 1_concatenated/${species}_pgenes.txt
echo "Removed $numHighOverlap pseudogenes with more than 30 nt overlap with exons!"

echo "Finished post-processing!"
