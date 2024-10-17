#!/bin/bash
module load bedtools

working_dir=$1
species=$(basename $working_dir)

# empty all_trans.out
> $working_dir/syri/all_trans.out

for out_file in $(ls $working_dir/syri/*syri.out)
do
    subg_comp=$(basename $out_file | grep -oE "[A-D]_[A-D]")
    echo $subg_comp
    # column 11 should be TRANS or INVTR
    awk -v subg_comp=$subg_comp '$11 == "TRANS" || $11 == "INVTR" {print $0"\t"$9"-"subg_comp}' $out_file >> $working_dir/syri/all_trans.out
done

# Get translocated regions
awk '{print $1"\t"$2"\t"$3"\t"$13}' $working_dir/syri/all_trans.out > transloc_regions_q_$species.bed
awk '{print $6"\t"$7"\t"$8"\t"$13}' $working_dir/syri/all_trans.out > transloc_regions_s_$species.bed

# Get genes in translocated regions
bedtools intersect -a transloc_regions_q_$species.bed -b $working_dir/i-ADHoRe-run/genes.bed -F 1.0 -wa -wb > transloc_genes_q_$species.bed
bedtools intersect -a transloc_regions_s_$species.bed -b $working_dir/i-ADHoRe-run/genes.bed -F 1.0 -wa -wb > transloc_genes_s_$species.bed

# Find homologous genes in same group that are in translocated regions
python3 find_translocated_genes_within_groups.py $working_dir $species

# Count number of genes in translocated regions (not looking at whether it has a homolog)
num_genes=$(cat transloc_genes_q_$species.bed transloc_genes_s_$species.bed | cut -f 8 | sort | uniq | wc -l)
echo "Number of genes in translocated regions: $num_genes"

# Count number of pairs of genes in translocated regions that have a homolog in the same group
num_transloc_pairs=$(cat translocated_genes_within_groups_$species.tsv | tail -n+2 | wc -l)
echo "Number of pairs of genes in translocated regions that have a homolog in the same group: $num_transloc_pairs"

# Clean up bed files
rm transloc_regions_q_$species.bed transloc_regions_s_$species.bed transloc_genes_q_$species.bed transloc_genes_s_$species.bed

# Update result table, such that translocated genes are counted if found within a group
python3 update_result_translocation.py $working_dir