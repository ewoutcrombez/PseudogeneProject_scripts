#!/bin/bash
# In this script I will visualize the hits of the lonely genes, i.e. the remnants of duplicated genes or pseudogenes, within the regions.

## Input in command line
multiplicon=$1 # number of multiplicon
working_dir=$2 # working directory

## Directories to use
iadhore_input=$(echo $working_dir"/data/iadhore_input")
iadhore_output=$(echo $working_dir"/results/i-ADHoRe-run2")
proc_output=$(echo $working_dir"/results/i-ADHoRe-run2/processing")
data_folder=$(echo $working_dir"/data")

## Get info about multiplicon
species1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $2}' $iadhore_output/multiplicons.txt)
chromosome1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $3}' $iadhore_output/multiplicons.txt)
species2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $5}' $iadhore_output/multiplicons.txt)
chromosome2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $6}' $iadhore_output/multiplicons.txt)

begin_1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $10}' $iadhore_output/multiplicons.txt)
begin_1=$(( $begin_1 + 2 ))
end_1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $11}' $iadhore_output/multiplicons.txt)
end_1=$(( $end_1 + 2 ))
begin_2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $12}' $iadhore_output/multiplicons.txt)
begin_2=$(( $begin_2 + 2 ))
end_2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $13}' $iadhore_output/multiplicons.txt)
end_2=$(( $end_2 + 2 ))

echo "Process files to get input for GenoPlotR!"
# Add haplotype and chromosome of each hit to beginning of name
awk -F "\t" -v species1=$species1 -v chromosome1=$chromosome1 '{print species1"_"chromosome1$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_filtered_hitname.tsv > $proc_output/$multiplicon/tmp_hits_${species1}_pairs.tsv
awk -F "\t" -v species2=$species2 -v chromosome2=$chromosome2 '{print species2"_"chromosome2$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $proc_output/$multiplicon/hits_${species1}_in_${species2}_${chromosome2}_filtered_hitname.tsv > $proc_output/$multiplicon/tmp_hits_${species2}_pairs.tsv
awk -F "\t" -v species1=$species1 -v chromosome1=$chromosome1 '{print species1"_"chromosome1$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_hits_hitname.tsv > $proc_output/$multiplicon/tmp_hits_${species1}.tsv
awk -F "\t" -v species2=$species2 -v chromosome2=$chromosome2 '{print species2"_"chromosome2$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $proc_output/$multiplicon/hits_${species1}_in_${species2}_${chromosome2}_hits_hitname.tsv > $proc_output/$multiplicon/tmp_hits_${species2}.tsv

# Convert hit TSV to BED format
awk -F "\t" '{print $4"\t"$5"\t"$6"\t"$1"\t"$3}' $proc_output/$multiplicon/tmp_hits_${species1}.tsv > $proc_output/$multiplicon/tmp_hits_${species1}.bed
awk -F "\t" '{print $4"\t"$5"\t"$6"\t"$1"\t"$3}' $proc_output/$multiplicon/tmp_hits_${species2}.tsv > $proc_output/$multiplicon/tmp_hits_${species2}.bed

# Add orientation to BED gene file
Rscript add_orientation_gene.R $proc_output/$multiplicon/region_${species1}_${chromosome1}.txt $proc_output/$multiplicon/genes_${species1}_${chromosome1}.bed
Rscript add_orientation_gene.R $proc_output/$multiplicon/region_${species2}_${chromosome2}.txt $proc_output/$multiplicon/genes_${species2}_${chromosome2}.bed

# Combine genes in the region with the hits (using same coordinate system) and sort by coordinate
cat $proc_output/$multiplicon/genes_${species1}_${chromosome1}.bed $proc_output/$multiplicon/tmp_hits_${species1}.bed | sort -k2,2n -k3,3n | uniq > $proc_output/$multiplicon/${species1}_${chromosome1}_genes_hits.bed
cat $proc_output/$multiplicon/genes_${species2}_${chromosome2}.bed $proc_output/$multiplicon/tmp_hits_${species2}.bed | sort -k2,2n -k3,3n | uniq > $proc_output/$multiplicon/${species2}_${chromosome2}_genes_hits.bed

# Convert combined BED file to table as input for GenoPlotR
python3 bed_to_table.py $proc_output/$multiplicon/${species1}_${chromosome1}_genes_hits.bed $proc_output/$multiplicon $species1 $proc_output/$multiplicon/lonely_genes.txt $working_dir"/results/SNPeff/snpeff_pseudogenes.txt" $proc_output/$multiplicon/result_table_${multiplicon}_MACSE.tsv
python3 bed_to_table.py $proc_output/$multiplicon/${species2}_${chromosome2}_genes_hits.bed $proc_output/$multiplicon $species2 $proc_output/$multiplicon/lonely_genes.txt $working_dir"/results/SNPeff/snpeff_pseudogenes.txt" $proc_output/$multiplicon/result_table_${multiplicon}_MACSE.tsv

cat $proc_output/$multiplicon/${species1}_lonelygenes.txt $proc_output/$multiplicon/${species2}_lonelygenes.txt | sed '/lonely_gene/d' > $proc_output/$multiplicon/${multiplicon}_lonely_genes.txt

# Add lonely gene - hit pair to multiplicon pair table as input for GenoPlotR
cp $proc_output/$multiplicon/multiplicon_${multiplicon}_pairs.txt $proc_output/$multiplicon/multiplicon_pairs_with_pseudogenes.txt
awk -F "\t" -v species1=$species1 -v chromosome1=$chromosome1 '{print $1"\t"$2"\tgene_pseudogene_pairs"}' $proc_output/$multiplicon/tmp_hits_${species1}_pairs.tsv >> $proc_output/$multiplicon/multiplicon_pairs_with_pseudogenes.txt
awk -F "\t" -v species1=$species2 -v chromosome1=$chromosome2 '{print $1"\t"$2"\tgene_pseudogene_pairs"}' $proc_output/$multiplicon/tmp_hits_${species2}_pairs.tsv >> $proc_output/$multiplicon/multiplicon_pairs_with_pseudogenes.txt

# Visualize using GenoPlotR
echo "Creating plot (with pseudogenes/remnants)!"
Rscript visualize_i-ADHoRe_multiplicon.R $multiplicon $species1 $chromosome1 $species2 $chromosome2 $proc_output "YES"
