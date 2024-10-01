#!/bin/bash
# In this script the results of i-ADHoRe are processed and it is checked whether we can find additional homologous pairs in the collinear segment

## Input in command line
multiplicon=$1 # number of multiplicon/collinear segment
working_dir=$2 # working directory

## Directories to use
iadhore_input=$(echo $working_dir"/data/iadhore_input")
iadhore_output=$(echo $working_dir"/results/i-ADHoRe-run")
proc_output=$(echo $working_dir"/results/i-ADHoRe-run/processing")
data_folder=$(echo $working_dir"/data")

## Get info about multiplicon/collinear segment
sub1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $2}' $iadhore_output/multiplicons.txt)
chromosome1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $3}' $iadhore_output/multiplicons.txt)
sub2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $5}' $iadhore_output/multiplicons.txt)
chromosome2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $6}' $iadhore_output/multiplicons.txt)

begin_1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $10}' $iadhore_output/multiplicons.txt)
begin_1=$(( $begin_1 + 2 ))
end_1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $11}' $iadhore_output/multiplicons.txt)
end_1=$(( $end_1 + 2 ))
begin_2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $12}' $iadhore_output/multiplicons.txt)
begin_2=$(( $begin_2 + 2 ))
end_2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $13}' $iadhore_output/multiplicons.txt)
end_2=$(( $end_2 + 2 ))

echo "_______________________________________________________________________"
echo "Info multiplicon:"
echo " "
echo "Sub-genome 1: $sub1 $chromosome1, Sub-genome 2: $sub2 $chromosome2"
echo "Start and end of genes for sub-genome 1: $begin_1 $end_1" 
echo "Start and end of genes for sub-genome 2: $begin_2 $end_2"
echo "_______________________________________________________________________"

echo "Start processing multiplicon $multiplicon!"

## Get positions of genes in the chromosomes of the respective multiplicon
### It uses the i-ADHoRe gene lists as input which should be present in the data folder
### For sub-genome 1
python3 obtain_multiplicon_gene_positions.py $iadhore_input $sub1 $chromosome1 $proc_output/$multiplicon
### For sub-genome 2
python3 obtain_multiplicon_gene_positions.py $iadhore_input $sub2 $chromosome2 $proc_output/$multiplicon

## Select region of chromosomes that are part of the multiplicon
head -n $end_1 $proc_output/$multiplicon/${chromosome1}_${sub1}.tab.txt | tail -n +${begin_1} > $proc_output/$multiplicon/region_${sub1}_${chromosome1}.txt
head -n $end_2 $proc_output/$multiplicon/${chromosome2}_${sub2}.tab.txt | tail -n +${begin_2} > $proc_output/$multiplicon/region_${sub2}_${chromosome2}.txt

## Get multiplicon anchorpairs from i-ADHoRe output
### It uses the i-ADHoRe anchorpoints.txt file as input which should be present in the i-ADHoRe results folder
awk -F "\t" -v multiplicon=$multiplicon '$2 == multiplicon' $iadhore_output/anchorpoints.txt > $proc_output/$multiplicon/anchorpairs_$multiplicon.txt

## Get genes that are part of region of multiplicon, but are not part of an anchorpair, i.e. has no counterpart in other region
### Select anchorpair genes
cut -f4,5 $proc_output/$multiplicon/anchorpairs_$multiplicon.txt | sed 's/\t/\n/' | sort | uniq > $proc_output/$multiplicon/paired_genes.txt
## Select genes not in anchorpair genes
## For these genes we will look for additional homologous genes
grep -f $proc_output/$multiplicon/paired_genes.txt -v $proc_output/$multiplicon/region_${sub1}_${chromosome1}.txt | cut -f1 > $proc_output/$multiplicon/lonely_g_${sub1}_${chromosome1}.txt
grep -f $proc_output/$multiplicon/paired_genes.txt -v $proc_output/$multiplicon/region_${sub2}_${chromosome2}.txt | cut -f1 > $proc_output/$multiplicon/lonely_g_${sub2}_${chromosome2}.txt
cat $proc_output/$multiplicon/lonely_g_${sub1}_${chromosome1}.txt $proc_output/$multiplicon/lonely_g_${sub2}_${chromosome2}.txt > $proc_output/$multiplicon/lonely_genes.txt

# If there are no lonely genes, exit
if [ ! -s ${proc_output}/${multiplicon}/lonely_genes.txt ]
then
    echo "No lonely genes found in multiplicon!"
    exit 1 # Exit the script
fi

## Get GFF information
## Depending on the working directory, the GFF files may be in different formats
## Need to extract the gene IDs from the GFF files
if [ $working_dir = "/scratch/recent_wgds/potato" ]
then
    cat $working_dir/data/GFF/$sub1.protein-coding.gene.gff3 | sed '/^#/d' | awk -F "\t" '$3 == "mRNA"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$sub1
    cat $working_dir/data/GFF/$sub2.protein-coding.gene.gff3 | sed '/^#/d' | awk -F "\t" '$3 == "mRNA"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$sub2
elif [ $working_dir = "/scratch/recent_wgds/sugarcane" ]
then
    zcat $working_dir/data/GFF/$sub1.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$sub1
    zcat $working_dir/data/GFF/$sub2.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$sub2
elif [ $working_dir = "/scratch/recent_wgds/kiwi" ]
then
    cat $working_dir/data/GFF/$sub1.protein-coding.gene.gff3 | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$sub1
    cat $working_dir/data/GFF/$sub2.protein-coding.gene.gff3 | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$sub2
else
    echo "Error: working_dir not recognized"
    exit 1
fi
## Combine GFF information of two sub-genomes
cat $proc_output/$multiplicon/tmp_gff_$sub1 $proc_output/$multiplicon/tmp_gff_$sub2 > $proc_output/$multiplicon/tmp_gff

## Get blastp results for "lonely" genes
zcat $working_dir/results/blastp/allvsallblastp_filtered_wlen.tsv.gz | grep -f $proc_output/$multiplicon/lonely_genes.txt > $proc_output/$multiplicon/tmp_blast

## Find "lonely" genes that are tandem repeats of a paired gene or genes that have match with lonely gene on other segment or on the same chromosome but not in region
## Remove these from lonely genes list
Rscript find_tandem_or_matched.R $working_dir $multiplicon $proc_output/$multiplicon/lonely_genes.txt $chromosome1 $chromosome2 $proc_output/$multiplicon/region_${sub1}_${chromosome1}.txt $proc_output/$multiplicon/region_${sub2}_${chromosome2}.txt $sub1 $sub2

## combine multiplicon pairs and "lonely gene pair"
### Multiplicon pairs
awk -F "\t" '{print $4"\t"$5"\tanchorpair"}' $proc_output/$multiplicon/anchorpairs_${multiplicon}.txt > $proc_output/$multiplicon/pairs_multiplicon_${multiplicon}.txt
### Lonely genes that pair with each other
cat $proc_output/$multiplicon/matched_pairs.txt | awk -F "\t" '{print $1"\t"$2"\tmatched_lonely_pairs"}' >> $proc_output/$multiplicon/pairs_multiplicon_${multiplicon}.txt
cat $proc_output/$multiplicon/tandem_withhom.txt | awk -F "\t" '{print $1"\t"$2"\ttandem_pairs"}' >> $proc_output/$multiplicon/pairs_multiplicon_${multiplicon}.txt
cat $proc_output/$multiplicon/matched_in_mult_pairs.txt | awk -F "\t" '{print $1"\t"$2"\tmultiplicon_pairs"}' >> $proc_output/$multiplicon/pairs_multiplicon_${multiplicon}.txt

## Create final lonely genes file
awk -v sub1=$sub1 -v sub2=$sub2 -v mult=$multiplicon '{print $1"\t"sub1"\t"sub2"\t"mult}' $proc_output/$multiplicon/lonely_genes.txt > $proc_output/$multiplicon/lonely_genes_final.txt
