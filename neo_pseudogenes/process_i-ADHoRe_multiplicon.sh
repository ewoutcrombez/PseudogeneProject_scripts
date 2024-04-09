#!/bin/bash
# Using this script I process the results of i-ADHoRe and find additional homologous gene pairs
# This is a preparatory step for searching for hits in non-coding region to identify "pseudogenes" (see find_remnants.sh)

## Input in command line
multiplicon=$1 # number of multiplicon
working_dir=$2 # working directory

## Directories to use
iadhore_input=$(echo $working_dir"/data/iadhore_input")
iadhore_output=$(echo $working_dir"/results/i-ADHoRe-run2")
proc_output=$(echo $working_dir"/results/i-ADHoRe-run2/processing")
data_folder=$(echo $working_dir"/data")

# COLLINEARITY ANALYSIS

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

echo "_______________________________________________________________________"
echo "Info multiplicon:"
echo " "
echo "Species 1: $species1 $chromosome1, Species 2: $species2 $chromosome2"
echo "Start and end for species 1: $begin_1 $end_1" 
echo "Start and end for species 2: $begin_2 $end_2"
echo "_______________________________________________________________________"

echo "Start processing multiplicon $multiplicon!"

## Get chromosomes involved in multiplicon
python3 obtain_i-ADHoRe_chrom_for_visualization_between_species.py $iadhore_input $species1 $chromosome1 $proc_output/$multiplicon

python3 obtain_i-ADHoRe_chrom_for_visualization_between_species.py $iadhore_input $species2 $chromosome2 $proc_output/$multiplicon


## Select region of chromosomes that are part of the multiplicon
head -n $end_1 $proc_output/$multiplicon/${chromosome1}_${species1}.tab.txt | tail -n +${begin_1} > $proc_output/$multiplicon/region_${species1}_${chromosome1}.txt
head -n $end_2 $proc_output/$multiplicon/${chromosome2}_${species2}.tab.txt | tail -n +${begin_2} > $proc_output/$multiplicon/region_${species2}_${chromosome2}.txt

## Get reverse of one of the two because this may be a better visualization (if there was a reversion in the past between the two segments)
# python3 obtain_i-ADHoRe_reversed_chrom_for_visualization_between_species.py $iadhore_input $species2 $chromosome2
# mv ${chromosome2}_${species2}.tab_reversed.txt $proc_output/$multiplicon
# len=$(wc -l $proc_output/$multiplicon/${chromosome2}_${species2}.tab_reversed.txt | cut -d " " -f 1)
# head -n $(( $len - $begin_2 + 2 )) $proc_output/$multiplicon/${chromosome2}_${species2}.tab_reversed.txt | tail -n +$(( $len - $end_2 + 2 )) > $proc_output/$multiplicon/region_${species2}_${chromosome2}_reversed.txt

## Get multiplicon anchorpoints
awk -F "\t" -v multiplicon=$multiplicon '$2 == multiplicon' $iadhore_output/anchorpoints.txt > $proc_output/$multiplicon/multiplicon_$multiplicon.txt

## Get genes that are part of region of multiplicon, but are not part of an anchorpair, i.e. has no counterpart in other region
cut -f4,5 $proc_output/$multiplicon/multiplicon_$multiplicon.txt | sed 's/\t/\n/' | sort | uniq > $proc_output/$multiplicon/paired_genes.txt
grep -f $proc_output/$multiplicon/paired_genes.txt -v $proc_output/$multiplicon/region_${species1}_${chromosome1}.txt | cut -f1 > $proc_output/$multiplicon/lonely_g_${species1}_${chromosome1}.txt
grep -f $proc_output/$multiplicon/paired_genes.txt -v $proc_output/$multiplicon/region_${species2}_${chromosome2}.txt | cut -f1 > $proc_output/$multiplicon/lonely_g_${species2}_${chromosome2}.txt
cat $proc_output/$multiplicon/lonely_g_${species1}_${chromosome1}.txt $proc_output/$multiplicon/lonely_g_${species2}_${chromosome2}.txt > $proc_output/$multiplicon/lonely_genes.txt

# If there are no lonely genes, exit
if [ ! -s ${proc_output}/${multiplicon}/lonely_genes.txt ]
then
    echo "No lonely genes!"
    exit 1
fi

## Select blastp results of "lonely" genes
#zcat ../results/blastp/otava_allvsallblastp_wlen.tsv.gz | grep -f $proc_output/$multiplicon/lonely_genes.txt > $proc_output/$multiplicon/tmp_blast

if [ $working_dir = "/scratch/recent_wgds/otava_pseudogenes/OtavaPseudogenes" ]
then
    zcat $working_dir/data/GFF/$species1.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "mRNA"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$species1
    zcat $working_dir/data/GFF/$species2.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "mRNA"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$species2
elif [ $working_dir = "/scratch/recent_wgds/sugarcane" ]
then
    zcat $working_dir/data/GFF/$species1.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$species1
    zcat $working_dir/data/GFF/$species2.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$species2
elif [ $working_dir = "/scratch/recent_wgds/napus" ]
then
    zcat $working_dir/data/GFF/$species1.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | grep "gene_biotype=protein_coding" | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$species1
    zcat $working_dir/data/GFF/$species2.protein-coding.gene.gff3.gz | sed '/^#/d' | awk -F "\t" '$3 == "gene"' | grep "gene_biotype=protein_coding" | sed -r 's/ID=([^;]*);.*/\1/' > $proc_output/$multiplicon/tmp_gff_$species2
else
    echo "Error: working_dir not recognized"
    exit 1
fi

zcat $working_dir/results/blastp/allvsallblastp_filtered_wlen.tsv.gz | grep -f $proc_output/$multiplicon/lonely_genes.txt > $proc_output/$multiplicon/tmp_blast

cat $proc_output/$multiplicon/tmp_gff_$species1 $proc_output/$multiplicon/tmp_gff_$species2 > $proc_output/$multiplicon/tmp_gff

## Find "lonely" genes that are tandem repeats of a paired gene or genes that have match with lonely gene on other segment or on the same chromosome but not in region
## Remove these from lonely genes list
Rscript find_tandem_or_matched.R $working_dir $multiplicon $proc_output/$multiplicon/lonely_genes.txt $chromosome1 $chromosome2 $proc_output/$multiplicon/region_${species1}_${chromosome1}.txt $proc_output/$multiplicon/region_${species2}_${chromosome2}.txt $proc_output/$multiplicon/region_${species2}_${chromosome2}_reversed.txt $species1 $species2

## combine multiplicon pairs and "lonely gene pair"
### Multiplicon pairs
awk -F "\t" '{print $4"\t"$5"\tanchorpair"}' $proc_output/$multiplicon/multiplicon_${multiplicon}.txt > $proc_output/$multiplicon/multiplicon_${multiplicon}_pairs.txt
### Lonely genes that pair with each other
cat $proc_output/$multiplicon/matched_pairs.txt | awk -F "\t" '{print $1"\t"$2"\tmatched_lonely_pairs"}' >> $proc_output/$multiplicon/multiplicon_${multiplicon}_pairs.txt
cat $proc_output/$multiplicon/tandem_withhom.txt | awk -F "\t" '{print $1"\t"$2"\ttandem_pairs"}' >> $proc_output/$multiplicon/multiplicon_${multiplicon}_pairs.txt
cat $proc_output/$multiplicon/matched_in_mult_pairs.txt | awk -F "\t" '{print $1"\t"$2"\tmultiplicon_pairs"}' >> $proc_output/$multiplicon/multiplicon_${multiplicon}_pairs.txt

cat $proc_output/$multiplicon/matched_chrom_genes.txt | awk -F "\t" '{print $1"\t"$2"\tmatched_in_chrom"}' > $proc_output/$multiplicon/${multiplicon}_matched_outside_mult.txt
cat $proc_output/$multiplicon/matched_scat_genes.txt | awk -F "\t" '{print $1"\t"$2"\tmatched_scattered"}' >> $proc_output/$multiplicon/${multiplicon}_matched_outside_mult.txt

## Create plot
#echo "Creating plot (without pseudogenes/remnants)!"
#Rscript visualize_i-ADHoRe_multiplicon.R $multiplicon $species1 $chromosome1 $species2 $chromosome2 $proc_output "NO"

# MASK GENIC REGIONS
echo "Masking genic region!"

## Obtain BED file for the two regions of a collinear block
Rscript obtain_bed_for_regions.R $working_dir $multiplicon $species1 $chromosome1 $proc_output/$multiplicon/tmp_gff_$species1 $proc_output/$multiplicon/region_${species1}_${chromosome1}.txt $species2 $chromosome2 $proc_output/$multiplicon/tmp_gff_$species2 $proc_output/$multiplicon/region_${species2}_${chromosome2}.txt

## Get the regions using `bedtools getfasta`
bedtools getfasta -fi $data_folder/genome/${species1}_haplotype_genome.fasta -bed $proc_output/$multiplicon/region_${species1}_${chromosome1}.bed -name -fo $proc_output/$multiplicon/region_${species1}_${chromosome1}.fa
bedtools getfasta -fi $data_folder/genome/${species2}_haplotype_genome.fasta -bed $proc_output/$multiplicon/region_${species2}_${chromosome2}.bed -name -fo $proc_output/$multiplicon/region_${species2}_${chromosome2}.fa

## Mask genic regions using `bedtools maskfasta`
bedtools maskfasta -fi $proc_output/$multiplicon/region_${species1}_${chromosome1}.fa -bed $proc_output/$multiplicon/genes_${species1}_${chromosome1}.bed -fo $proc_output/$multiplicon/region_${species1}_${chromosome1}_gmasked.fa
bedtools maskfasta -fi $proc_output/$multiplicon/region_${species2}_${chromosome2}.fa -bed $proc_output/$multiplicon/genes_${species2}_${chromosome2}.bed -fo $proc_output/$multiplicon/region_${species2}_${chromosome2}_gmasked.fa

## Get CDS and prot sequences of the lonely genes
### If lonely genes are present on that segment, get fasta sequences
if [ -s $proc_output/$multiplicon/lonely_g_${species1}_${chromosome1}.txt ]
then
    samtools faidx $data_folder/CDS/${species1}.CDS.fa -r $proc_output/$multiplicon/lonely_g_${species1}_${chromosome1}.txt > $proc_output/$multiplicon/lonely_genes_cds_${species1}.fa
    samtools faidx $data_folder/prot/${species1}.prot.fa -r $proc_output/$multiplicon/lonely_g_${species1}_${chromosome1}.txt > $proc_output/$multiplicon/lonely_genes_prot_${species1}.fa
fi

if [ -s $proc_output/$multiplicon/lonely_g_${species2}_${chromosome2}.txt ]
then
    samtools faidx $data_folder/CDS/${species2}.CDS.fa -r $proc_output/$multiplicon/lonely_g_${species2}_${chromosome2}.txt > $proc_output/$multiplicon/lonely_genes_cds_${species2}.fa
    samtools faidx $data_folder/prot/${species2}.prot.fa -r $proc_output/$multiplicon/lonely_g_${species2}_${chromosome2}.txt > $proc_output/$multiplicon/lonely_genes_prot_${species2}.fa
else
    echo "No lonely genes on segment ${species2}-${chromosome2}!"
fi

cp $proc_output/$multiplicon/lonely_g_${species1}_${chromosome1}.txt $proc_output/$multiplicon/lonely_genes_segment_1.txt
cp $proc_output/$multiplicon/lonely_g_${species2}_${chromosome2}.txt $proc_output/$multiplicon/lonely_genes_segment_2.txt