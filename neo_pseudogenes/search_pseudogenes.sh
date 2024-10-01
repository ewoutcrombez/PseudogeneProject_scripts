#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 4G # memory pool for all cores
#SBATCH -o out2/slurm_psg_search.%j.out # STDOUT
#SBATCH -e out2/slurm_psg_search.%j.err # STDERR

#set -e # error leads to script abort
set -u # abort if variable not set
set -o pipefail # error in any program leads to script abort
# Date of start script
date
start=$(date +%s)

module load exonerate/x86_64/2.4.0
module load bedtools
module load samtools
module load python/x86_64/3.6.5
module load R/x86_64/4.1.3

echo $SLURM_ARRAY_TASK_ID
# running as an array job
multiplicon=$SLURM_ARRAY_TASK_ID
# running as a single job
#multiplicon=$2

## Input in command line
working_dir=$1

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

# MASK GENIC REGIONS
echo "Masking genic region!"

## Obtain BED file for the two regions of a collinear block
Rscript obtain_bed_for_regions.R $working_dir $multiplicon $sub1 $chromosome1 $proc_output/$multiplicon/tmp_gff_$sub1 $proc_output/$multiplicon/region_${sub1}_${chromosome1}.txt $sub2 $chromosome2 $proc_output/$multiplicon/tmp_gff_$sub2 $proc_output/$multiplicon/region_${sub2}_${chromosome2}.txt

## Get the regions using `bedtools getfasta`
bedtools getfasta -fi $data_folder/genome/${sub1}_haplotype_genome.fasta -bed $proc_output/$multiplicon/region_${sub1}_${chromosome1}.bed -name -fo $proc_output/$multiplicon/region_${sub1}_${chromosome1}.fa
bedtools getfasta -fi $data_folder/genome/${sub2}_haplotype_genome.fasta -bed $proc_output/$multiplicon/region_${sub2}_${chromosome2}.bed -name -fo $proc_output/$multiplicon/region_${sub2}_${chromosome2}.fa

## Mask genic regions using `bedtools maskfasta`
bedtools maskfasta -fi $proc_output/$multiplicon/region_${sub1}_${chromosome1}.fa -bed $proc_output/$multiplicon/genes_${sub1}_${chromosome1}.bed -fo $proc_output/$multiplicon/region_${sub1}_${chromosome1}_gmasked.fa
bedtools maskfasta -fi $proc_output/$multiplicon/region_${sub2}_${chromosome2}.fa -bed $proc_output/$multiplicon/genes_${sub2}_${chromosome2}.bed -fo $proc_output/$multiplicon/region_${sub2}_${chromosome2}_gmasked.fa

## Get lonely genes that are in this multiplicon
echo "Getting lonely genes!"
awk -F "\t" -v multiplicon=$multiplicon '$6 == multiplicon {print $3}' $iadhore_output/mult_with_missing.tsv | grep "^$sub1" | sort | uniq > $proc_output/$multiplicon/lonely_genes_segment_1.txt
awk -F "\t" -v multiplicon=$multiplicon '$6 == multiplicon {print $3}' $iadhore_output/mult_with_missing.tsv | grep "^$sub2" | sort | uniq > $proc_output/$multiplicon/lonely_genes_segment_2.txt

### If lonely genes are present on that segment, get fasta sequences
echo "Getting fasta sequences for lonely genes!"
if [ -s $proc_output/$multiplicon/lonely_genes_segment_1.txt ]
then
    samtools faidx $data_folder/CDS/${sub1}.CDS.fa -r $proc_output/$multiplicon/lonely_genes_segment_1.txt > $proc_output/$multiplicon/lonely_genes_cds_${sub1}.fa
    samtools faidx $data_folder/prot/${sub1}.prot.fa -r $proc_output/$multiplicon/lonely_genes_segment_1.txt > $proc_output/$multiplicon/lonely_genes_prot_${sub1}.fa
else
    echo "No lonely genes on segment ${sub1}-${chromosome1}!"
fi

if [ -s $proc_output/$multiplicon/lonely_genes_segment_2.txt ]
then
    samtools faidx $data_folder/CDS/${sub2}.CDS.fa -r $proc_output/$multiplicon/lonely_genes_segment_2.txt > $proc_output/$multiplicon/lonely_genes_cds_${sub2}.fa
    samtools faidx $data_folder/prot/${sub2}.prot.fa -r $proc_output/$multiplicon/lonely_genes_segment_2.txt > $proc_output/$multiplicon/lonely_genes_prot_${sub2}.fa
else
    echo "No lonely genes on segment ${sub2}-${chromosome2}!"
fi

# # if there are lonely genes, search for pseudogenes
if [ -s $proc_output"/"$multiplicon"/lonely_genes_segment_1.txt" ]
then
    # find remnants of pseudogenes
    bash find_remnants_segment1.sh $multiplicon $working_dir
    # clean up intermediate files
    #bash clean_up.sh $multiplicon $working_dir
else
    echo "No lonely genes in segment 1 found!"
fi

if [ -s $proc_output"/"$multiplicon"/lonely_genes_segment_2.txt"  ]
then
    # find remnants of pseudogenes
    bash find_remnants_segment2.sh $multiplicon $working_dir
    # clean up intermediate files
    #bash clean_up.sh $multiplicon $working_dir
else
    echo "No lonely genes in segment 2 found!"
fi

# Date of finishing script
date
end=$(date +%s)
runtime=$(($end-$start))
echo "Script executed in ${runtime} seconds or $((($runtime)/60)) minutes!"