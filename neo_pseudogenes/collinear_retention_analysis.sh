#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 15G # memory pool for all cores
#SBATCH -o slurm_ret_analysis.%j.out # STDOUT
#SBATCH -e slurm_ret_analysis.%j.err # STDERR

set -e # error leads to script abort
set -u # abort if variable not set
set -o pipefail # error in any program leads to script abort
# Date of start script
date
start=$(date +%s)

module load bedops/x86_64/2.4.37
module load R/x86_64/4.1.3
module load python/x86_64/3.8.0

working_dir=$1
script_dir=$2
echo $script_dir

cd $working_dir/results/i-ADHoRe-run
# Combine results of i-ADHoRe processing of all multiplicons/collinear segments
for folder in $(ls processing)
do
    echo $folder
    if [ -f processing/$folder/pairs_multiplicon_${folder}.txt ]
    then 
        awk -v mult=$folder '{print $1"\t"$2"\t"mult}' processing/$folder/pairs_multiplicon_${folder}.txt >> hom_pairs_all.tsv
    elif [ -f processing/$folder/anchorpairs_${folder}.txt ]
    then 
        echo "No additional homologous pairs were found for this multiplicon so using only the i-ADHoRe anchor pairs"
        awk -v mult=$folder '{print $4"\t"$5"\t"mult}' processing/$folder/anchorpairs_${folder}.txt >> hom_pairs_all.tsv
    else 
        echo "ERROR: No homologous pairs found for this multiplicon!"
        exit 1 # Exit the script
    fi
done

# Get a list of all genes
echo "Getting a list of all genes..."
cd $working_dir/data/iadhore_input/genelists
# get all files in gene lists folders and concatenate them, then remove the last character of each line
cat */* | sed 's/.$//' | sort | uniq > $working_dir/results/i-ADHoRe-run/all_genes.txt

# Get a list of all lonely genes linked to the multiplicon
echo "Getting a list of all lonely genes linked to multiplicon..."
cd $working_dir/results/i-ADHoRe-run/processing
cat */lonely_genes_final.txt > ../lonely_genes_all.txt

cd $script_dir
# Assess level of collinearity retention across the different sub-genomes
echo "Starting collinear retention analysis..."
Rscript collinear_retention_analysis.R $working_dir

cd $working_dir/results/i-ADHoRe-run
# Get file with all locations of genes in the genome (all sub-genomes)
echo "Getting file with all locations of genes in the genome..."
# for potato better to work with mRNA names
if [ $working_dir = "/scratch/recent_wgds/potato" ]
then
    gff2bed < $working_dir/data/GFF/A.protein-coding.gene.gff3 | awk '$8 == "mRNA"' > genes.bed
    gff2bed < $working_dir/data/GFF/B.protein-coding.gene.gff3 | awk '$8 == "mRNA"' >> genes.bed
    gff2bed < $working_dir/data/GFF/C.protein-coding.gene.gff3 | awk '$8 == "mRNA"' >> genes.bed
    gff2bed < $working_dir/data/GFF/D.protein-coding.gene.gff3 | awk '$8 == "mRNA"' >> genes.bed
else
    gff2bed < $working_dir/data/GFF/A.protein-coding.gene.gff3 | awk '$8 == "gene"' > genes.bed
    gff2bed < $working_dir/data/GFF/B.protein-coding.gene.gff3 | awk '$8 == "gene"' >> genes.bed
    gff2bed < $working_dir/data/GFF/C.protein-coding.gene.gff3 | awk '$8 == "gene"' >> genes.bed
    gff2bed < $working_dir/data/GFF/D.protein-coding.gene.gff3 | awk '$8 == "gene"' >> genes.bed
fi

# Search for the collinear segments that have lonely genes and that correspond to genes that are not present in the other sub-genome
cd $script_dir
echo "Starting search for translocated genes and potential regions for pseudogenes..."
python3 search_missing_homologues.py $working_dir

# Date of finishing script
date
end=$(date +%s)
runtime=$(($end-$start))
echo "Script executed in ${runtime} seconds or $((($runtime)/60)) minutes!"