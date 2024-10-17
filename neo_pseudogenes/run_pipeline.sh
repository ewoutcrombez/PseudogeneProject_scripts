#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 10G # memory pool for all cores
#SBATCH -o slurm_pipeline.%j.out # STDOUT
#SBATCH -e slurm_pipeline.%j.err # STDERR

set -e # error leads to script abort
set -u # abort if variable not set
set -o pipefail # error in any program leads to script abort
# Date of start script
date
start=$(date +%s)

working_dir="/scratch/recent_wgds/kiwi"

# Required software
module load i-adhore/x86_64/3.0.01
module load R/x86_64/4.1.3

# STEP 0: prepare data
# 0.1: create data folder and result folder with all necessary files (see README, gene names should start with A, B, C, D indicating sub-genome, and protein names should correspond to gene names)
# 0.2: create input for i-ADHoRe
mkdir -p $working_dir/data/iadhore_input
bash genelist_iadhore_creation.sh $working_dir
bash settingsfile_iadhore_creation.sh $working_dir
sbatch run_blastp.sh $working_dir

# STEP 1: filter blastp results and run i-ADHoRe
# 1.1: filter blastp results:
## - E-val < 1e-10 
## - #aligned AAs >= 0.3*length_query 
## - abs(length_query - length_subject) <= length_query & abs(length_query - length_subject) <= length_subject 
##   (thus the difference between the length of the genes may not be bigger than the length of the query and not bigger than the length of the subject) 
## - Identity score >= 40%
awk -F "\t" 'function abs(v) {return v < 0 ? -v : v} $3 >= 40 && $4 >= 0.3*$13 && abs($13 - $14) <= $13 && abs($13 -$14) <= $14 && $11 <= 1e-10 {print}' $working_dir/results/blastp/allvsallblastp_wlen.tsv > $working_dir/results/blastp/allvsallblastp_filtered_wlen.tsv
gzip $working_dir/results/blastp/allvsallblastp_wlen.tsv
cut -f 1,2 $working_dir/results/blastp/allvsallblastp_filtered_wlen.tsv | sort -u > $working_dir/data/iadhore_input/blastp_table.txt
gzip $working_dir/results/blastp/allvsallblastp_filtered_wlen.tsv
# 1.2: run i-ADHoRe
mkdir -p $working_dir/results/i-ADHoRe-run
cd $working_dir/data/iadhore_input
# might need to change the settings of the i-ADHoRe run!!
i-adhore settingsfile_iadhore.txt

# STEP 2: Process i-ADHoRe results
# Run as array job with the number of multiplicons as the number of jobs
num_multiplicons=$(wc -l $working_dir/results/i-ADHoRe-run/multiplicons.txt | awk '{print $1-1}')
echo "Number of multiplicons: $num_multiplicons"
sbatch --array=1-${num_multiplicons}%10 process_i-ADHoRe_wrap.sh $working_dir

# STEP 3: Syri translocation detection
3.1: run syri
sbatch syri_translocation_detection.sh $working_dir

# STEP 4: Collinear retention analysis
# 4.1: get translocated regions
for out_file in $(ls $working_dir/results/syri/*syri.out)
do
    # column 11 should be TRANS or INVTR
    awk '$11 == "TRANS" || $11 == "INVTR" {print $0}' $out_file >> $working_dir/results/syri/all_trans.out
done
# 4.2: filter blastp results based on translocated gene filters
zcat $working_dir/results/blastp/allvsallblastp_filtered_wlen.tsv.gz | awk '$3 >= 60 && ($4/$13  >= 0.7) {print $0"\t"$4/$13}' > $working_dir/results/blastp/allvsallblastp_transloc_filters.tsv
# 4.3: run collinear retention analysis
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
sbatch collinear_retention_analysis.sh $working_dir $script_dir

# STEP 5: Search for pseudogenes
num_multiplicons=$(wc -l $working_dir/results/i-ADHoRe-run/multiplicons.txt | awk '{print $1-1}')
echo "Number of multiplicons: $num_multiplicons"
sbatch --array=1-${num_multiplicons}%20 search_pseudogenes.sh $working_dir

# STEP 6: Combine results and summarize
# 6.1: combine results of pseudogene search
bash get_final_pseudogene_df.sh $working_dir
Rscript combine_all_results.R $working_dir
(head -n 1 $working_dir/results/i-ADHoRe-run/table_level_all.tsv && tail -n +2 $working_dir/results/i-ADHoRe-run/table_level_all.tsv | sort -n -k2,2r -k1,1 ) > tmp && mv tmp $working_dir/results/i-ADHoRe-run/table_level_all.tsv
# 6.2: for the genes that still lack a homologous (pseudo)gene, check if they are in a collinear segment for the missing sub-genome(s)
python3 search_collinear_missing.py $working_dir
(head -n 1 $working_dir/results/i-ADHoRe-run/lonely_between_missing_subg.txt && tail -n +2 $working_dir/results/i-ADHoRe-run/lonely_between_missing_subg.txt | sort | uniq ) > tmp && mv tmp $working_dir/results/i-ADHoRe-run/lonely_between_missing_subg.txt
# 6.3: Get translocated genes (within and between groups) and update results table based on this
bash overlap_transloc_regions_genes.sh $working_dir
# 6.4: Get final results table and plot results
Rscript summarize_and_plot_results.R $working_dir
# Finished!

# Date of finishing script
date
end=$(date +%s)
runtime=$(($end-$start))
echo "Script executed in ${runtime} seconds or $((($runtime)/60)) minutes!"