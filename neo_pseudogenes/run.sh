#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 6G # memory pool for all cores
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.out # STDERR

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
#multiplicon=$1

mkdir /scratch/recent_wgds/sugarcane/results/i-ADHoRe-run2/processing/$multiplicon
awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print}' /scratch/recent_wgds/sugarcane/results/i-ADHoRe-run2/multiplicons.txt > /scratch/recent_wgds/sugarcane/results/i-ADHoRe-run2/processing/$multiplicon/check_multiplicon.txt

if [ ! -s /scratch/recent_wgds/sugarcane/results/i-ADHoRe-run2/processing/$multiplicon/check_multiplicon.txt ]
then
    echo "Multiplicon does not exist!"
    exit 1
fi

## Sugarcane
working_dir="/scratch/recent_wgds/sugarcane"

# process i-ADHoRe results and filter lonely genes that still have a counterpart gene
bash process_i-ADHoRe_multiplicon.sh $multiplicon $working_dir
# if there are lonely genes, search for pseudogenes
if [ -s $working_dir"/results/i-ADHoRe-run2/processing/"$multiplicon"/lonely_genes_segment_1.txt" ]
then
    # find remnants of pseudogenes
    bash find_remnants_segment1.sh $multiplicon $working_dir
    # visualize multiplicon
    #bash visualize_with_remnants.sh $multiplicon $working_dir
    # clean up intermediate files
    #bash clean_up.sh $multiplicon $working_dir
else
    echo "No lonely genes in segment 1 found!"
fi

if [ -s $working_dir"/results/i-ADHoRe-run2/processing/"$multiplicon"/lonely_genes_segment_2.txt" ]
then
    # find remnants of pseudogenes
    bash find_remnants_segment2.sh $multiplicon $working_dir
    # visualize multiplicon
    #bash visualize_with_remnants.sh $multiplicon $working_dir
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