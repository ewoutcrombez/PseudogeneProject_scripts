#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 2G # memory pool for all cores
#SBATCH -o out2/slurm_process.%j.out # STDOUT
#SBATCH -e out2/slurm_process.%j.err # STDERR

# This script is used to go over the i-ADHoRe multiplicons and process them to find additional homologous within the same multiplicon
# It can be run as an array job or as a single job whereby the multiplicon number is passed as an argument

set -e # error leads to script abort
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
# multiplicon=$1

## working directory should contain a data and results folder
## data folder should contain the genome, GFF and i-ADHoRe input folder
## results folder should contain the i-ADHoRe results folder
working_dir=$1
# working_dir=$2

res_dir="$working_dir/results/i-ADHoRe-run"

mkdir -p $res_dir/processing/$multiplicon
awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print}' $res_dir/multiplicons.txt > $res_dir/processing/$multiplicon/check_multiplicon.txt

if [ ! -s $res_dir/processing/$multiplicon/check_multiplicon.txt ]
then
    echo "Multiplicon does not exist!"
    exit 1 # Exit the script
fi

# process i-ADHoRe results and find additional homologous genes within the same multiplicon
bash process_i-ADHoRe_multiplicon.sh $multiplicon $working_dir

# Date of finishing script
date
end=$(date +%s)
runtime=$(($end-$start))
echo "Script executed in ${runtime} seconds or $((($runtime)/60)) minutes!"
