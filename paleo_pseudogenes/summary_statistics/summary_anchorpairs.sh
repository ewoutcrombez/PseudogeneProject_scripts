#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 6G # memory pool for all cores
#SBATCH -o slurm.%j_anchorpoint_pid.out # STDOUT
#SBATCH -e slurm.%j_anchorpoint_pid.err # STDERR

# In this script I assess the divergence between the paralogous anchorpair genes
# First I need to split up the anchorpoint file in smaller files to run in parallel
# To do this I use "split_file $dir/anchorpoints.txt 10"
# Pseudogene anchorpoints need to be removed (for file in $(ls anchorpoints.txt.*); do echo $file; grep -v "Pseudogene" $file > tmp && mv tmp $file ; done)
# Then I use "sbatch --array=0-9 summary_anchorpairs.sh" to run the script in parallel

# start date of script
date
start=`date +%s`

module load KaKs_Calculator/x86_64/2.0

echo "Task ID: $SLURM_ARRAY_TASK_ID"

species=$1
plaza_name=$2

python3 summary_anchorpairs.py $species $plaza_name $SLURM_ARRAY_TASK_ID

# remove temporary files
#rm tmp*

# end date of script
date
end=`date +%s`

runtime=$((end-start))
echo "Script executed in $runtime seconds or $(($runtime/60)) minutes!"
