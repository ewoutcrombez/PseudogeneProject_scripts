#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem-per-cpu 8G # memory per core
#SBATCH -o slurm.%j_ath.out # STDOUT
#SBATCH -e slurm.%j_ath.err # STDERR

# In this script I (re)assess the summary statistics of the pseudogenes for species as array jobs
# Use e.g. sbatch --array=0-4 for Amborella trichopoda
# Pseudogene table is split up per ~10,000 to run in parallel

# start date of script
date
start=`date +%s`

echo "Task ID: $SLURM_ARRAY_TASK_ID"

module load KaKs_Calculator/x86_64/2.0
module load exonerate

species=$1
plaza_name=$2

python3 summary_pseudogenes.py $species $plaza_name $SLURM_ARRAY_TASK_ID

# remove temporary files
#rm gene_$species.fa pseudogene_$species.fa alignment*$species*

# end date of script
date
end=`date +%s`

runtime=$((end-start))
echo "Script executed in $runtime seconds or $(($runtime/60)) minutes!"
