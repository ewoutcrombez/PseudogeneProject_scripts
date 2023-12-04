#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem-per-cpu 2G # memory per core
#SBATCH -o slurm.%j_ath_rate_indel.out # STDOUT
#SBATCH -e slurm.%j_ath_rate_indel.err # STDERR

# start date of script
date
start=`date +%s`

module load iqtree
module load python/x86_64/3.8.0
module load java/x86_64/1.8.0_60
module load fasta

echo "Task ID: $SLURM_ARRAY_TASK_ID"

# With Arabidopsis thaliana estimated mutation rate (7 x 10-9 substitution per site per generation)
#     simulate per step of 1 MY of evolution
#     Arabidopsis thaliana mutation rate is 7 x 10-9 substitution per site per generation
#     and generation time is 1 year (annual plant), so 0.007 substitutions per site per 1 million year
#     The estimated rates of 1- to 3-bp deletions and insertions are 0.6 × 10−9 ± 0.2 × 10−9 (0.0857 deletions relative to substitution rates)
#     and 0.3 × 10−9 ± 0.1 × 10−9 per site per generation (0.0429 insertions relative to substitution rates), respectively.

python3 alisim_ath_rate_indel_tfasty.py $SLURM_ARRAY_TASK_ID

# Remove the temporary files
#rm tmp_*

# end date of script
date
end=`date +%s`

runtime=$((end-start))
echo "Script executed in $runtime seconds or $(($runtime/60)) minutes!"
