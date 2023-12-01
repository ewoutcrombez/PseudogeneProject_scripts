#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 3 # number of cores
#SBATCH --mem-per-cpu 15G # memory per core
#SBATCH -o slurm.%j_Malus_17_pp_run.out # STDOUT
#SBATCH -e slurm.%j_Malus_17_pp_run.err # STDERR

# In this script I run PseudoPipe for species as array jobs

## Input in command line
species=$1
plaza_code=$2

echo "Task ID: $SLURM_ARRAY_TASK_ID"

bash /scratch/recent_wgds/recentWGDs/scripts/helper/pseudopipe_command.sh $species $plaza_code $SLURM_ARRAY_TASK_ID
