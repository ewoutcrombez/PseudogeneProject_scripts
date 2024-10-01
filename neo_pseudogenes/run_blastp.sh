#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 8 # number of cores
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -o slurm_blastp.%j.out # STDOUT
#SBATCH -e slurm_blastp.%j.err # STDERR

set -e # error leads to script abort
set -u # abort if variable not set
set -o pipefail # error in any program leads to script abort
# Date of start script
date
start=$(date +%s)

working_dir=$1
cd $working_dir

mkdir -p results/blastp

# Load BLAST+ module
module load blast+/x86_64/2.11.0+

# Concatenate all protein sequences
cat data/prot/*.prot.fa > data/prot/all_prot.fa

# Create database
makeblastdb -in data/prot/all_prot.fa -dbtype prot -out data/prot/proteins -parse_seqids

# perform all-vs-all BLASTP search
## with lengths of query and subject length included in result table
blastp -db data/prot/proteins -query data/prot/all_prot.fa -evalue 0.0001 -out results/blastp/allvsallblastp_wlen.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads 8

# Date of finishing script
date
end=$(date +%s)
runtime=$(($end-$start))
echo "Script executed in ${runtime} seconds or $((($runtime)/60)) minutes!"