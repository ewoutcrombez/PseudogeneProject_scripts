#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 4 # number of cores
#SBATCH --mem-per-cpu 10G # memory per core
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

module load blast+

species=$1
species_short=$2

# create database
makeblastdb -in /scratch/recent_wgds/data/$species/fasta/proteome.selected_transcript.${species_short}.fasta -dbtype prot -out /scratch/recent_wgds/data/$species/fasta/${species}_prot -parse_seqids

# perform all-vs-all BLASTP search
## with lengths of query and subject length included in result table
blastp -db /scratch/recent_wgds/data/$species/fasta/${species}_prot -query /scratch/recent_wgds/data/$species/fasta/proteome.selected_transcript.${species_short}.fasta -evalue 0.0001 -out ${species}_allvsallblastp_wlen.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads 8

# filter
awk -F "\t" \
  'function abs(v) {return v < 0 ? -v : v} $4 >= 0.3*$13 && abs($13 - $14) <= $13 && abs($13 -$14) <= $14 && $11 <= 1e-5 {print}' ${species}_allvsallblastp_wlen.tsv > ${species}_allvsallblastp_filtered_wlen.tsv 
