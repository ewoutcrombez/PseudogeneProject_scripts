# In this script I clean up the files that are not needed anymore

## Input in command line
multiplicon=$1 # number of multiplicon
working_dir=$2 # working directory

## Directories to use
iadhore_input=$(echo $working_dir"/data/iadhore_input")
iadhore_output=$(echo $working_dir"/results/i-ADHoRe-run")
proc_output=$(echo $working_dir"/results/i-ADHoRe-run/processing")
data_folder=$(echo $working_dir"/data")

## Get info about multiplicon
species1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $2}' $iadhore_output/multiplicons.txt)
chromosome1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $3}' $iadhore_output/multiplicons.txt)
species2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $5}' $iadhore_output/multiplicons.txt)
chromosome2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $6}' $iadhore_output/multiplicons.txt)

## Remove files that are not needed anymore
cd $proc_output/$multiplicon
rm tmp* *.tab.txt *_hitname.tsv *_results* 
rm matched_chrom_table.txt *_selected.gff *_cds.fa* 
rm hits_${species1}_in_${species2}_${chromosome2}.*
rm hits_${species2}_in_${species1}_${chromosome1}.*
rm multiplicon_${multiplicon}.txt multiplicon_${multiplicon}_pairs.txt
rm gene_seq_AA.fasta gene_seq_NT_NT.fasta gene_seq_NT_AA.fasta
rm lonely*
