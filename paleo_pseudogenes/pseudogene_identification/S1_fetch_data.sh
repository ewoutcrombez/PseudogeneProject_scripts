#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 2 # number of cores
#SBATCH --mem 30G # memory pool for all cores
#SBATCH -o slurm.%j_fetch.out # STDOUT
#SBATCH -e slurm.%j_fetch.err # STDERR
# In this script I fetch the required data for the pseudogene analysis from PLAZA and do an initial processing on these files

## Required modules
module load python
module load RepeatMasker

## Input in command line
species=$1 # e.g. Arabidopsis_thaliana
plaza_name=$2 # species abbreviation used in PLAZA, e.g. ath
species_RepeatMasker=$3 # species abbreviation used in RepeatMasker, e.g. Arabidopsis
clade=$4 # monocot or dicot? -> fill in "monocots" or "dicots"

### Create data directory for species
mkdir /scratch/recent_wgds/data/$species
cd /scratch/recent_wgds/data/$species

### Download fasta files
echo "Downloading fasta files..."
mkdir fasta
cd fasta

### Full dna fasta file
wget --no-check-certificate https://ftp.psb.ugent.be/pub/plaza/plaza_public_${clade}_05/Genomes/${plaza_name}.fasta.gz
file=$(echo "${plaza_name}.fasta.gz")
gunzip --force $file
file=$(echo "$file" | sed 's/.gz//')
### Repeat mask dna file
RepeatMasker -norna -species $species_RepeatMasker $file

### Split per fasta header (thus per chromosome/scaffold/...) for both repeat masked and unmasked dna fasta file
### This is needed for the input of PseudoPipe
mkdir ../ppipe_input
cd ../ppipe_input
awk '/^>/{filename=substr($1,2)".fa"}{print > filename}' ../fasta/$file
awk '/^>/{filename=substr($1,2)".fa"}{print > "dna_rm"filename}' ../fasta/${file}.masked
for RMfile in $(ls dna_rm*.fa)
do
    name=$(echo "$RMfile" | sed 's/.fa//' | sed 's/dna_rm//')
    mkdir -p $name/dna
    mkdir -p $name/mysql
    mkdir -p $name/pep
    mv $RMfile $name/dna/dna_rm.fa
    mv $name".fa" $name/dna
done

# ## cds fasta file
cd ../fasta
wget --no-check-certificate https://ftp.psb.ugent.be/pub/plaza/plaza_public_${clade}_05/Fasta/cds.selected_transcript.${plaza_name}.fasta.gz
file=$(echo "cds.selected_transcript.${plaza_name}.fasta.gz")
gunzip --force $file
file=$(echo "$file" | sed 's/.gz//')

# ## protein fasta file
wget --no-check-certificate https://ftp.psb.ugent.be/pub/plaza/plaza_public_${clade}_05/Fasta/proteome.selected_transcript.${plaza_name}.fasta.gz
file=$(echo "proteome.selected_transcript.${plaza_name}.fasta.gz")
gunzip --force $file
file=$(echo "$file" | sed 's/.gz//')
### Copy the proteome file to every PseudoPipe input pep directory
for folder in $(ls ../ppipe_input)
do
    cp $file ../ppipe_input/$folder/pep/$file
done

## Download gff files
echo "Downloading gff files..."
mkdir ../gff
cd ../gff

## exons gff file
wget --no-check-certificate https://ftp.psb.ugent.be/pub/plaza/plaza_public_${clade}_05/GFF/${plaza_name}/annotation.selected_transcript.exon_features.${plaza_name}.gff3.gz
file=$(echo "annotation.selected_transcript.exon_features.${plaza_name}.gff3.gz")
gunzip --force $file
file=$(echo "$file" | sed 's/.gz//')
### Remove headers starting with hashtag
cat $file | sed '/^#/d' > tmp && mv tmp $file
### Remove features annotated with "pseudogene" or are "non-coding RNA"
awk -F "\t" '$3 != "pseudogene" {print}' $file > tmp && mv tmp $file
#### Remember which genes are annotated with small RNAs
awk -F "\t" '$3 == "lncRNA" || $3 == "miRNA" || $3 == "ncRNA" || $3 == "snoRNA" || $3 == "snRNA" || $3 == "antisense lncRNA" || $3 == "novel_transcribed_region" || $3 == "tRNA" || $3 == "rRNA" {print}' $file | sed -r 's/.*gene_id=(.*)/\1/' > smallRNAs.txt
#### And remove these from the gff file (required for i-ADHoRe)
grep -vif smallRNAs.txt $file > tmp && mv tmp $file
### Select only exons
awk -F "\t" '$3 == "exon" {print}' $file > tmp && mv tmp $file
### Select only coding exons (CDS)
awk -F "\t" '$3 == "CDS" {print}' $file > tmp && mv tmp $file
### Create exLocs files required for PseudoPipe
cd ../ppipe_input
awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5 > $1"_exLocs"}' ../gff/$file
for file in $(ls *_exLocs)
do
     name=$(echo "$file" | sed 's/_exLocs//')
     mv $file $name/mysql
done

cd ../gff
## all features gff file
wget --no-check-certificate https://ftp.psb.ugent.be/pub/plaza/plaza_public_${clade}_05/GFF/${plaza_name}/annotation.selected_transcript.all_features.${plaza_name}.gff3.gz
file=$(echo "annotation.selected_transcript.all_features.${plaza_name}.gff3.gz")
gunzip --force $file
file=$(echo "$file" | sed 's/.gz//')
### Remove headers starting with hashtag
cat $file | sed '/^#/d' > tmp && mv tmp $file
### Remove features annotated with "pseudogene" or are "non-coding RNA"
awk -F "\t" '$3 != "pseudogene" && $3 != "pseudogenic_transcript" && $3 != "pseudogenic_exon" {print}' $file > tmp && mv tmp $file
#### Remember which genes are annotated with small RNAs
awk -F "\t" '$3 == "lncRNA" || $3 == "miRNA" || $3 == "ncRNA" || $3 == "snoRNA" || $3 == "snRNA" || $3 == "antisense lncRNA" || $3 == "novel_transcribed_region" || $3 == "tRNA" || $3 == "rRNA" || $3 == "antisense RNA" {print}' $file | sed -r 's/.*gene_id=(.*)/\1/' > smallRNAs.txt
#### And remove these from the gff file (required for i-ADHoRe)
grep -vif smallRNAs.txt $file > tmp && mv tmp $file
rm smallRNAs.txt

awk -F "\t" '$3 == "CDS" {print}' $file > tmp && mv tmp annotation.selected_transcript.exon_features.${plaza_name}.gff3
### Create exLocs files required for PseudoPipe
cd ../ppipe_input
awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5 > $1"_exLocs"}' ../gff/annotation.selected_transcript.exon_features.${plaza_name}.gff3
for file in $(ls *_exLocs)
do
    name=$(echo "$file" | sed 's/_exLocs//')
    mv $file $name/mysql
done

echo "Fetched all data! Hooray!"