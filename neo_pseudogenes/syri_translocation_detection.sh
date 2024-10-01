#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 12 # number of cores
#SBATCH --mem 160G # memory pool for all cores
#SBATCH -o slurm_syri.%j.out # STDOUT
#SBATCH -e slurm_syri.%j.err # STDERR

set -e # error leads to script abort
set -u # abort if variable not set
set -o pipefail # error in any program leads to script abort
# Date of start script
date
start=$(date +%s)

# Required modules
module load minimap/x86_64/2.19
module load samtools/x86_64/1.20
module load gcc
module load syri/x86_64/1.6.3
module load plotsr/x86_64/0.5.4
module load MUMmer/x86_64/4.0.0beta2

working_dir=$1

# Genome directory
genome_dir="$working_dir/data/genome"
result_dir="$working_dir/results/syri"
gff_dir="$working_dir/data/GFF"

mkdir -p $result_dir
cd $result_dir

for genome1 in A B C D
do
    for genome2 in A B C D
    do
        if [[ $genome1 != $genome2 && ! -f ${genome2}_${genome1}.bam && ! -f ${genome1}_${genome2}.bam ]]
        then
            echo "$genome1 vs $genome2"
            # Step 1: Align the subgenomes
            minimap2 -ax asm5 --eqx -t 12 $genome_dir/${genome1}_haplotype_genome.fasta $genome_dir/${genome2}_haplotype_genome.fasta > ${genome1}_${genome2}.sam
            samtools sort -o ${genome1}_${genome2}.bam ${genome1}_${genome2}.sam
            rm ${genome1}_${genome2}.sam
            samtools index ${genome1}_${genome2}.bam

            # Step 2: Finding structural annotations between genomes
            syri -c ${genome1}_${genome2}.bam -r $genome_dir/${genome1}_haplotype_genome.fasta -q $genome_dir/${genome2}_haplotype_genome.fasta -F B --prefix ${genome1}_${genome2}
            
            # # Step 3: Running plotsr
            # ## Create a file with the genomes
            # echo -e "#file\tname\ttags" > genomes.txt
            # echo -e "$genome_dir/${genome1}_haplotype_genome.fasta\t$genome1\tlw:1.5" >> genomes.txt
            # echo -e "$genome_dir/${genome2}_haplotype_genome.fasta\t$genome2\tlw:1.5" >> genomes.txt
            # ## Create a file with the tracks
            # echo -e "#file\tname\ttags" > tracks.txt
            # echo -e "$gff_dir/${genome1}.protein-coding.gene.gff3\tGenes\tft:gff" >> tracks.txt
            # ## Run plotsr
            # plotsr \
            #     --sr ${genome1}_${genome2}syri.out \
            #     --genomes genomes.txt \
            #     --tracks tracks.txt
            #     -o output_plot_${genome1}_${genome2}.png
        fi
    done
done

# Date of finishing script
date
end=$(date +%s)
runtime=$(($end-$start))
echo "Script executed in ${runtime} seconds or $((($runtime)/60)) minutes!"
