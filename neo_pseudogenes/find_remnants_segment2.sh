#!/bin/bash
# Using this script I search for remnants of genes (i.e. pseudogenes) in the non-coding region

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

# HIT SEARCH

## search hits of protein in the gene-masked regions of the other segment in the collinear block
### GFF result file as output
echo "Find hits of lonely genes in the other segment!"
exonerate -m protein2genome --maxintron 6000 -q $proc_output/$multiplicon/lonely_genes_prot_${species2}.fa -t $proc_output/$multiplicon/region_${species1}_${chromosome1}_gmasked.fa --showtargetgff yes --showquerygff no --verbose 0 --showalignment no --showvulgar no --ryo %qi\\t%ti\\t%qS\\t%tS\\t%ql\\t%tab\\t%tae\\t%tal\\t%qab\\t%qae\\t%qal\\t%pi\\t%s\\t%et\\t%ei\\t%es\\t%em\\n > $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exonerate_result.txt

## Process exonerate result file
echo "Process exonerate result file!"
python3 process_exonerate_result.py $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exonerate_result.txt $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_results.tsv

# # Select the best hit for each gene if they match the same region of the gene

Rscript combine_and_get_best_hits.R $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_results.tsv $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_selected.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_results_proc.tsv

## Process GFF file by merging hits that are within 6000 bp from each other and match the same gene
echo "Process resulting GFF file!"
python3 process_gff.py $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_selected.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exons.gff

## Get hit table  (columns: hit name, hit length, hit strand, hit start, hit end, #exons)
echo "Get hit table!"
Rscript get_hit_table.R $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exons.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_hits.tsv

## GFF to BED
awk -F "\t" '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$10}' $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exons.gff > $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.bed
# make new column with the hit name
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$4}' $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.bed | sed -r 's/_exon[0-9]+$//' > tmp && mv tmp $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.bed

# sort column based on start position if hits have the same name in column 5
sort -k5,5 -k2,2n $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4}'> tmp && mv tmp $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.bed

## Get sequences of hits
echo "Get sequences of hits!"
### Obtain sequences
bedtools getfasta -fi $proc_output/$multiplicon/region_${species1}_${chromosome1}_gmasked.fa -bed $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.bed -name > $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.fa

### Combine exons together per hit --> CDS
echo "Combine exons together per hit!"
python3 combine_exons.py $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.fa $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_cds.fa

## Align sequences
### Remove files if already exists
if [ -e $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.tsv ]; then
	rm $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.tsv # if tsv file already exists, it should be removed, otherwise append to this
fi
if [ -e $proc_output/$multiplicon/hits_${species1}_in_${species2}_${chromosome2}.tsv ]; then
	rm $proc_output/$multiplicon/hits_${species1}_in_${species2}_${chromosome2}.tsv # if tsv file already exists, it should be removed, otherwise append to this
fi

if [ -s $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_cds.fa ]
then
	### Exonerate alignment
	echo "Align pseudogene/remnant with lonely genes!"
	while read lonely_gene
	do 
		samtools faidx $proc_output/$multiplicon/lonely_genes_cds_${species2}.fa $lonely_gene > $proc_output/$multiplicon/tmp_lonely_gene_cds
		grep $lonely_gene $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_cds.fa | sed 's/>//' > $proc_output/$multiplicon/tmp_hits
		samtools faidx $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_cds.fa -r $proc_output/$multiplicon/tmp_hits > $proc_output/$multiplicon/tmp_hits_cds
		if [ -s $proc_output/$multiplicon/tmp_hits_cds ]
		then
			exonerate -m affine:local $proc_output/$multiplicon/tmp_hits_cds $proc_output/$multiplicon/tmp_lonely_gene_cds --ryo %qi\\t%ti\\t%qS\\t%tS\\t%ql\\t%tl\\t%tab\\t%tae\\t%tal\\t%qab\\t%qae\\t%qal\\t%pi\\t%s\\t%et\\t%ei\\t%es\\t%em\\n --verbose 0 --showalignment no --showvulgar no >> $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.tsv
		fi
	done < $proc_output/$multiplicon/lonely_genes_segment_2.txt

	if [ -s $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.tsv ]
	then
		### Filter hits based on exonerate output
		Rscript filter_hits.R $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}.tsv $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_hits.tsv $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exons.gff $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_filtered.tsv

		if [ -s $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_filtered.tsv ]
		then
			## Get result table based on exonerate results
			Rscript get_result_table_segment.R $proc_output/$multiplicon/region_${species1}_${chromosome1}_start.txt $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_hits.tsv $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_filtered.tsv $proc_output/$multiplicon/result_table_${multiplicon}.tsv $proc_output/$multiplicon/hits_${species2}_in_${species1}_${chromosome1}_exons.gff $proc_output $multiplicon

			## Get result table based on MACSE results and classify pseudogenes
			echo "Run MACSE and classify pseudogenes!"
			python3 get_g_psg_alignments.py $working_dir $multiplicon $species2 $species1 $chromosome1

			#cat $proc_output/$multiplicon/result_table_${species1}_MACSE.tsv $proc_output/$multiplicon/result_table_${species2}_MACSE.tsv > $proc_output/$multiplicon/result_table_${multiplicon}_MACSE.tsv
			if [ -s $proc_output/$multiplicon/result_table_${multiplicon}_MACSE.tsv ]
			then
				tail -n +2 $proc_output/$multiplicon/result_table_${species2}_MACSE.tsv >> $proc_output/$multiplicon/result_table_${multiplicon}_MACSE.tsv
			else
				cp $proc_output/$multiplicon/result_table_${species2}_MACSE.tsv $proc_output/$multiplicon/result_table_${multiplicon}_MACSE.tsv
			fi
			rm $proc_output/$multiplicon/tmp_*_cds
			rm $proc_output/$multiplicon/tmp_hits
			rm $proc_output/$multiplicon/gene_seq.fasta
			rm $proc_output/$multiplicon/psg_seq_full.fasta
			rm $proc_output/$multiplicon/gene_seq_NT.fasta
		fi
	fi
else
	echo "No hits found in segment 2!"
fi


