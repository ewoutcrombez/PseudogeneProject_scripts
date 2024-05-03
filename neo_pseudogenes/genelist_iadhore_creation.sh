#!/bin/bash
# Process GFF3 files
for file in $(ls -d GFF/*)
do
	# get name of haplotype
	name=$(echo $file | sed 's/.gff3.gz//g' | sed -r 's=GFF/==' | sed 's/.protein-coding.gene//')
	echo $name
	# make folder to store gene lists of that haplotype
	mkdir iadhore_input/genelists/$name
	# obtain list of genes
	#zcat $file | awk -F "\t" '$3 == "gene"' | sed -r 's/ID=([^;]*);.*/\1/' | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' > iadhore_input/genelists/$name/$name.lst
	# For Brassica napus, only protein-coding genes are considered
	zcat $file | awk -F "\t" '$3 == "gene"' | grep "gene_biotype=protein_coding" | sed -r 's/ID=([^;]*);.*/\1/' | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' > iadhore_input/genelists/$name/$name.lst
	# parse list of genes per chromosome in separate folder
	cd iadhore_input/genelists/$name	
       	cat $name.lst | sort -k1,1 -k2,2n | awk '{print $1"\t"$5$4}' | awk -F "\t" '{print $2>$1}'
	rm $name.lst
	cd -
done
