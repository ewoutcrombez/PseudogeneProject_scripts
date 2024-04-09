#!/bin/bash

# Obtain TE data from APTEdb
for species in amborella_trichopoda vitis_vinifera sorghum_bicolor oryza_sativa solanum_lycopersicum arabidopsis_thaliana populus_trichocarpa brassica_rapa brassica_oleracea zea_mays glycine_max
do
	wget http://apte.cp.utfpr.edu.br/te-annotation/$species/TEAnnotationFinal.gff3
	wget http://apte.cp.utfpr.edu.br/te-annotation/$species/TEAnnotationFinal_LTR.gff3
	wget http://apte.cp.utfpr.edu.br/te-annotation/$species/TEAnnotationFinal_LINE.gff3
	wget http://apte.cp.utfpr.edu.br/te-annotation/$species/TEAnnotationFinal_SINE.gff3
	mv TEAnnotationFinal.gff3 ${species}_TE.gff3
	grep -v "ARTEFACT" ${species}_TE.gff3 > tmp && mv tmp ${species}_TE.gff3
	mv TEAnnotationFinal_LTR.gff3 ${species}_LTR.gff3
	grep -v "ARTEFACT" ${species}_LTR.gff3 > tmp && mv tmp ${species}_LTR.gff3
	mv TEAnnotationFinal_LINE.gff3 ${species}_LINE.gff3
	grep -v "ARTEFACT" ${species}_LINE.gff3 > tmp && mv tmp ${species}_LINE.gff3
	mv TEAnnotationFinal_SINE.gff3 ${species}_SINE.gff3
	grep -v "ARTEFACT" ${species}_SINE.gff3 > tmp && mv tmp ${species}_SINE.gff3
	sort -t$'\t' -k 1,1 -k 4,4n ${species}_TE.gff3 | bedtools merge > ${species}_TE.bed
	sort -t$'\t' -k 1,1 -k 4,4n ${species}_LTR.gff3 | bedtools merge > ${species}_LTR.bed
	sort -t$'\t' -k 1,1 -k 4,4n ${species}_LINE.gff3 | bedtools merge > ${species}_LINE.bed
	sort -t$'\t' -k 1,1 -k 4,4n ${species}_SINE.gff3 | bedtools merge > ${species}_SINE.bed
done

# Malus domestica was not part of APTEdb, so TE numbers were obtained from the genome paper (DOI: 10.1038/ng.3886)

python3 count_TE_size.py