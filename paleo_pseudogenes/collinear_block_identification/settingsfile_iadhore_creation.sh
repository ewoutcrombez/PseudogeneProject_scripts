#!/bin/bash

species=$1

cd i-adhore/no_introns2/${species}/genome_lists
echo "genome=$species" >> ../settingsfile_iadhore.txt
# Create settingsfile
for chrom in $(ls)
do
	echo $chrom
	path=$(realpath $chrom)
	echo "$chrom $path" >> ../settingsfile_iadhore.txt
done
cd -
cat settingsiADHoRe.txt >> i-adhore/no_introns2/${species}/settingsfile_iadhore.txt
