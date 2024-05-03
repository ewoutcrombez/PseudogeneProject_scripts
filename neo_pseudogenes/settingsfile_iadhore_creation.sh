#!/bin/bash
cd iadhore_input/genelists

# Create settingsfile
for haplotype in $(ls)
do
	echo "genome=$haplotype" >> ../settingsfile_iadhore.txt
	cd $haplotype
	for file in $(ls)
	do
		path=$(realpath $file)
		echo "$file $path" >> ../../settingsfile_iadhore.txt
	done
	cd -
done
cat /scratch/recent_wgds/general_scripts/settingsiADHoRe.txt >> ../settingsfile_iadhore.txt
