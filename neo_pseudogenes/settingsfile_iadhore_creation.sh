#!/bin/bash

working_dir=$1

cd $working_dir/data/iadhore_input/genelists

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
	cd $working_dir/data/iadhore_input/genelists
done
cat /scratch/recent_wgds/general_scripts/settingsiADHoRe.txt >> ../settingsfile_iadhore.txt