#!/bin/bash
# Using this script I obtain a plot of a selected multiplicon obtained from i-ADHoRe
# At this moment only level 2 multiplicons!

## Input in command line
species=$1
multiplicon=$2 # number of multiplicon

## Directories to use
iadhore_folder="/scratch/recent_wgds/paleo_polyploids/collinearity/i-adhore/with_middle_filtering_level2_false/$species"

## Get info about multiplicon
species1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $2}' $iadhore_folder/output/multiplicons.txt)
chromosome1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $3}' $iadhore_folder/output/multiplicons.txt)
species2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $5}' $iadhore_folder/output/multiplicons.txt)
chromosome2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $6}' $iadhore_folder/output/multiplicons.txt)

begin_1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $10}' $iadhore_folder/output/multiplicons.txt)
begin_1=$(( $begin_1 + 2 ))
end_1=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $11}' $iadhore_folder/output/multiplicons.txt)
end_1=$(( $end_1 + 2 ))
begin_2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $12}' $iadhore_folder/output/multiplicons.txt)
begin_2=$(( $begin_2 + 2 ))
end_2=$(awk -F "\t" -v multiplicon=$multiplicon '$1 == multiplicon {print $13}' $iadhore_folder/output/multiplicons.txt)
end_2=$(( $end_2 + 2 ))

echo "_______________________________________________________________________"
echo "Info multiplicon:"
echo " "
echo "Species 1: $species1 $chromosome1, Species 2: $species2 $chromosome2"
echo "Start and end for species 1: $begin_1 $end_1" 
echo "Start and end for species 2: $begin_2 $end_2"
echo "_______________________________________________________________________"

## Required modules
module load R/x86_64/4.1.3
module load python

## Get chromosomes involved in multiplicon
python3 obtain_i-ADHoRe_chrom_for_visualization.py $iadhore_folder $species1 $chromosome1
python3 obtain_i-ADHoRe_chrom_for_visualization.py $iadhore_folder $species2 $chromosome2

## Select region of chromosomes that are part of the multiplicon
head -n $end_1 ${chromosome1}_${species1}.tab.txt | tail -n +${begin_1} > region_${species1}_${chromosome1}.txt
head -n $end_2 ${chromosome2}_${species2}.tab.txt | tail -n +${begin_2} > region_${species2}_${chromosome2}.txt

## Get reverse of one of the two because this may be a better visualization (if there was a reversion in the past between the two segments)
python3 obtain_i-ADHoRe_reversed_chrom_for_visualization.py $iadhore_folder $species2 $chromosome2
len=$(wc -l ${chromosome2}_${species2}.tab_reversed.txt | cut -d " " -f 1)
head -n $(( $len - $begin_2 + 2 )) ${chromosome2}_${species2}.tab_reversed.txt | tail -n +$(( $len - $end_2 + 2 )) > region_${species2}_${chromosome2}_reversed.txt

## Get multiplicon anchorpoints
awk -F "\t" -v multiplicon=$multiplicon '$2 == multiplicon' $iadhore_folder/output/anchorpoints.txt > multiplicon_$multiplicon.txt

## Create plot
Rscript visualize_i-ADHoRe_multiplicon.R $multiplicon $species1 $chromosome1 $species2 $chromosome2

## Clean up intermediate files
rm region_* multiplicon_${multiplicon}.txt *tab*
