#!/bin/bash
module load R

species=$1
plaza_name=$2

mkdir i-adhore/with_middle_filtering_level2_false/$species
echo "Get blast table"
Rscript get_blast_table.R $species
echo "Get genelists"
bash get_genelists.sh $species $plaza_name
echo "Create settingsfile"
bash settingsfile_iadhore_creation.sh $species
