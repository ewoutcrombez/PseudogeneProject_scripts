#!/bin/bash
module load R

species=$1
plaza_name=$2

echo "Working on $species"

mkdir i-adhore/no_introns2/${species}
echo "Get blast table"
Rscript get_blast_table.R $species
echo "Get genelists"
bash get_genelists.sh $species $plaza_name
echo "Create settingsfile"
bash settingsfile_iadhore_creation.sh $species
