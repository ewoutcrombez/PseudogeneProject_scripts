#!/bin/bash

# Required modules
module load R/x86_64/4.1.3

working_dir=$1

for multiplicon in $(ls $working_dir/results/i-ADHoRe-run/processing/)
do
    echo $multiplicon
    # if there are lonely genes, search for pseudogenes
    if [ -s $working_dir"/results/i-ADHoRe-run/processing/"$multiplicon"/result_table_${multiplicon}_MACSE.tsv" ]
    then
        tail -n+2 $working_dir"/results/i-ADHoRe-run/processing/"$multiplicon"/result_table_${multiplicon}_MACSE.tsv" | awk -v mult=$multiplicon '{print $0"\t"mult}' >> $working_dir"/results/i-ADHoRe-run/result_table_MACSE.tsv"
        cat $working_dir"/results/i-ADHoRe-run/processing/"$multiplicon"/result_table_${multiplicon}.tsv" | awk -v mult=$multiplicon '{print $0"\t"mult}' >> $working_dir"/results/i-ADHoRe-run/result_table.tsv"
    else
        echo "No hits found for multiplicon ${multiplicon}!"
    fi
done

# Filter the result table
Rscript get_final_pseudogene_df.R $working_dir