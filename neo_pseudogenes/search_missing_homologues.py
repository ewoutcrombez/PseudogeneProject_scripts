import sys
import pandas as pd
import bioframe as bf

working_dir = sys.argv[1]
iadhore_dir = working_dir + "/results/i-ADHoRe-run"

# load in blastp results as pandas dataframe
# This consists of the BLASTP results filtered based on the translocation filters 
# (at least 60% protein percent identity and at least 70% alignment coverage)
# such as, `awk '$3 >= 60 && ($4/$13  >= 0.7) {print $0"\t"$4/$13}' allvsallblastp_filtered_wlen.tsv > allvsallblastp_transloc_filters.tsv`
print("Loading in BLASTP results...")
blastp_results = pd.read_csv(working_dir + "/results/blastp/allvsallblastp_transloc_filters.tsv", sep="\t", header=None)

# load in lonely genes file as pandas dataframe
print("Loading lonely genes file...")
lonely_genes = pd.read_csv(iadhore_dir + "/lonely_genes_all.txt", sep="\t", header=None)

# load in translocated regions file from syri as pandas dataframe
print("Loading translocated regions file...")
translocated_regions = pd.read_csv(working_dir + "/results/syri/all_trans.out", header = None, sep = "\t",
    names = ["ref_chr", "ref_start", "ref_end", "ref_seq", "q_seq", "q_chr", "q_start", "q_end", "unq_id", 
            "parent_id", "annot", "copy_stat", "unq_id_full"])
print(translocated_regions)

transloc_regions_ref = translocated_regions[["ref_chr", "ref_start", "ref_end"]]
transloc_regions_q = translocated_regions[["q_chr", "q_start", "q_end"]]
# change column names
transloc_regions_ref.columns = ["chrom", "start", "end"]
transloc_regions_q.columns = ["chrom", "start", "end"]

# create dictionary with gene name as key and chromosome and position as value
# this is used to obtain gene coordinates and next check whether a gene is in a translocated region
gene_coords = {}
with(open(iadhore_dir + "/genes.bed")) as f:
    for line in f:
        line_split = line.split("\t")
        gene_name = line_split[3]
        gene_chr = line_split[0]
        gene_start = line_split[1]
        gene_end = line_split[2]
        gene_coords[gene_name] = [gene_chr, gene_start, gene_end]

out_translocated = iadhore_dir + "/translocated_genes.tsv"
out_mult_with_missing = iadhore_dir + "/mult_with_missing.tsv"

# check if the output files already exist, if so, remove them
try:
    open(out_translocated, "x")
except FileExistsError:
    open(out_translocated, "w").close()
try:
    open(out_mult_with_missing, "x")
except FileExistsError:
    open(out_mult_with_missing, "w").close()

with(open(iadhore_dir + "/anchorpair_groups_all.tsv")) as f:
    # Go through the anchorpair groups
    for line in f:
        line_split = line.split("\t")
        # skip header
        if line_split[0] == "group_id":
            continue
        else:
            # Get the group id, genes and subgenomes it is present in
            group_id = line_split[0]
            genes = line_split[1].split(",")
            subg = line_split[2].split(",")
            subg_level = line_split[3].strip()
            
            # SEARCH FOR TRANSLOCATED GENES
            put_translocated = []
            # Find translocated homologous genes to the genes in the anchorpair group
            for gene in genes:
                # Search for all hits to the gene in the blastp results
                hit_results = blastp_results[blastp_results[0] == gene][1].tolist()
                for gene_match in hit_results:
                    # If the matching gene is part of the anchorpair group, skip it as it is not translocated
                    if gene_match in genes:
                        continue
                    # If the matching gene is not part of the anchorpair group, check if it is already in the list
                    else:
                        # Check if the gene is already in the list
                        if gene_match in put_translocated:
                            continue
                        else:
                            # Get the gene coordinates for the gene in the anchor pair group and the putative translocated gene
                            gene_coords_hit = gene_coords[gene_match]
                            gene_coords_hom = gene_coords[gene]
                            ## Convert to dataframe
                            df_put_translocated = pd.DataFrame(gene_coords_hit).T
                            df_hom_gene = pd.DataFrame(gene_coords_hom).T
                            df_put_translocated.columns = ["chrom", "start", "end"]
                            df_hom_gene.columns = ["chrom", "start", "end"]
                            ## Change start and end to integer
                            df_put_translocated.start = df_put_translocated.start.astype(int)
                            df_put_translocated.end = df_put_translocated.end.astype(int)
                            df_hom_gene.start = df_hom_gene.start.astype(int)
                            df_hom_gene.end = df_hom_gene.end.astype(int)

                            # Check if the putative translocated gene and homologous gene are in translocated regions
                            overlap_trans = bf.overlap(df_put_translocated, transloc_regions_q, how="inner")
                            overlap_hom = bf.overlap(df_hom_gene, transloc_regions_ref, how="inner")
                            
                            # Check if dataframe is empty
                            # If not empty, the gene is in a translocated region
                            if not overlap_trans.empty and not overlap_hom.empty:
                                print("Both genes are in translocated regions")
                                print(overlap_trans)
                                print(overlap_hom)
                                # loop over the rows of both dataframes and check if the gene is in the same translocated region
                                for index, row in overlap_trans.iterrows():
                                    for index2, row2 in overlap_hom.iterrows():
                                        # get the translocated region the gene is in
                                        region_trans = translocated_regions[((translocated_regions['q_chr'] == row.chrom_) & (translocated_regions['q_start'] == row.start_) & (translocated_regions['q_end'] == row.end_)) | ((translocated_regions['ref_chr'] == row.chrom_) & (translocated_regions['ref_start'] == row.start_) & (translocated_regions['ref_end'] == row.end_))]
                                        region_hom = translocated_regions[((translocated_regions['q_chr'] == row2.chrom_) & (translocated_regions['q_start'] == row2.start_) & (translocated_regions['q_end'] == row2.end_)) | ((translocated_regions['ref_chr'] == row2.chrom_) & (translocated_regions['ref_start'] == row2.start_) & (translocated_regions['ref_end'] == row2.end_))] 
                                        print(region_trans)
                                        print(region_hom)
                                        if region_trans["unq_id"].values[0] == region_hom["unq_id"].values[0]:
                                            print("The gene is in the same translocated region!!!!")
                                            put_translocated.append(gene_match)
                                        else:
                                            print("The gene is in different translocated regions")

            # If translocated genes are found, print them
            if put_translocated != []:  
                out = str(group_id) + "\t" + str(subg_level) + "\t" + str(put_translocated)
                # Write to file
                write_to_file = open(out_translocated, "a")
                write_to_file.write(out + "\n")
                write_to_file.close()

            # Check if a homologue is found in all subgenomes
            # If so, no need to search for translocated genes or collinear segments with potential pseudogenes
            if subg_level == "4":
                continue
            # If not, search for translocated genes or collinear segments with potential pseudogenes
            else:
                # SEARCH FOR MULTIPLICONS WITH A LONELY GENE AND A SUBGENOME WITHOUT HOMOLOGUE
                # Check which subgenome(s) are missing in the anchorpair group
                # by taking the difference between the set of subgenomes present and the set of all subgenomes
                missing_subg = set(["A", "B", "C", "D"]) - set(subg)
                # convert set to string
                missing_subg = list(missing_subg)
                # Check if we can find the anchorpair group genes in the lonely genes file
                for gene in genes:
                    lonely_results = lonely_genes[lonely_genes[0] == gene]
                    # Check whether it is lonely between the subgenome(s) that have no homologue in the anchorpair group
                    for missing in missing_subg:
                        # Check if one of the two subgenomes corresponds to a missing subgenome
                        lonely_between_missing_subg = lonely_results[(lonely_results[1] == missing) | (lonely_results[2] == missing)]
                        if lonely_between_missing_subg.empty:
                            continue
                        else:
                            for index, row in lonely_between_missing_subg.iterrows():
                                out = str(group_id) + "\t" + str(subg_level) + "\t" + str(row[0]) + "\t" + str(row[1]) + "\t" + str(row[2]) + "\t" + str(row[3])
                                # Write to file
                                write_to_file = open(out_mult_with_missing, "a")
                                write_to_file.write(out + "\n")
                                write_to_file.close()