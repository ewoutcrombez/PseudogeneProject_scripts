import sys
import pandas as pd

working_dir = sys.argv[1]

# load in lonely genes file as pandas dataframe
print("Loading lonely genes file...")
lonely_genes = pd.read_csv(working_dir + "/results/i-ADHoRe-run/lonely_genes_all.txt", sep="\t", header=None)

# initialize output file
with open(working_dir + "/results/i-ADHoRe-run/lonely_between_missing_subg.txt", "w") as out:
    out.write("gene\tsubg_gene\tcomparison\tgroup_id\n")
    out.close()

with(open(working_dir + "/results/i-ADHoRe-run/table_level_all.tsv")) as f:
    for line in f:
        # skip header
        if line.startswith("group_id"):
            continue
        else:
            line_split = line.split("\t")
            group_id = line_split[0] # id of the anchor pair group
            num_aps = int(line_split[1]) # number of sub-genomes where anchor pair is present
            num_psgs = int(line_split[2]) # number of sub-genomes where pseudogene is present
            num_trans = int(line_split[3]) # number of translocated genes
            num_unan = int(line_split[4]) # number of unannotated genes

            # Check which sub-genomes are represented in the anchor pair group
            subg_aps = line_split[5].split(",")
            subg_psgs = line_split[6].split(",")
            subg_unan = line_split[7].split(",")
            # Represented sub-genomes
            subg_all = subg_aps + subg_psgs + subg_unan
            # Remove "NA"s
            subg_all = [x for x in subg_all if x != "NA"]

            # Number of represented sub-genomes plus number of translocated genes
            num_rep = len(subg_all) + num_trans

            if num_rep >= 4:
                print(group_id + ": All sub-genomes represented")
            else:
                print(group_id + ": Not all sub-genomes represented")
                # Get the genes in the anchor pair group
                genes = line_split[8].split(",")
                # Check which sub-genomes are not represented
                missing_subg = [x for x in ["A", "B", "C", "D"] if x not in subg_all]

                # Check whether there are collinear segments between an anchor pair gene and the missing sub-genome(s)
                for missing in missing_subg:
                    in_collinear = False
                    for gene in genes:
                        # Check if gene is in lonely genes file
                        lonely_results = lonely_genes[lonely_genes[0] == gene]

                        # Check if one of the two subgenomes corresponds to a missing subgenome
                        lonely_between_missing_subg = lonely_results[(lonely_results[1] == missing) | (lonely_results[2] == missing)]
                        if lonely_between_missing_subg.empty:
                            continue
                        else:
                            in_collinear = True
                            # Append lonely_between_missing_subg to file
                            lonely_between_missing_subg = lonely_between_missing_subg.values.tolist()
                            with open(working_dir + "/results/i-ADHoRe-run/lonely_between_missing_subg.txt", "a") as out:
                                for row in lonely_between_missing_subg:
                                    gname = row[0]
                                    subgs = sorted(row[1:3])
                                    comparison = str(subgs[0]) + "-" + str(subgs[1])
                                    subg_gene = gname[0]
                                    out.write(gname + "\t" + subg_gene + "\t" + comparison + "\t" + str(group_id) + "\n")
                                out.close()
                    print("Missing sub-genome " + missing + " is in collinear segment: " + str(in_collinear))
                    result = group_id + "\t" + missing + "\t" + str(in_collinear) + "\n"
                    with open(working_dir + "/results/i-ADHoRe-run/collinear_missing.txt", "a") as out:
                        out.write(result)
                    out.close()