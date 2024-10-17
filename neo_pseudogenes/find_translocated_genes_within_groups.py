import pandas as pd
import sys

working_dir = sys.argv[1]
species = sys.argv[2]

level_file = "%s/results/i-ADHoRe-run/table_level_all.tsv"%working_dir

print("Loading translocated regions file...")
transloc_genes_q = pd.read_csv("transloc_genes_q_%s.bed"%species, header = None, sep = "\t",
    names = ["chr", "start", "end", "group_id", "g_chr", "g_start", "g_end", "gene",
            "x", "strand", "source", "type", "y", "id"])
transloc_genes_s = pd.read_csv("transloc_genes_s_%s.bed"%species, header = None, sep = "\t",
    names = ["chr", "start", "end", "group_id", "g_chr", "g_start", "g_end", "gene",
            "x", "strand", "source", "type", "y", "id"])

# initialize output file
with open("%s/results/i-ADHoRe-run/transloc_genes.txt"%working_dir, "w") as f:
    f.write("gene\tgroup_id\tcol_group_id\tchr\n")

with open("%s/results/i-ADHoRe-run/translocated_genes_within_groups.tsv"%working_dir, "w") as f:
    f.write("group_id\tgene_q\tgene_s\ttrans_id\n")

with open(level_file, "r") as f:
    for line in f:
        if line.startswith("group_id"):
            continue
        line_split = line.split("\t")
        col_group_id = line_split[0]
        genes = line_split[8].split(",")
        for gene in genes:
            # check if gene is in translocated region
            if gene in transloc_genes_q["gene"].values:
                # print row of translocated gene
                # print(transloc_genes_q[transloc_genes_q["gene"] == gene])
                # select group_id of translocated gene
                group_id = transloc_genes_q[transloc_genes_q["gene"] == gene]["group_id"].values[0]
                chr_q = transloc_genes_q[transloc_genes_q["gene"] == gene]["chr"].values[0]
                with open("%s/results/i-ADHoRe-run/transloc_genes.txt"%working_dir, "a") as f:
                    f.write(gene + "\t" + group_id + "\t" + col_group_id + "\t" + chr_q + "\n")

                # select all genes in the same group on the subject chromosome
                corresponding_genes = transloc_genes_s[(transloc_genes_s["group_id"] == group_id) & (transloc_genes_s["g_chr"] != chr_q)]
                # check if any of the corresponding genes are in the same group
                for cor_gene in corresponding_genes["gene"].values:
                    if cor_gene in genes:
                        print(str(col_group_id) + "\t" + gene + "\t" + cor_gene + "\t" + group_id)
                        result = str(col_group_id) + "\t" + gene + "\t" + cor_gene + "\t" + group_id + "\n"
                        with open("%s/results/i-ADHoRe-run/translocated_genes_within_groups.tsv"%working_dir, "a") as f:
                            f.write(result)

            if gene in transloc_genes_s["gene"].values:
                group_id = transloc_genes_s[transloc_genes_s["gene"] == gene]["group_id"].values[0]
                chr_s = transloc_genes_s[transloc_genes_s["gene"] == gene]["chr"].values[0]
                with open("%s/results/i-ADHoRe-run/transloc_genes.txt"%working_dir, "a") as f:
                    f.write(gene + "\t" + group_id + "\t" + col_group_id + "\t" + chr_s + "\n")
