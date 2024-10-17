
import sys
import pandas as pd

working_dir = sys.argv[1]

# load in within-group translocated genes file as pandas dataframe
g_transloc_genes = pd.read_csv("%s/results/i-ADHoRe-run/translocated_genes_within_groups.tsv"%working_dir, sep="\t")
# change group_id to string
g_transloc_genes["group_id"] = g_transloc_genes["group_id"].astype(str)

# Initialize dictionary
dict_transloc = {}

# Loop over each group
with open("%s/results/i-ADHoRe-run/table_level_all.tsv"%working_dir, "r") as f:
    for line in f:
        # Initialize redundancy and parent group
        redundant = False
        parent_group = "NA"

        if line.startswith("group_id"):
            #initialize output file
            with open("%s/results/i-ADHoRe-run/table_level_all_upd.tsv"%working_dir, "w") as f:
                f.write(line.strip("\n") + "\tredundant\tparent_group\n")
        else:
            line_split = line.strip("\n").split("\t")
            col_group_id = line_split[0] # ID of the collinear gene group
            num_trans = int(line_split[3]) # Number of translocated genes
            genes = line_split[8].split(",") # Genes in the collinear gene group
            transloc_genes = line_split[9].strip("][").replace("'", "").split(", ") # Translocated genes in the collinear gene group

            # STEP 1: Annotate redudant groups due to between-group translocated genes
            # If there is currently a translocated gene found for the group, add it to the dictionary
            if transloc_genes != ["NA"]:
                for tg in transloc_genes:
                    dict_transloc[tg] = col_group_id # key is the translocated gene, value is the group_id
            
            # Check for all genes whether they are in the dictionary of translocated genes
            # If they are all found in the dictionary, it means that the group is redundant
            check = [gene in dict_transloc for gene in genes]
            if check == [True]*len(genes):
                redundant = True
                parent_group = dict_transloc[genes[0]]
            
            # STEP 2: Move within-group translocated genes from collinear genes to translocated genes
            within_group_transloc = g_transloc_genes[g_transloc_genes["group_id"] == col_group_id]

            while not within_group_transloc.empty:
                # Combine genes in a list
                within_transloc_genes = within_group_transloc["gene_q"].values.tolist() + within_group_transloc["gene_s"].values.tolist()
                # Select the most frequent gene
                within_transloc_gene = max(set(within_transloc_genes), key = within_transloc_genes.count)
                # Add the translocated gene to the other translocated genes
                if transloc_genes != ["NA"]:
                    transloc_genes.append(within_transloc_gene)
                    num_trans += 1
                else:
                    transloc_genes = [within_transloc_gene]
                    num_trans += 1
                # Filter out pairs where newly added gene is present
                within_group_transloc = within_group_transloc[(within_group_transloc["gene_q"] != within_transloc_gene) & (within_group_transloc["gene_s"] != within_transloc_gene)]
    
            # Write the updated translocated genes to the file
            with open("%s/results/i-ADHoRe-run/table_level_all_upd.tsv"%working_dir, "a") as f:
                f.write("\t".join(line_split[0:3]) + "\t" + str(num_trans) + "\t" + "\t".join(line_split[4:8]) + "\t" + ",".join(genes) + \
                "\t" + ",".join(transloc_genes) + "\t" + "\t".join(line_split[10:]) + "\t" + str(redundant) + "\t" + parent_group + "\n")