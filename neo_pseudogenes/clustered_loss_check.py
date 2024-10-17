import pandas as pd
import random
import statistics
import matplotlib.pyplot as plt
import numpy as np
import sys

working_dir = sys.argv[1]

# Read in annotation of each gene (whether it has a homologue in the other subgenome or not)
annotation = pd.read_csv("%s/results/i-ADHoRe-run/lonely_between_missing_subg.txt"%working_dir, sep = "\t")
print(annotation)

# Load in GFF file
gff_A = pd.read_csv("%s/data/GFF/A.protein-coding.gene.gff3"%working_dir,
                    sep = "\t", comment = "#", header = None,
                    names = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
gff_B = pd.read_csv("%s/data/GFF/B.protein-coding.gene.gff3"%working_dir,
                    sep = "\t", comment = "#", header = None,
                    names = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
gff_C = pd.read_csv("%s/data/GFF/C.protein-coding.gene.gff3"%working_dir,
                    sep = "\t", comment = "#", header = None,
                    names = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
gff_D = pd.read_csv("%s/data/GFF/D.protein-coding.gene.gff3"%working_dir,
                    sep = "\t", comment = "#", header = None,
                    names = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

# type of interest, for potato this is "mRNA", others are "gene"
toi = "gene"

def count_neighbouring_lonely_genes(annotation, comparison, gff):
    gff = gff[gff["type"] == toi]
    # Create a column with the gene name
    gff["gene"] = gff["attributes"].str.extract("ID=([^;]+);")
    # Select the annotation for the specific comparison
    annotation = annotation[annotation["comparison"] == comparison]
    # Merge annotation and GFF
    df = pd.merge(gff, annotation, on = "gene", how = "left")

    # Loop over dataframe
    cnt_list = []
    genes_list = []
    for i in range(len(df)):
        # If gene is a lonely gene for that comparison, check 10 genes upstream and downstream
        # (if it is not lonely, the comparison will be NaN)
        if df.iloc[i]["comparison"] == comparison:
            # Check how many genes are annotated as lonely upstream and downstream of that gene
            # Check upstream
            cnt = 0
            genes = []
            genes.append(df.iloc[i]["gene"])
            j = i - 1
            while df.iloc[j]["comparison"] == comparison:
                gene = df.iloc[j]["gene"]
                genes.append(gene)
                j -= 1
                cnt += 1
            # Check downstream
            k = i + 1
            while df.iloc[k]["comparison"] == comparison:
                gene = df.iloc[k]["gene"]
                genes.append(gene)
                k += 1
                cnt += 1
            cnt_list.append(cnt)
            genes_list.append(genes)

    # Get mean of number of adjacent lonely genes
    mean_adjacent = statistics.mean(cnt_list)
    return(mean_adjacent, cnt_list, genes_list)

A_vs_B_A = count_neighbouring_lonely_genes(annotation, "A-B", gff_A)
A_vs_C_A = count_neighbouring_lonely_genes(annotation, "A-C", gff_A)
A_vs_D_A = count_neighbouring_lonely_genes(annotation, "A-D", gff_A)
A_vs_B_B = count_neighbouring_lonely_genes(annotation, "A-B", gff_B)
A_vs_C_C = count_neighbouring_lonely_genes(annotation, "A-C", gff_C)
A_vs_D_D = count_neighbouring_lonely_genes(annotation, "A-D", gff_D)
B_vs_C_B = count_neighbouring_lonely_genes(annotation, "B-C", gff_B)
B_vs_C_C = count_neighbouring_lonely_genes(annotation, "B-C", gff_C)
B_vs_D_B = count_neighbouring_lonely_genes(annotation, "B-D", gff_B)
B_vs_D_D = count_neighbouring_lonely_genes(annotation, "B-D", gff_D)
C_vs_D_C = count_neighbouring_lonely_genes(annotation, "C-D", gff_C)
C_vs_D_D = count_neighbouring_lonely_genes(annotation, "C-D", gff_D)

# plot histogram of cnt_list
# plt.hist(A_vs_B_A[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(A_vs_B_B[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(A_vs_C_A[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(A_vs_C_C[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(A_vs_D_A[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(A_vs_D_D[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(B_vs_C_B[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(B_vs_C_C[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(B_vs_D_B[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(B_vs_D_D[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(C_vs_D_C[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.hist(C_vs_D_D[1], bins = np.arange(0, 20, 1), alpha = 0.5)
# plt.show()


def count_neighbouring_lonely_genes_random(annotation, comparison, gff):
    gff = gff[gff["type"] == toi]
    # Create a column with the gene name
    gff["gene"] = gff["attributes"].str.extract("ID=([^;]+);")
    # Select the annotation for the specific comparison
    annotation = annotation[annotation["comparison"] == comparison]
    # Merge annotation and GFF
    df = pd.merge(gff, annotation, on = "gene", how = "left")

    # Get number of lonely genes
    length = len(df[df["comparison"] == comparison])

    # Get indices of genes that are not lonely
    indices = df[df["comparison"] != comparison].index.tolist()
    # create random list of indices that are not lonely
    random_num = random.sample(indices, length)
    cnt_list_random = []
    for i in random_num:
        cnt = 0
        genes = []
        j = i - 1
        if j >= 0:
            while df.iloc[j]["comparison"] == comparison:
                j -= 1
                cnt += 1
        # Check downstream
        k = i + 1
        if k < len(df):
            while df.iloc[k]["comparison"] == comparison:
                k += 1
                cnt += 1
        cnt_list_random.append(cnt)
    return(cnt_list_random)

random_mean_list = []
for x in range(1000):
    print(x)
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "A-B", gff_A)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "A-B", gff_B)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "A-C", gff_A)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "A-C", gff_C)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "A-D", gff_A)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "A-D", gff_D)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "B-C", gff_B)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "B-C", gff_C)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "B-D", gff_B)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "B-D", gff_D)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "C-D", gff_C)
    random_mean_list.append(statistics.mean(cnt_list_random))
    cnt_list_random = count_neighbouring_lonely_genes_random(annotation, "C-D", gff_D)
    random_mean_list.append(statistics.mean(cnt_list_random))

mean_real = statistics.mean([A_vs_B_A[0], A_vs_B_B[0], A_vs_C_A[0], A_vs_C_C[0], 
                            A_vs_D_A[0], A_vs_D_D[0], B_vs_C_B[0], B_vs_C_C[0], B_vs_D_B[0], 
                            B_vs_D_D[0], C_vs_D_C[0], C_vs_D_D[0]])

# plot histogram of random_mean_list and plot mean of real values
plt.hist(random_mean_list, bins = np.arange(0, 2, 0.01), alpha = 0.5)
plt.axvline(mean_real, color = "red")
plt.show()

print("Average of real number of adjacent lonely genes for a lonely gene: " + str(mean_real))
print("Average of random number of adjacent lonely genes for a lonely gene: " + str(statistics.mean(random_mean_list)))