# In this script I check the overlap between lncRNAs and pseudogenes and
# reshuffle lncRNAs randomly and subsequently calculate the overlap between the pseudogenes and the random lncRNAs
# The strategy is to randomly shuffle the lncRNA many times, each time doing an intersection with another file of interest and counting the number of 
# intersections (or computing some other statistic on the overlap). Upon doing this many times, an empirical distribution is constructed, and the
# number of intersections between the original, un-shuffled file is compared to this empirical distribution to obtain a p-value
# SEE: https://daler.github.io/pybedtools/topical-random.html

import sys
import os
import pybedtools
import matplotlib.pyplot as plt

# Get command line arguments
## There should be 2 arguments: the plaza id and the species name, if not, exit
if len(sys.argv) != 3:
    print("ERROR: wrong number of arguments")
    print("USAGE: python3 overlap_lncRNAs.py <plaza_id> <species_name>")
    exit(1)

plaza_id = sys.argv[1]
species_name = sys.argv[2]

# Load in lncRNA GFF file obtained from PLncDB (obtained on 2024-01-31)
lncRNA_bed = pybedtools.BedTool('plncdb/%s.lncRNA.gff'%(plaza_id))
## Filter for transcript features only (type: transcript)
lncRNA_bed_transcript = lncRNA_bed.filter(lambda feature: feature[2] == "transcript").saveas("tmp_%s.gff"%(plaza_id))
## Load in mRNA GFF file 
gene_bed = pybedtools.BedTool("gff/annotation.selected_transcript.mRNA.%s.gff3"%(species_name, plaza_id))
## intersect the lncRNA with the genes
results = lncRNA_bed_transcript.intersect(gene_bed, wo = True)
## keep only the lncRNAs that do not overlap with mRNA
lncRNA_bed_transcript = lncRNA_bed_transcript.intersect(gene_bed, v = True)

# Load in pseudogene BED file
pseudogene_bed = pybedtools.BedTool('%s/pseudogenes_%s.bed'%(plaza_id, plaza_id))

# get number of overlaps between (shuffled or real) lncRNA and pseudogenes (> 50 bp overlap)
def count_overlaps(results):
    count = 0
    count_wo_frag = 0
    count_dup = 0
    count_frag = 0
    count_wgm = 0
    count_pssd = 0
    for line in results:
        if (int(line[-1]) >= (0.5 * int(line[-2]))) and (int(line[-1]) > 30) and (int(line[-1]) > (0.5 * (int(line[4]) - int(line[3]) + 1))):
            count += 1
            # check if the pseudogene is fragmented
            if line[-4] != "Fragmented Ψ's":
                count_wo_frag += 1
            # check the type of pseudogene
            if line[-4] == "SSD-derived Ψ's":
                count_dup += 1
            elif line[-4] == "Fragmented Ψ's":
                count_frag += 1
            elif line[-4] == "WGM-derived Ψ's":
                count_wgm += 1
            elif line[-4] == "Retro-transposed Ψ's":
                count_pssd += 1
            else:
                print("ERROR")
    assert count == (count_dup + count_frag + count_wgm + count_pssd)
    return {"count": count,
            "count_wo_frag": count_wo_frag,
            "count_dup": count_dup, 
            "count_frag": count_frag, 
            "count_wgm": count_wgm, 
            "count_pssd": count_pssd}

# intersect the lncRNA with the pseudogenes
results = lncRNA_bed_transcript.intersect(pseudogene_bed, wo = True)
print(results.head())
real_count = count_overlaps(results)
print("real count: " + str(real_count.get("count")))
print("real count_wo_frag: " + str(real_count.get("count_wo_frag")))
print("real count_dup: " + str(real_count.get("count_dup")))
print("real count_frag: " + str(real_count.get("count_frag")))
print("real count_wgm: " + str(real_count.get("count_wgm")))
print("real count_pssd: " + str(real_count.get("count_pssd")))

# shuffle the lncRNA 1000 times and intersect with the pseudogenes
distribution_count = []
distribution_count_wo_frag = []
distribution_count_dup = []
distribution_count_frag = []
distribution_count_wgm = []
distribution_count_pssd = []
for i in range(0, 1000):
    print("iteration: " + str(i))
    # shuffle the lncRNA
    lncRNA_bed_shuffled = lncRNA_bed_transcript.shuffle(g= "%s/%s.chromsizes"%(plaza_id, plaza_id),
                                            excl = "gff/annotation.selected_transcript.mRNA.%s.gff3"%(species_name, plaza_id))

    # intersect the shuffled lncRNA with the pseudogenes
    results = lncRNA_bed_shuffled.intersect(pseudogene_bed, wo = True)
    print(results.head())
    count = count_overlaps(results)

    # append the count to the distribution
    distribution_count.append(count.get("count"))
    distribution_count_wo_frag.append(count.get("count_wo_frag"))
    distribution_count_dup.append(count.get("count_dup"))
    distribution_count_frag.append(count.get("count_frag"))
    distribution_count_wgm.append(count.get("count_wgm"))
    distribution_count_pssd.append(count.get("count_pssd"))

pval_count = sum([1 for i in distribution_count if i > real_count.get("count")])/len(distribution_count)
pval_count_wo_frag = sum([1 for i in distribution_count_wo_frag if i > real_count.get("count_wo_frag")])/len(distribution_count_wo_frag)
pval_count_dup = sum([1 for i in distribution_count_dup if i > real_count.get("count_dup")])/len(distribution_count_dup)
pval_count_frag = sum([1 for i in distribution_count_frag if i > real_count.get("count_frag")])/len(distribution_count_frag)
pval_count_wgm = sum([1 for i in distribution_count_wgm if i > real_count.get("count_wgm")])/len(distribution_count_wgm)
pval_count_pssd = sum([1 for i in distribution_count_pssd if i > real_count.get("count_pssd")])/len(distribution_count_pssd)
print("pval_count: " + str(pval_count))
print("pval_count_wo_frag: " + str(pval_count_wo_frag))
print("pval_count_dup: " + str(pval_count_dup))
print("pval_count_frag: " + str(pval_count_frag))
print("pval_count_wgm: " + str(pval_count_wgm))
print("pval_count_pssd: " + str(pval_count_pssd))

# Write p-values to file
with open("%s/results.txt"%(plaza_id), "w") as f:
    f.write("name" + "\t" + "pval" + "\t" + "real" + "\n")
    f.write("count" + "\t" + str(pval_count) + "\t" + str(real_count.get("count")) + "\n")
    f.write("count_wo_frag" + "\t" + str(pval_count_wo_frag) + "\t" + str(real_count.get("count_wo_frag")) + "\n")
    f.write("count_dup" + "\t" + str(pval_count_dup) + "\t" + str(real_count.get("count_dup")) + "\n")
    f.write("count_frag" + "\t" + str(pval_count_frag) + "\t" + str(real_count.get("count_frag")) + "\n")
    f.write("count_wgm" + "\t" + str(pval_count_wgm) + "\t" + str(real_count.get("count_wgm")) + "\n")
    f.write("count_pssd" + "\t" + str(pval_count_pssd) + "\t" + str(real_count.get("count_pssd")) + "\n")

# plot the distribution
def plot_distr(distribution_count, name):
    plt.hist(distribution_count)
    plt.axvline(x=real_count.get(name), color='red', linestyle='--', label='Vertical Line at x=1582')
    # save the plot
    plt.savefig("%s/%s_%s_overlap_distribution_plncdb.png"%(plaza_id, name, plaza_id))
    plt.close()

plot_distr(distribution_count, "count")
plot_distr(distribution_count_wo_frag, "count_wo_frag")
plot_distr(distribution_count_dup, "count_dup")
plot_distr(distribution_count_frag, "count_frag")
plot_distr(distribution_count_wgm, "count_wgm")
plot_distr(distribution_count_pssd, "count_pssd")

# Remove temporary file
os.remove("tmp_%s.gff"%(plaza_id))