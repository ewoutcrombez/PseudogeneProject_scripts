import sys

def converttoDnaSeg(bed_file = sys.argv[1], folder = sys.argv[2], species = sys.argv[3], lonely_genes = sys.argv[4],
                    snpeff_result = sys.argv[5], annotation_psg = sys.argv[6]):
    string = 'name\tstart\tend\tstrand\tcol\tfill\n'
    lonely = 'lonely_gene\n'
    coordinate1 = 0
    coordinate2 = 150
    hits = []
    with open(lonely_genes) as f:
        lonelygenes = f.read().splitlines()

    with open(snpeff_result) as f:
        snpeff = f.read().splitlines()
    
    annot = {}
    linked_genes = []
    with open(annotation_psg) as f:
        i = 1
        for line in f:
            if i == 1:
                i += 1
                continue
            else:
                fields = line.split("\t")
                annot[fields[1]] = fields[-1].rstrip()
                linked_genes.append(fields[2])
        
    with open(bed_file) as f:
        for line in f:
            if line == '\n':
                next
            else:
                fields = line.split("\t")
                if "hit" in fields[3]:
                    # if fields[3] in annot.keys():
                    #     if annot[fields[3]] == "intact_unannotated_gene" or annot[fields[3]] == "likely_unannotated_gene":
                    #         col = '#686f7a' # grey
                    #         fill = '#686f7a' # grey
                    #     elif annot[fields[3]] == "intact_pseudogene":
                    #         col = 'blue' # blue
                    #         fill = 'blue' # blue
                    #     elif annot[fields[3]] == "slightly_fragmented_pseudogene" or annot[fields[3]] == "fragmented_pseudogene":
                    #         col = 'cyan'
                    #         fill = 'cyan'
                    #     elif annot[fields[3]] == "dubious_pseudogene":
                    #         col = 'yellow'
                    #         fill = 'yellow'
                    #     elif annot[fields[3]] == "partial_pseudogene":
                    #         col = 'red'
                    #         fill = 'red'
                    # else:
                    col = 'blue'
                    fill = 'blue'
                elif fields[3] in lonelygenes and not fields[3] in linked_genes:
                    col = 'black'
                    fill = 'black'
                    lonely += fields[3] + "\n"
                elif fields[3] in linked_genes:
                    col = 'green'
                    fill = 'green'
                elif fields[3] in snpeff:
                    col = 'purple'
                    fill = 'purple'
                else:
                    col = 'grey'
                    fill = 'grey'
                # if the hit is not in the list yet, add it to the list and add the hit to the string
                # this is to prevent the same hit to be added multiple times because they overlapped
                if fields[3].rstrip() not in hits:
                    hits.append(fields[3].rstrip()) # add hit to list
                    string += fields[3].rstrip() + '\t' + str(coordinate1) + '\t' + str(coordinate2) + '\t' + str(fields[4].rstrip()) + '\t' + col + '\t' + fill + '\n'
                    coordinate1 += 400
                    coordinate2 += 400
    with open("%s/%s.tab.txt"%(folder, species), "w") as f:
        f.write(string)

    with open("%s/%s_lonelygenes.txt"%(folder, species), "w") as f:
        f.write(lonely)

converttoDnaSeg()