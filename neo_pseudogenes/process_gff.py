# In this script, the gff file is processed to merge hits that are within 6000 bp of each other and 
# a name is given to each hit and exon
import sys
import re

# input from command line
gff_file = sys.argv[1]
gff_select_file = sys.argv[2] # file with the selected hits
out_file_exon = sys.argv[3]

# create list with the gff lines of the selected hits from `select_best_hit_gene_region.R`
selected_list = []
with open(gff_select_file) as f:
    for line in f:
        selected_hits = "\t".join(line.split("\t")[0:9]) + "\n"
        selected_list.append(selected_hits)
# create a dictionary that has the hit name as key and all exon lines of that hit as value in a list
dictio = {}
hits = []
hit_number = 0
with open(gff_file) as f: # open file
    for line in f: # iterate over lines
        if line.startswith("#"): # skip comment lines
            continue
        fields = line.split("\t")
        feature_type = fields[2] # get feature type
        if feature_type == "gene":
            if line in selected_list: # if the line is in the selected list
                selected = True
                hit_id = ""
                hit_number += 1
                target_id = fields[0]
                start = int(fields[3]) # get start position 
                end = int(fields[4]) # get end position
                key = "".join(fields[8].split(";")[0:2]).strip() + " " + str(start) + "-" + str(end) # create a key from the gene name and the start and end position
                query = fields[8].split(";")[1].strip().split(" ")[1] # get query name
                orientation = fields[6] # get orientation
                # if there is already a hit that overlaps with the start and end position, give it the same name
                for x in hits:
                    # check if the start and end position overlap
                    if (start >= x['start'] and start <= x['end'] and orientation == x['orientation']) or (end >= x['start'] and end <= x['end'] and orientation == x['orientation']):
                        # get the hit name of the hit that overlaps with the start and end position,
                        # remove the query name from the hit name and add the current query name
                        hit_id = query + re.sub(x['query'], '', x['hit_id'])
                        x['start'] = min(start, x['start']) # get the minimum start position
                        x['end'] = max(end, x['end']) # get the maximum end position
                # if there is no hit that overlaps, create a new hit name
                if hit_id == "":
                    hit_id = query + "_hit" + str(hit_number)
                    hit = { # create a dictionary with the hit name, start, end and query name
                    'query': query,
                    'start': start,
                    'end': end,
                    'orientation': orientation,
                    'hit_id': hit_id,
                    'target_id': target_id,
                    'key': key
                    }
                    hits.append(hit)
            else: # if the line is not in the selected list, skip it
                selected = False
        elif feature_type == "exon":
            if selected: # if the gene line was in the selected list
                if key in dictio.keys(): # if the hit name is already in the dictionary, append the line to the list
                    i += 1
                    dictio[key].append(line[:-2] + "\t" + hit_id + "_exon" + str(i) + "\n")

                else: # if the hit name is not in the dictionary, create a new key and add the line to the list
                    i = 1
                    dictio[key] = list()
                    dictio[key].append(line[:-2] + "\t" + hit_id + "_exon" + str(i) + "\n")
# sort hits by start position
hits = sorted(hits, key=lambda k: k['start'])

# merge hits that are within 6000 bp of each other
merge_range = 6000
merged_hits = dict()
merged_num = 0
for i in range(len(hits)): # iterate over all hits
    for j in range(i+1, len(hits)): # iterate over all hits after the current hit
        if hits[i]['query'] != hits[j]['query']: # if the query names are not the same
            continue
        # if the hits are within 6000 bp of each other, merge them
        if (abs(hits[i]['start'] - hits[j]['end']) <= merge_range or abs(hits[j]['start'] - hits[i]['end']) <= merge_range) and (hits[i]['orientation'] == hits[j]['orientation']):
            merged_num += 1
            if (hits[i]['key'] in merged_hits.keys()) or (hits[j]['key'] in merged_hits.keys()): # if one of the hits is already merged
                already_merged = hits[i]['key'] if hits[i]['key'] in merged_hits.keys() else hits[j]['key'] # get the already merged hit

                merged_start = min(hits[i]['start'], hits[j]['start'], merged_hits[already_merged]['start']) # get the start position of the merged hit
                merged_end = max(hits[i]['end'], hits[j]['end'], merged_hits[already_merged]['end']) # get the end position of the merged hit
                merged_hit = { # create a dictionary with the merged hit name, start, end and query name
                    'query': hits[i]['query'],
                    'start': merged_start,
                    'end': merged_end,
                    'orientation': hits[i]['orientation'],
                    'hit_id': hits[i]['hit_id'],
                    'target_id': hits[i]['target_id'],
                    'key': "gene_id merged " + hits[i]['query'] + " " + str(merged_start) + "-" + str(merged_end),
                    # old keys of the hits that are merged
                    'old_keys': [merged_hits[already_merged]['key'], hits[j]['key'] if hits[i]['key'] in merged_hits.keys() else hits[i]['key']]
                }
                merged_hits[merged_hits[already_merged]['key']] = merged_hit # add the merged hit that is merged further to the merged hits dictionary as key
                merged_hits[hits[j]['key']] = merged_hit # add the other hit that is merged with the already merged hit to the merged hits dictionary as key
                # Obtain lines and add new hit and exon names
                merged_line = ["\t".join(exon[1].split("\t")[0:9]) + "\t" + exon[1].split("\t")[9].split("_")[0] + "_hitmerged" + str(merged_num) + "_exon" + str(exon[0] + 1) + "\n" for exon in enumerate(dictio[merged_hit['old_keys'][0]] + dictio[merged_hit['old_keys'][1]])]
                dictio[merged_hit['key']] = merged_line # add the merged line to the dictionary
            else: # if none of the hits is already merged
                merged_start = min(hits[i]['start'], hits[j]['start'])
                merged_end = max(hits[i]['end'], hits[j]['end'])
                merged_hit = { # create a dictionary with the merged hit name, start, end and query name
                    'query': hits[i]['query'],
                    'start': merged_start,
                    'end': merged_end,
                    'orientation': hits[i]['orientation'],
                    'hit_id': hits[i]['hit_id'],
                    'target_id': hits[i]['target_id'],
                    'key': "gene_id merged " + hits[i]['query'] + " " + str(merged_start) + "-" + str(merged_end),
                    # old keys of the hits that are merged
                    'old_keys': [hits[i]['key'], hits[j]['key']]
                }
                # Obtain lines and add new hit and exon names
                merged_lines = ["\t".join(exon[1].split("\t")[0:9]) + "\t" + exon[1].split("\t")[9].split("_")[0] + "_hitmerged" + str(merged_num) + "_exon" + str(exon[0] + 1) + "\n" for exon in enumerate(dictio[merged_hit['old_keys'][0]] + dictio[merged_hit['old_keys'][1]])]
                dictio[merged_hit['key']] = merged_lines # add the merged line to the dictionary

                merged_hit['hit_id'] = "_".join(merged_lines[0].split("\t")[9].split("_")[0:2])
                merged_hits[hits[i]['key']] = merged_hit # add the first hit that is merged to the merged hits dictionary as key
                merged_hits[hits[j]['key']] = merged_hit # add the second hit that is merged to the merged hits dictionary as key
for old_key in merged_hits.keys(): # iterate over all old keys
    del dictio[old_key] # remove the old keys from the dictionary

# Combine overlapping exons
lines = []
previous_line = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0"
for key in dictio.keys():
        for line in dictio[key]:
            # if there is a line with overlap with the previous line, remove the previous line
            # and use the line with the start and end positions of the merged hit
            start = int(line.split("\t")[3])
            end = int(line.split("\t")[4])
            orientation = line.split("\t")[6]
            start_previous = int(previous_line.split("\t")[3])
            end_previous = int(previous_line.split("\t")[4])
            orientation_previous = previous_line.split("\t")[6]
            if (start >= start_previous and start <= end_previous and orientation == orientation_previous) or (end >= start_previous and end <= end_previous and orientation == orientation_previous):
                lines.pop() # remove the previous line
                lines.append("\t".join(line.split("\t")[0:3]) + "\t" + str(min(start, start_previous)) + "\t" + str(max(end, end_previous)) + "\t" + "\t".join(line.split("\t")[5:]))
            else:
                lines.append(line)
            previous_line = line

# print all lines of the dictionary containing the exon lines to a new file
with open(out_file_exon, "w") as f:
    for line in lines:
        # if line is not empty, write it to the file
        if line != "\n":
            f.write(line)