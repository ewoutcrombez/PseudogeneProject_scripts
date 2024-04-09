# In this script I will process the exonerate result file and split it up in a CSV and GFF file.
import sys

# initialize variables
gff_result = ""
tsv_result = ""

# input from command line
exonerate_result_file = sys.argv[1]
gff_file = sys.argv[2]
tsv_file = sys.argv[3]

with(open(exonerate_result_file)) as f:
    # iterate over lines
    for line in f:
        # check if line is a comment
        if line.startswith("#"):
            # check if line is the start of the GFF dump
            if line.startswith("# --- START OF GFF DUMP ---"):
                gff = True # set gff to True
                tsv = False # set csv to False
            # check if line is the end of the GFF dump
            elif line.startswith("# --- END OF GFF DUMP ---"):
                gff = False # set gff to False
                tsv = True # set csv to True
            else:   
                continue
        else:
            # write line to GFF or CSV file
            if gff:
                gff_result += line
            elif tsv:
                tsv_result += line

# write GFF file
with(open(gff_file, "w")) as f:
    f.write(gff_result)
# write CSV file
with(open(tsv_file, "w")) as f:
    f.write(tsv_result)
