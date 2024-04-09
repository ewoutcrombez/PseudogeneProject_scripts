from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import re

# specify input and output filenames
input_file = sys.argv[1]
output_file = sys.argv[2]

# initialize
full_seq = Seq("")
prev_parent_hit_nr = ""
seq_list = list()
# read in the input fasta file and combine sequences if exons from same gene
for record in SeqIO.parse(input_file, "fasta"):
    # get info about the sequence
    list_descr = str(record.description).split("_")
    parent_hit_nr = "_".join(list_descr[0:2])
    seq = record.seq
    # if first sequence, then set full_seq to the sequence
    if full_seq == Seq(""):
        full_seq = seq
        prev_parent_hit_nr = parent_hit_nr
    # if the current sequence is from same hit as previous one, combine the sequences together
    elif parent_hit_nr == prev_parent_hit_nr:
        full_seq += seq
    # if the current sequence is new, append the previous one to list
    else:
        seq_list.append(SeqRecord(full_seq, id = prev_parent_hit_nr, description = ""))
        # re-initialize
        full_seq = seq
        prev_parent_hit_nr = parent_hit_nr
# append last sequence
seq_list.append(SeqRecord(full_seq, id = prev_parent_hit_nr, description = ""))

# If seq list only contains an empty sequence
if len(seq_list) == 1 and str(seq_list[0].seq) == "":
    print(" No hits")
else:
    # write the combined records to the output fasta file
    with open(output_file, "w") as out_handle:
        SeqIO.write(seq_list, out_handle, "fasta")

