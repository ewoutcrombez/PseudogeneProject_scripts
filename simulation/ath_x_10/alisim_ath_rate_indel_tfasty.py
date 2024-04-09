import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Initialize the result string
result_str = ""

split = sys.argv[1]

# Function to calculate percentages of two sequences
def calculate_percentages(seq1, seq2):
    identity = 0
    i = 0
    for x in range(len(seq1)):
        if (str(seq1[x]) != "-") & (str(seq2[x]) != "-"):
            i += 1
            if (seq1[x] == seq2[x]):
                identity += 1
    percent_identity = identity/i
    fraction = i/(len(seq1) - seq1.count("-"))
    length_aligned = i
    return [percent_identity, fraction, length_aligned]

# Load in the CDS sequences of Arabidopsis thaliana
cds_seq = SeqIO.to_dict(SeqIO.parse("../split_2/cds_ath_0MY.fasta.%s"%(str(split).zfill(3)), "fasta"))

# Loop over each CDS sequence and evolve it for 100 generations
for cds in cds_seq:
    # Translate the ancestral sequence and get length of protein
    prot_len = len(cds_seq[cds].translate().seq)

    if prot_len <= 5:
      print("Very small protein, ignore...")
      continue

    else:
        # Write the ancestral sequence to a file
        with open("tmp/tmp_0_%s.fa"%split, "w") as f:
            f.write(">A\n%s\n"%(cds_seq[cds].seq))
        
        # Write the ancestral protein sequence to a file
        with open("tmp/tmp_0_%s_prot.fa"%split, "w") as f:
            f.write(">A\n%s\n"%(cds_seq[cds].translate().seq))

        print("Working on %s"%cds_seq[cds].id, flush=True)

        # Run alisim for 100 generations
        for i in range(1, 101):
            j = i - 1
            command = ["iqtree2", "--alisim", "tmp/tmp_%s_%s"%(i, split), "alignSequences", "--seqtype", "DNA",
                                "-m", "GTR+F{0.32/0.18/0.18/0.32}", "--indel", "0.0429,0.0857",
                                "--indel-size", "LAV{0.5/3},LAV{0.5/3}",
                                "-t", "tree.nwk", "--root-seq",
                                "tmp/tmp_%s_%s.fa,A"%(j, split), "--out-format", "fasta"]
            result = subprocess.run(command, stdout=subprocess.PIPE)

            # Load in the evolved sequence of alisim
            evolved = SeqIO.to_dict(SeqIO.parse(f"tmp/tmp_{i}_{split}.unaligned.fa", "fasta"))
            # Change the header of the evolved sequence to "A_evolved"
            evolved["A"].id = "A_evolved"
            # Write the evolved sequence to a file
            with open("tmp/tmp_%s_%s.evolved.fa"%(i, split), "w") as f:
                f.write(">A_evolved\n%s\n"%(evolved["A"].seq))

            # Align ancestral and evolved CDS sequence using tfasty
            command = ["tfasty", "-m", "3" , "-Q", "-3",
                        "tmp/tmp_0_%s_prot.fa"%split, "tmp/tmp_%s_%s.evolved.fa"%(i, split),
                        "-O", "tmp/tmp_alignment_AA_%s_%s.txt"%(i, split)]        
            result = subprocess.run(command, stdout=subprocess.PIPE)

            # Parse the results of tfasty
            ancestor = False
            ancestor_prot_seq = ""
            evolved = False
            evolved_prot_seq = ""
            stats = {}
            ancestor_seqs = {}
            evolved_seqs = {}
            match_nr = 1
            match = "match" + str(match_nr)
            for line in open("tmp/tmp_alignment_AA_%s_%s.txt"%(i, split)):
                # Get statistics of the alignment
                if line.startswith("A_evolved"):
                    e_val = line.strip().split(" ")[-1]
                if line.startswith("Smith-Waterman score"):
                    scores = line.strip().split(" ")
                    scores = [x for x in scores if x != ""]
                    print(scores, flush=True)
                    ident = float(scores[3].replace("%", "")) / 100
                    similarity = float(scores[5].replace("(", "").replace("%", "")) / 100
                    overlap = int(scores[8])
                    start_overlap = int(scores[11].split(":")[1].split("-")[0])
                    end_overlap = int(scores[11].split(":")[1].split("-")[1].replace(")", ""))
                    stats[match] = [ident, similarity, overlap, start_overlap, end_overlap]
                # Get the ancestral aligned protein sequence
                if line.startswith(">A"):
                    ancestor = True
                # Get the evolved aligned protein sequence
                if line.startswith(">A_evol"):
                    ancestor = False
                    evolved = True
                if line.startswith("\n"):
                    ancestor = False
                    evolved = False
                # New match starts, initialize again
                if line.startswith(">--"):
                    ancestor = False
                    evolved = False
                    ancestor_prot_seq = ancestor_prot_seq.replace(">A ..", "")
                    ancestor_seqs[match] = ancestor_prot_seq
                    evolved_prot_seq = evolved_prot_seq.replace(">A_evol ..", "")
                    evolved_seqs[match] = evolved_prot_seq
                    ancestor_prot_seq = ""
                    evolved_prot_seq = ""
                    match_nr += 1
                    match = "match" + str(match_nr)
                if ancestor:
                    ancestor_prot_seq += line.strip()
                if evolved:
                    evolved_prot_seq += line.strip()
            # Save the last match
            ancestor_prot_seq = ancestor_prot_seq.replace(">A ..", "")
            ancestor_seqs[match] = ancestor_prot_seq
            evolved_prot_seq = evolved_prot_seq.replace(">A_evol ..", "")
            evolved_seqs[match] = evolved_prot_seq
        
            # Sort the different matches by their start position
            print(stats)
            stats = dict(sorted(stats.items(), key = lambda x: x[1][3]))
            ancestor_seqs = dict(sorted(ancestor_seqs.items(), key = lambda x: list(stats.keys()).index(x[0])))
            evolved_seqs = dict(sorted(evolved_seqs.items(), key = lambda x: list(stats.keys()).index(x[0])))

            # If there are matches that overlap, pick the longest match
            if len(stats) > 1:
                for match in list(stats.keys()):
                    print(match, flush=True)
                    # if match is not the first index
                    if list(stats.keys()).index(match) > 0:
                        # if the start position of the match is smaller than the end position of the previous match
                        if stats[match][3] < stats[list(stats.keys())[list(stats.keys()).index(match) - 1]][4]:
                            # if the length of the match is smaller than the length of the previous match
                            if stats[match][2] < stats[list(stats.keys())[list(stats.keys()).index(match) - 1]][2]:
                                # remove the match
                                print("Removed match %s"%match, flush=True)
                                del stats[match]
                                del ancestor_seqs[match]
                                del evolved_seqs[match]
                            
                            else:
                                # remove the previous match
                                print("Removed match %s"%list(stats.keys())[list(stats.keys()).index(match) - 1], flush=True)
                                del stats[list(stats.keys())[list(stats.keys()).index(match) - 1]]
                                del ancestor_seqs[list(stats.keys())[list(stats.keys()).index(match) - 1]]
                                del evolved_seqs[list(stats.keys())[list(stats.keys()).index(match) - 1]]

            # Overlap is the sum of the overlap of the different matches
            overlap = sum([stats[match][2] for match in stats])
            # Calculate the weighted average statistics of the different matches
            ident = sum((stats[match][2] / overlap) * stats[match][0] for match in stats)
            similarity = sum((stats[match][2] / overlap) * stats[match][1] for match in stats)

            # Combine the different matches into one sequence
            ancestor_prot_seq = "".join(list(ancestor_seqs.values()))
            evolved_prot_seq = "".join(list(evolved_seqs.values()))

            # Convert the aligned ancestral and evolved protein sequence to a Biopython sequence object
            ancestor_prot_seq = Seq(ancestor_prot_seq)
            evolved_prot_seq = Seq(evolved_prot_seq)

            # Calculate percent identity of ancestral and evolved protein
            result_prot = calculate_percentages(ancestor_prot_seq, evolved_prot_seq)
            percent_identity_prot = result_prot[0]
            fraction_prot = result_prot[1]
            length_aligned_prot = result_prot[2]

            # Get number of frameshifts
            num_fs = evolved_prot_seq.count("\\")
            num_fs += evolved_prot_seq.count("/")
            # Get number of stop codons
            num_stop = evolved_prot_seq.count("*")
            
            # Print the results
            result_str = cds + "\t" + str(prot_len) + "\t" + str(i) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(num_fs) + "\t" + str(num_stop) + "\t" + str(ident) + "\t" + str(similarity) +  "\t" + str(overlap) + "\t" + str(length_aligned_prot / prot_len) + "\t" + e_val + "\n"
            print(result_str, flush=True)
            # Write the result to a file
            with open("result_sim_%s.tsv"%split, "a") as f:
                f.write(result_str)
