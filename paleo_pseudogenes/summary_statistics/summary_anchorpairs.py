import sys
import os
import subprocess
from Bio import SeqIO

# Input: cds file
species = sys.argv[1]
plaza_name = sys.argv[2]
split = sys.argv[3]

cds = "/scratch/recent_wgds/data/%s/fasta/cds.selected_transcript.%s.fasta"%(species, plaza_name)

# Get anchor points
anchor_points = "%s/anchorpoints.txt.%s"%(species, str(split).zfill(3))

# Load in the CDS sequences
cds_seq = SeqIO.to_dict(SeqIO.parse(cds, "fasta"))

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
    fraction = i/max((len(seq1) - seq1.count("-")), (len(seq2) - seq2.count("-")))
    length_aligned = i
    return [percent_identity, fraction, length_aligned]

result_str = ""
# Loop over each anchor pair
for line in open(anchor_points):
    if line.strip().startswith('id'):
        continue  # Skip this line and move to the next one
    line = line.strip().split("\t")
    # Get the anchor points
    anchor_1 = line[3]
    anchor_2 = line[4]
    # Get the anchor point sequence
    anchor_1_seq = cds_seq[anchor_1].seq
    anchor_2_seq = cds_seq[anchor_2].seq
    # Write the anchor point sequence to a file
    with open("tmp/tmp_anchor_%s_%s.fa"%(species, split), "w") as f:
        f.write(">%s\n%s\n"%(anchor_1, anchor_1_seq))
    with open("tmp/tmp_anchor_%s_%s.fa"%(species, split), "a") as f:
        f.write(">%s\n%s\n"%(anchor_2, anchor_2_seq))
    
    # Align anchor points using MACSE
    command = ["java", "-jar", "/home/ewcro/tools/macse_v2.06.jar", "-prog", "alignSequences", "-seq",
                    "tmp/tmp_anchor_%s_%s.fa"%(species, split), "-out_NT", 
                    "tmp/tmp_alignment_NT_%s_%s.fa"%(species, split), "-out_AA", 
                    "tmp/tmp_alignment_AA_%s_%s.fa"%(species, split)]
    result = subprocess.run(command, stdout=subprocess.PIPE)

    # Load in the aligned anchorpoint CDS sequences
    aligned_cds = SeqIO.to_dict(SeqIO.parse("tmp/tmp_alignment_NT_%s_%s.fa"%(species, split), "fasta"))

    if anchor_1 not in aligned_cds.keys() or anchor_2 not in aligned_cds.keys():
        percent_identity_cds = "NA"
        fraction_cds = "NA"
        length_aligned_cds = "NA"
        percent_identity_prot = "NA"
        fraction_prot = "NA"
        length_aligned_prot = "NA"
        Ka = "NA"
        Ks = "NA"
        Ka_Ks = "NA"
        p_value = "NA"
        continue
    else:
        # Extract the aligned CDS sequences
        ap1_cds = aligned_cds[anchor_1].seq
        ap2_cds = aligned_cds[anchor_2].seq
        
        # Calculate percent identity of ancestral and evolved DNA
        results_cds = calculate_percentages(ap1_cds, ap2_cds)
        percent_identity_cds = results_cds[0]
        fraction_cds = results_cds[1]
        length_aligned_cds = results_cds[2]

        # Translate the anchorpoint CDS sequences and write it to a file
        prot_len1 = len(anchor_1_seq.translate())
        prot_len2 = len(anchor_2_seq.translate())

        # Load in the aligned anchorpoint protein sequences
        aligned_prot = SeqIO.to_dict(SeqIO.parse("tmp/tmp_alignment_AA_%s_%s.fa"%(species, split), "fasta"))
        # Extract the aligned protein sequences
        ap1_prot = aligned_prot[anchor_1].seq
        ap2_prot = aligned_prot[anchor_2].seq

        # Calculate percentages of anchorpoint proteins
        results_prot = calculate_percentages(ap1_prot, ap2_prot)
        percent_identity_prot = results_prot[0]
        fraction_prot = results_prot[1]
        length_aligned_prot = results_prot[2]

        # Change fasta to axt
        command = ["fasta_to_axt.pl", "tmp/tmp_alignment_NT_%s_%s.fa"%(species, split), "tmp/tmp_alignment_%s_%s"%(species, split)]
        result = subprocess.run(command, stdout=subprocess.PIPE)

        # Run KaKs_Calculator
        command = ["KaKs_Calculator", "-i", "tmp/tmp_alignment_%s_%s.axt"%(species, split), "-o", "tmp/tmp_alignment_%s_%s.kaks"%(species, split), "-m", "MA"]
        result = subprocess.run(command, stdout=subprocess.PIPE)

        # Read in KaKs_Calculator results
        i = 1
        for line in open("tmp/tmp_alignment_%s_%s.kaks"%(species, split)):
            if i == 1:
                i += 1
                continue
            else:
                line = line.split("\t")
                Ka = line[2]
                Ks = line[3]
                Ka_Ks = line[4]
                p_value = line[5]

    # Write results to a string
    # Print the results
    result_str += anchor_1 + "\t" + str(prot_len1) + "\t" + anchor_2 + "\t" + str(prot_len2) + "\t" + str(percent_identity_cds) + "\t" + str(fraction_cds) + "\t" + str(length_aligned_cds) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(Ka) + "\t" + str(Ks) + "\t" + str(Ka_Ks) + "\t" + str(p_value) + "\n"
    print(anchor_1 + "\t" + str(prot_len1) + "\t" + anchor_2 + "\t" + str(prot_len2) + "\t" + str(percent_identity_cds) + "\t" + str(fraction_cds) + "\t" + str(length_aligned_cds) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(Ka) + "\t" + str(Ks) + "\t" + str(Ka_Ks) + "\t" + str(p_value) + "\n", flush=True)

# Write the results to a file
with open("%s/anchorpoint_result_%s_%s.txt"%(species, species, split), "w") as f:
    f.write(result_str)