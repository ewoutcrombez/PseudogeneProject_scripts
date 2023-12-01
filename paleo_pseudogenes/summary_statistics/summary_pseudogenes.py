import sys
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

# Initialize the result string
result_str = ""

# Function to calculate percent identity and fraction
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

# Input: 1) pseudogenes file 2) gene file
species = sys.argv[1]
plaza_name = sys.argv[2]
split = sys.argv[3]

# pseudogenes file
pseudogene_exon_file = "/scratch/recent_wgds/paleo_polyploids/filtering/exons/split/exons_%s.gff.%s"%(species, str(split))
# gene file
genes = "/scratch/recent_wgds/data/%s/fasta/cds.selected_transcript.%s.fasta"%(species, plaza_name)
genome = "/scratch/recent_wgds/data/%s/fasta/%s.fasta"%(species, plaza_name)

# Read genome file
genome = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

# Read gene file
genes = SeqIO.to_dict(SeqIO.parse(genes, "fasta"))

pseudogene_seq_full = []
# Get pseudogene sequences based on coordinates
for pseudogene_exons in open(pseudogene_exon_file):
    # Get pseudogene information (name, chromosome, start, end, strand, functional paralogue gene)
    pseudogene_exon = pseudogene_exons.split("\t")
    print(pseudogene_exon, flush=True)
    pseudogene_chr_exon = pseudogene_exon[0]
    pseudogene_start_exon = int(pseudogene_exon[1])
    pseudogene_end_exon = int(pseudogene_exon[2])
    pseudogene_strand_exon = pseudogene_exon[3]
    gene = pseudogene_exon[4]
    pseudogene = pseudogene_exon[22].strip()
    exon_num = int(pseudogene_exon[6])
    total_exons = int(pseudogene_exon[7])
    print(pseudogene, flush=True)
    # Get pseudogene sequence
    pseudogene_seq = genome[pseudogene_chr_exon][pseudogene_start_exon-1:pseudogene_end_exon]
    print(pseudogene_seq, flush=True)
    if pseudogene_strand_exon == "-":
        pseudogene_seq = pseudogene_seq.reverse_complement()

    if exon_num < total_exons:
        pseudogene_seq_full.append(str(pseudogene_seq.seq))
    elif exon_num == total_exons:
        print("Last exon", flush=True)
        # Get gene sequence
        gene_seq = genes[gene]
        if pseudogene_strand_exon == "-":
            pseudogene_seq_full.append(str(pseudogene_seq.seq))
            pseudogene_seq_full = "".join(pseudogene_seq_full[::-1])
            pseudogene_seq_full = Seq(pseudogene_seq_full)
            print(pseudogene_seq_full, flush=True)
        else:
            pseudogene_seq_full.append(str(pseudogene_seq.seq))
            pseudogene_seq_full = "".join(pseudogene_seq_full)
            pseudogene_seq_full = Seq(pseudogene_seq_full)
        print(pseudogene_seq_full, flush=True)
        # Write pseudogene and gene sequences to file
        with open("pseudogene_%s_%s.fa"%(species, split), "w") as f:
            f.write(">%s\n%s\n"%(pseudogene, pseudogene_seq_full))
        with open("gene_%s_%s.fa"%(species, split), "w") as f:
            f.write(">%s\n%s\n"%(gene, gene_seq.seq))
        print(pseudogene_seq_full, flush=True)
        pseudogene_seq_full = []

        # Align pseudogene and gene sequences with MACSE
        command = ["java", "-jar", "/home/ewcro/tools/macse_v2.06.jar", "-prog", "alignSequences", "-seq",
                        "gene_%s_%s.fa"%(species, split), "-seq_lr", "pseudogene_%s_%s.fa"%(species, split), "-out_NT", 
                        "alignment_NT_%s_%s.fa"%(species, split), "-out_AA", 
                        "alignment_AA_%s_%s.fa"%(species, split), "-fs", "10000", "-fs_lr", "10", "-stop_lr", "15"]
        result = subprocess.run(command, stdout=subprocess.PIPE)

        # Load in the aligned pseudogene and gene sequences
        aligned_cds = SeqIO.to_dict(SeqIO.parse("alignment_NT_%s_%s.fa"%(species, split), "fasta"))
        # Extract the aligned pseudogene and gene sequences
        # If there is no alignment, set percent identity and fraction to NA
        if pseudogene not in aligned_cds.keys():
            print("Error: No alignment", flush=True)
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
        else:
            gene_cds = aligned_cds[gene].seq
            pseudogene_cds = aligned_cds[pseudogene].seq

            # Calculate percent identity of pseudogene and gene DNA
            result_cds = calculate_percentages(gene_cds, pseudogene_cds)
            percent_identity_cds = result_cds[0]
            fraction_cds = result_cds[1]
            length_aligned_cds = result_cds[2]

            # Load in the aligned pseudogene and gene protein sequences
            aligned_prot = SeqIO.to_dict(SeqIO.parse("alignment_AA_%s_%s.fa"%(species, split), "fasta"))

            # Extract the aligned pseudogene and gene protein sequences
            gene_prot = aligned_prot[gene].seq
            print(gene_prot, flush=True)
            pseudogene_prot = aligned_prot[pseudogene].seq
            print(pseudogene_prot, flush=True)
            # Get number of frameshifts
            num_fs = pseudogene_prot.count("!")
            # Get number of stop codons
            num_stop = pseudogene_prot.count("*")

            # Calculate percent identity of pseudogene and gene protein
            result_prot = calculate_percentages(gene_prot, pseudogene_prot)
            percent_identity_prot = result_prot[0]
            fraction_prot = result_prot[1]
            length_aligned_prot = result_prot[2]

            # Change fasta to axt
            command = ["fasta_to_axt.pl", "alignment_NT_%s_%s.fa"%(species, split), "alignment_%s_%s"%(species, split)]
            result = subprocess.run(command, stdout=subprocess.PIPE)
            
            # Run KaKs_Calculator
            command = ["KaKs_Calculator", "-i", "alignment_%s_%s.axt"%(species, split), "-o", "alignment_%s_%s.kaks"%(species, split), "-m", "MA"]
            result = subprocess.run(command, stdout=subprocess.PIPE)

            # Read in KaKs_Calculator results
            i = 1
            for line in open("alignment_%s_%s.kaks"%(species, split)):
                if i == 1:
                    i += 1
                    continue
                else:
                    line = line.split("\t")
                    Ka = line[2]
                    Ks = line[3]
                    Ka_Ks = line[4]
                    p_value = line[5]

        # Print results of percent identity, fraction, Ka, Ks, Ka/Ks, and p-value
        result_str += pseudogene + "\t" + gene + "\t" + str(percent_identity_cds) + "\t" + str(fraction_cds) + "\t" + str(length_aligned_cds) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(Ka) + "\t" + str(Ks) + "\t" + str(Ka_Ks) + "\t" + str(p_value) + "\t" + str(num_fs) + "\t" + str(num_stop) + "\n"
        print(pseudogene + "\t" + gene + "\t" + str(percent_identity_cds) + "\t" + str(fraction_cds) + "\t" + str(length_aligned_cds) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(Ka) + "\t" + str(Ks) + "\t" + str(Ka_Ks) + "\t" + str(p_value) + "\t" + str(num_fs) + "\t" + str(num_stop) + "\n", flush=True)
    
    
    else:
        print("Error: exon_num > total_exons", flush=True)
    

# Write results to file
with open("result_%s_%s.tsv"%(species, split), "w") as f:
    f.write(result_str)
