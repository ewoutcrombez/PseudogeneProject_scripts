from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import re
import sys

working_dir = sys.argv[1]
multiplicon = sys.argv[2]
haplotype1 = sys.argv[3]
haplotype2 = sys.argv[4]
chrom = sys.argv[5]

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

# Load in exon result table
result_table_exons = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/hits_{haplotype1}_in_{haplotype2}_{chrom}_exons.gff", "r")

# Load in genome fasta files
haplotype_1 = SeqIO.to_dict(SeqIO.parse(f"{working_dir}/data/genome/{haplotype1}_haplotype_genome.fasta", "fasta"))
haplotype_2 = SeqIO.to_dict(SeqIO.parse(f"{working_dir}/data/genome/{haplotype2}_haplotype_genome.fasta", "fasta"))

# Load in CDS fasta files
haplotype_1_CDS = SeqIO.to_dict(SeqIO.parse(f"{working_dir}/data/CDS/{haplotype1}.CDS.fa", "fasta"))
haplotype_2_CDS = SeqIO.to_dict(SeqIO.parse(f"{working_dir}/data/CDS/{haplotype2}.CDS.fa", "fasta"))

# Header for result table
string = "hit_name\tpseudogene_len\tparent_gene\tgene_len\tpercent_identity_cds\tfraction_cds\tlength_aligned_cds\tpercent_identity_prot\tfraction_prot\tlength_aligned_prot\tnum_fs\tnum_stop\tstart_codon\tstop_codon\n"
"hit_name\thit_synonym\tparent_gene\tgene_length\tpseudogene_length\tinternal_FS_psg\tinternal_stop_psg\tinternal_del_psg\tgaps_in_psg_aln\tFS_in_psg\tmatches\tmismatches\tlongest_aligned_stretch\tlongest_gap_in_pseudogene\tstart_present\tstop_present\tinternal_FS_g\tinternal_stop_g\tinternal_del_g\tgaps_in_g_aln\tFS_in_g\tlongest_gap_in_gene\tannotation\n"
# Write to file
output = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/result_table_{haplotype1}_MACSE.tsv", "w")
output.write(string)
output.close()

# Read in the gff table
exon_coordinates = []
i = 1
for line in result_table_exons:
    # initialize the table
    if i == 1:
        line = line.strip().split("\t")
        hit = line[9].split("_exon")[0] # hit name
        hit_start = int(line[3]) # hit start coordinate
        hit_end = int(line[4]) # hit end coordinate
        strand = line[6] # hit strand
        chrom = line[0] # hit chromosome
        gene = line[9].split("_hit")[0] # corresponding gene match
        hit_syn = chrom + "_" + hit.split("_")[-1]
        previous_hit = hit # initialize previous hit
        previous_gene = gene # initialize previous gene
        previous_chrom = chrom # initialize previous chromosome
        previous_strand = strand # initialize previous strand
        previous_hit_syn = hit_syn # initialize previous hit synonym
        exon_coordinates.append((hit_start, hit_end, chrom, strand)) # add coordinates to list
        i += 1

    # read in the rest of the table
    else:
        line = line.strip().split("\t")
        hit = line[9].split("_exon")[0] # hit name
        hit_start = int(line[3]) # hit start coordinate
        hit_end = int(line[4]) # hit end coordinate
        strand = line[6] # hit strand
        chrom = line[0] # hit chromosome
        gene = line[9].split("_hit")[0] # corresponding gene match
        hit_syn = chrom + "_" + hit.split("_")[-1]

        # if the exon belongs to the same hit as the previous exon
        # add the coordinates to the list
        if hit == previous_hit:
            exon_coordinates.append((hit_start, hit_end, chrom, strand))
            previous_hit = hit
            previous_gene = gene
            previous_chrom = chrom
            previous_strand = strand
            previous_hit_syn = hit_syn

        # if the exon belongs to a new hit
        # align the gene to the exon coordinates
        # of the previous hit
        else:            
            # Obtain gene sequence
            if previous_gene.startswith(haplotype1):
                gene_seq = haplotype_1_CDS[previous_gene].seq
            else:
                gene_seq = haplotype_2_CDS[previous_gene].seq
            
            # save gene sequence to fasta file
            gene_seq_fasta = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq.fasta", "w")
            gene_seq_fasta.write(">gene_seq\n" + str(gene_seq))
            gene_seq_fasta.close()

            gene_len = len(gene_seq)

            # Obtain pseudogene sequence
            psg_seqs = []
            if previous_chrom.startswith(haplotype1):
                previous_chrom = previous_chrom.split("_")[1]
                # Get exon sequences
                for exon in exon_coordinates:
                    psg_seq = haplotype_1[previous_chrom].seq[exon[0]-1:exon[1]]
                    psg_seqs.append(str(psg_seq))
            else:
                previous_chrom = previous_chrom.split("_")[1]
                # Get exon sequences
                for exon in exon_coordinates:
                    psg_seq = haplotype_2[previous_chrom].seq[exon[0]-1:exon[1]]
                    psg_seqs.append(str(psg_seq))
            
            # Combine exon sequences into one sequence
            # Reverse complement if on reverse strand
            if previous_strand == "-":
                reversed_psg_seq_full = [str(Seq(psg_seq).reverse_complement()) for psg_seq in psg_seqs]
                psg_seq_full = "".join(reversed_psg_seq_full[::-1])
            else:
                psg_seq_full = Seq("".join(psg_seqs))

            # Save pseudogene sequence to fasta file
            psg_seq_full_fasta = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/psg_seq_full.fasta", "w")
            psg_seq_full_fasta.write(">psg_seq_full\n" + str(psg_seq_full))
            psg_seq_full_fasta.close()

            pseudogene_len = len(psg_seq_full)
            
            print("Aligning " + previous_gene + " to " + previous_hit + "..." + "\n")
            # Run MACSE to align gene to pseudogene
            # Frameshift penalty = 10000 for the gene and much lower penalty for the pseudogene (10)
            command = ["java", "-jar", "/home/ewcro/tools/macse_v2.06.jar", "-prog", "alignSequences", "-seq",
                       f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq.fasta", "-seq_lr", 
                       f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/psg_seq_full.fasta",
                       "-fs", "10000", "-fs_lr", "10", "-stop_lr", "15"]
            result = subprocess.run(command, stdout=subprocess.PIPE)

            # Read in the alignment
            aligned_cds = list(SeqIO.parse(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq_NT.fasta", "fasta"))
            gene_cds = aligned_cds[0].seq
            print(gene_cds, flush=True)
            pseudogene_cds = aligned_cds[1].seq
            print(pseudogene_cds, flush=True)

            # Calculate percent identity of pseudogene and gene DNA
            result_cds = calculate_percentages(gene_cds, pseudogene_cds)
            percent_identity_cds = result_cds[0]
            fraction_cds = result_cds[1]
            length_aligned_cds = result_cds[2]

            # Load in the aligned pseudogene and gene protein sequences
            aligned_prot = list(SeqIO.parse(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq_AA.fasta", "fasta"))
            print(aligned_prot, flush=True)
            # Extract the aligned pseudogene and gene protein sequences
            gene_prot = aligned_prot[0].seq
            print(gene_prot, flush=True)
            pseudogene_prot = aligned_prot[1].seq
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

            # check if start and end are present in pseudogene
            if pseudogene_cds[0:3].upper() == gene_cds[0:3].upper():
                start_codon = "present"
            else:
                start_codon = "absent"

            if pseudogene_cds[-3:].upper() in ["TAA", "TAG", "TGA"]:
                stop_codon = "present"
            else:
                stop_codon = "absent"

            string = previous_hit + "\t" + str(pseudogene_len) + "\t" + previous_gene + "\t" + str(gene_len) + "\t" + str(percent_identity_cds) + "\t" + str(fraction_cds) + "\t" + str(length_aligned_cds) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(num_fs) + "\t" + str(num_stop) + "\t" + start_codon + "\t" + stop_codon + "\n"

            # Write to file
            output = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/result_table_{haplotype1}_MACSE.tsv", "a")
            output.write(string)
            output.close()

            # Reset exon coordinates for new hit
            exon_coordinates = [] # reset exon coordinates
            exon_coordinates.append((hit_start, hit_end, chrom, strand))
            previous_hit = hit
            previous_gene = gene
            previous_chrom = chrom
            previous_strand = strand
            previous_hit_syn = hit_syn

# Obtain gene sequence
if previous_gene.startswith(haplotype1):
    gene_seq = haplotype_1_CDS[previous_gene].seq
else:
    gene_seq = haplotype_2_CDS[previous_gene].seq

# save gene sequence to fasta file
gene_seq_fasta = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq.fasta", "w")
gene_seq_fasta.write(">gene_seq\n" + str(gene_seq))
gene_seq_fasta.close()

gene_len = len(gene_seq)

# Obtain pseudogene sequence
psg_seqs = []
if previous_chrom.startswith(haplotype1):
    previous_chrom = previous_chrom.split("_")[1]
    # Get exon sequences
    for exon in exon_coordinates:
        psg_seq = haplotype_1[previous_chrom].seq[exon[0]-1:exon[1]]
        psg_seqs.append(str(psg_seq))
else:
    previous_chrom = previous_chrom.split("_")[1]
    # Get exon sequences
    for exon in exon_coordinates:
        psg_seq = haplotype_2[previous_chrom].seq[exon[0]-1:exon[1]]
        psg_seqs.append(str(psg_seq))

# Combine exon sequences into one sequence
# Reverse complement if on reverse strand
if previous_strand == "-":
    reversed_psg_seq_full = [str(Seq(psg_seq).reverse_complement()) for psg_seq in psg_seqs]
    psg_seq_full = "".join(reversed_psg_seq_full[::-1])
else:
    psg_seq_full = Seq("".join(psg_seqs))

# Save pseudogene sequence to fasta file
psg_seq_full_fasta = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/psg_seq_full.fasta", "w")
psg_seq_full_fasta.write(">psg_seq_full\n" + str(psg_seq_full))
psg_seq_full_fasta.close()

pseudogene_len = len(psg_seq_full)

print("Aligning " + previous_gene + " to " + previous_hit + "..." + "\n")
# Run MACSE to align gene to pseudogene
# Frameshift penalty = 10000 for the gene and much lower penalty for the pseudogene (10)
command = ["java", "-jar", "/home/ewcro/tools/macse_v2.06.jar", "-prog", "alignSequences", "-seq",
            f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq.fasta", "-seq_lr", 
            f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/psg_seq_full.fasta",
            "-fs", "10000", "-fs_lr", "10", "-stop_lr", "15"]
result = subprocess.run(command, stdout=subprocess.PIPE)

# Read in the alignment
aligned_cds = list(SeqIO.parse(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq_NT.fasta", "fasta"))
gene_cds = aligned_cds[0].seq
print(gene_cds, flush=True)
pseudogene_cds = aligned_cds[1].seq
print(pseudogene_cds, flush=True)

# Calculate percent identity of pseudogene and gene DNA
result_cds = calculate_percentages(gene_cds, pseudogene_cds)
percent_identity_cds = result_cds[0]
fraction_cds = result_cds[1]
length_aligned_cds = result_cds[2]

# Load in the aligned pseudogene and gene protein sequences
aligned_prot = list(SeqIO.parse(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/gene_seq_AA.fasta", "fasta"))

# Extract the aligned pseudogene and gene protein sequences
gene_prot = aligned_prot[0].seq
print(gene_prot, flush=True)
pseudogene_prot = aligned_prot[1].seq
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

# check if start and end are present in pseudogene
if pseudogene_cds[0:3].upper() == gene_cds[0:3].upper():
    start_codon = "present"
else:
    start_codon = "absent"

if pseudogene_cds[-3:].upper() in ["TAA", "TAG", "TGA"]:
    stop_codon = "present"
else:
    stop_codon = "absent"

string = previous_hit + "\t" + str(pseudogene_len) + "\t" + previous_gene + "\t" + str(gene_len) + "\t" + str(percent_identity_cds) + "\t" + str(fraction_cds) + "\t" + str(length_aligned_cds) + "\t" + str(percent_identity_prot) + "\t" + str(fraction_prot) + "\t" + str(length_aligned_prot) + "\t" + str(num_fs) + "\t" + str(num_stop) + "\t" + start_codon + "\t" + stop_codon + "\n"

# Write to file
output = open(f"{working_dir}/results/i-ADHoRe-run/processing/{multiplicon}/result_table_{haplotype1}_MACSE.tsv", "a")
output.write(string)
output.close()