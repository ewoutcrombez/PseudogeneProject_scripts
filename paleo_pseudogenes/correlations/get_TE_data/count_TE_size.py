import os

# I need to first remove overlapping TEs in the gff file
# and combine the overlapping TEs into one TE

# combine overlapping TEs into one TE


def count_TE_size(species):
    """
    count the size of TE in the gff file
    """

    print("Processing species: ", species)

    TE_size = []
    TE_number = 0
    LINE_size = []
    LINE_number = 0
    SINE_size = []
    SINE_number = 0
    LTR_size = []
    LTR_number = 0

    with open("%s_TE.bed"%species, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                line = line.split('\t')
                TE_size += [int(line[2]) - int(line[1])]
                TE_number += 1
    
    with open("%s_LTR.bed"%species, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                line = line.split('\t')
                LTR_size += [int(line[2]) - int(line[1])]
                LTR_number += 1

    with open("%s_LINE.bed"%species, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                line = line.split('\t')
                LINE_size += [int(line[2]) - int(line[1])]
                LINE_number += 1
    
    with open("%s_SINE.bed"%species, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                line = line.split('\t')
                SINE_size += [int(line[2]) - int(line[1])]
                SINE_number += 1

    retro_size = sum(LTR_size) + sum(LINE_size) + sum(SINE_size)
    retro_number = LTR_number + LINE_number + SINE_number
    string = species + "\t" + str(sum(TE_size) / 10**6) + "\t" + str(TE_number) + "\t" + str(sum(LTR_size) / 10**6) + "\t" + str(LTR_number) + "\t" + str(sum(LINE_size) / 10**6) + "\t" + str(LINE_number) + "\t" + str(sum(SINE_size) / 10**6) + "\t" + str(SINE_number) + "\t" + str(retro_size / 10**6) + "\t" + str(retro_number) + "\n"
    return string 

with open("TE_size.txt", 'w') as f:
        f.write("species\tTE_size\tTE_number\tLTR_size\tLTR_number\tLINE_size\tLINE_number\tSINE_size\tSINE_number\tretro_size\tretro_number\n")

for species in ["amborella_trichopoda", "arabidopsis_thaliana",
                "brassica_oleracea", "brassica_rapa", "glycine_max",
                "oryza_sativa", "populus_trichocarpa", "solanum_lycopersicum",
                "sorghum_bicolor", "vitis_vinifera", "zea_mays"]:
    string = count_TE_size(species)
    with open("TE_size.txt", 'a') as f:
        f.write(string)