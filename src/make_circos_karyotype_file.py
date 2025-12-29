### make circos karyotype file

## load these modules to run on uppmax:
# module load bioinfo-tools python3 biopython/1.80-py3.10.8


## Formatting:

# name, label, start and end position and a color
# chr - ID LABEL START END COLOR
# for example:
# chr - hs1 1 0 249250621 chr1
# chr - hs2 2 0 243199373 chr2
# chr - hs3 3 0 198022430 chr3


## Summary:

# This file makes the above explained format directly from a fasta file.
# it assumes a good assembly, at least approximating chromosome level since it does not contain any stringent filtering criteria for the contigs

# the script makes individual karyotype files per species, but for the synteny stuff you probably need one that contains two species
# you can just merge them with stuff like "cat file1 file2 > merged_file"


from Bio import SeqIO
import pandas as pd
import os
from collections import Counter


##########################################
######## lists of X and Y contigs ########
##########################################

cmac_x_contigs = ["utg000057l_1", "utg000114l_1", "utg000139l_1", "utg000191l_1", "utg000326l_1", "utg000359l_1", "utg000532l_1", "utg000602l_1"]
cmac_X_superscaffold = "scaffold_10"
obtectus_X = "CAVLJG010000002.1_Acanthoscelides_obtectus_genome_assembly,_contig:_chr_10,_whole_genome_shotgun_sequence"
siliquastri_X = "X_softmasked:chromosome_primary_assembly:icBruSili1.1:X:1:22912655:1 chr_X_siliquastri"

X_list = [cmac_X_superscaffold, obtectus_X, siliquastri_X]


# contig outfile
circos_outfile = "circos_karyotype_chinensis.txt"
# add labels to the contigs
circos_labels_outfile = "circos_karyotype_chinensis_labels.txt"


def make_karyotype_file_no_colors(assembly, min_contig_length = 5000000, color = "lorange"):
    outfile_prefix = assembly.split(".fna")[0].split("/")[-1]
    with open(f"{outfile_prefix}_karyotype.txt", "w") as circos_outfile:
        karyotype_entries = []
        for record in SeqIO.parse(assembly, "fasta"):
            if len(record.seq) > min_contig_length:
             
                label = record.id
                entry = ["chr - "+record.id.split(":")[0], label, "0", str(len(record.seq)), color] # split the record.id at ":" because that confuses circos
                
                karyotype_entries.append(entry)
                circos_outfile.write(" ".join(entry)+"\n")

    print(f"karyotype saved in current working directory as: {outfile_prefix}_karyotype.txt")


def assemblies(username = "miltr339"):
    # ass_dir = f"/Users/{username}/work/assemblies_masked_uniform"
    ass_dir = "/proj/naiss2023-6-65/Milena/chapter2/assemblies"
    ass_dict = {
        "A_obtectus" : f"{ass_dir}/A_obtectus_assembly.fna",
        "B_siliquastri" : f"{ass_dir}/B_siliquastri_assembly.fna",
        "C_chinensis" : f"{ass_dir}/C_chinensis_assembly.fna",
        "C_maculatus" : f"{ass_dir}/C_maculatus_superscaffolded_assembly.fna",
        "D_carinulata" : f"{ass_dir}/D_carinulata_assembly.fna",
        "D_sublineata" : f"{ass_dir}/D_sublineata_assembly.fna",
    }
    return ass_dict


def X_lists():
    lists_dict = {
        "A_obtectus" : ["CAVLJG010000002.1","CAVLJG010003236.1","CAVLJG010003544.1","CAVLJG010000099.1","CAVLJG010000155.1","CAVLJG010000244.1","CAVLJG010000377.1","CAVLJG010000488.1",],
        "B_siliquastri" : ['X'],
        "C_chinensis" : ["1092_quiver","1124_quiver","1080_quiver","1105_quiver","1148_quiver","1339_quiver","1435_quiver","1482_quiver","1501_quiver","1565_quiver","1694_quiver","1688_quiver","1758_quiver","1618_quiver","1786_quiver","1816_quiver","1815_quiver","1817_quiver","1826_quiver","1889_quiver","1898_quiver","1908_quiver","1911_quiver","1933_quiver","2046_quiver","2054_quiver","2056_quiver","5713_quiver","2194_quiver","2226_quiver","2306_quiver","2357_quiver","2381_quiver","2392_quiver","2400_quiver","2435_quiver","2453_quiver","2513_quiver","2524_quiver","2569_quiver","2576_quiver","2580_quiver","2599_quiver","2693_quiver","2733_quiver","1210_quiver","2935_quiver","2958_quiver","2964_quiver","3034_quiver","3068_quiver","3080_quiver","3091_quiver"],
        "C_maculatus": ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
        "D_carinulata" : ['NC_079472.1'],
        "D_sublineata" : ['NC_079485.1']
    }
    return lists_dict

def make_karyotype_file(assembly, outfile_name = "", min_contig_length = 0, X_list = []):
    if outfile_name == "":
        outfile_prefix = assembly.split(".fna")[0].split("/")[-1]
        outfile_name = f"{outfile_prefix}_karyotype.txt"
    with open(outfile_name, "w") as circos_outfile:
        karyotype_entries = []
        colors = {
            "A" : "lorange",
            "X" : "acen"
        }
        for record in SeqIO.parse(assembly, "fasta"):
            if len(record.seq) > min_contig_length:
                label = record.id.split(":")[0]

                # match for X_contig in label (so if X = scaffold_10  then scaffold_100 etc. are also mtched)
                # if any(record.id in X_contig for X_contig in X_list):

                # match X names exactly
                if label in X_list:
                    entry = [f"chr - {label}", record.id, "0", str(len(record.seq)), colors["X"]] # split the record.id at ":" because that confuses circos
                else:
                    entry = [f"chr - {label}", record.id, "0", str(len(record.seq)), colors["A"]] # split the record.id at ":" because that confuses circos
                karyotype_entries.append(entry)
                circos_outfile.write(" ".join(entry)+"\n")

    print(f"karyotype saved in current working directory as: {outfile_prefix}_karyotype.txt")


if __name__ == "__main__":

    minlen = 5e6
    assemblies_dict = assemblies()
    X_lists_dict = X_lists()    
    outdir = "/proj/naiss2023-6-65/Milena/chapter2/pairwise_wga/circos"
    
    for species, assembly in assemblies_dict.items():
        outfile_name = f"{outdir}/{species}_karyotype.txt"
        print(outfile_name)
        make_karyotype_file(assembly=assembly, outfile_name=outfile_name,min_contig_length=minlen, X_list=X_lists_dict[species])
        print(f"\t--> done!")