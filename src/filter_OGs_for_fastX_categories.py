#### Parse the hierarchical orthogroups output from orthofinder to get the transcript IDs of orthologs of different categories

import parse_orthogroups as OGs
import parse_gff as gff
from Bio import SeqIO

def filepaths(username = "miltr339"):
    """
    get all filepaths with the correct username
    """
    tree = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder/species_tree.nw"
    orthogroups = f"/Users/{username}/work/chapter3/orthofinder/N0.tsv"

    proteins_dir=f"/Users/{username}/work/native_proteinseqs/"
    proteins_dict = {
        "A_obtectus" : f"{proteins_dir}A_obtectus.faa",
        "B_siliquastri" : f"{proteins_dir}B_siliquastri.faa",
        "C_chinensis" : f"{proteins_dir}C_chinensis.faa",
        "C_maculatus" : f"{proteins_dir}C_maculatus.faa",
        "D_carinulata" : f"{proteins_dir}D_carinulata.faa",
        "D_sublineata" : f"{proteins_dir}D_sublineata.faa",
    }

    nucleotides_dir=f"/Users/{username}/work/native_nucleotides/"
    nucleotides_dict = {
        "A_obtectus" : f"{nucleotides_dir}A_obtectus_transcripts.fna",
        "B_siliquastri" : f"{nucleotides_dir}B_siliquastri_transcripts.fna",
        "C_chinensis" : f"{nucleotides_dir}C_chinensis_transcripts.fna",
        "C_maculatus" : f"{nucleotides_dir}C_maculatus_transcripts.fna",
        "D_carinulata" : f"{nucleotides_dir}D_carinulata_transcripts.fna",
        "D_sublineata" : f"{nucleotides_dir}D_sublineata_transcripts.fna"
    }   

    annotations_dir=f"/Users/{username}/work/chapter3/isoform_filtered_native_annotations/"
    annotations_dict = {
        "A_obtectus" : f"{annotations_dir}A_obtectus.gff",
        "B_siliquastri" : f"{annotations_dir}B_siliquastri.gff",
        "C_chinensis" : f"{annotations_dir}C_chinensis.gff",
        "C_maculatus" : f"{annotations_dir}C_maculatus_superscaffolded.gff",
        "D_carinulata" : f"{annotations_dir}D_carinulata.gff",
        "D_sublineata" : f"{annotations_dir}D_sublineata.gff",
    }

    return tree, orthogroups, proteins_dict, nucleotides_dict, annotations_dict

def sex_chromosome_names():
    contig_names_dict = {
        "A_obtectus" : {
            "X" : ['CAVLJG010000002.1'],
            "Y" : ['scaffold_13', 'scaffold_86']
        },
        "B_siliquastri" : {
            "X" : ['X'],
            "Y" : ['Y']
        },
        "C_chinensis" : {
            "X" : ["1092_quiver","1124_quiver","1080_quiver","1105_quiver","1148_quiver","1339_quiver","1435_quiver","1482_quiver","1501_quiver","1565_quiver","1694_quiver","1688_quiver","1758_quiver","1618_quiver","1786_quiver","1816_quiver","1815_quiver","1817_quiver","1826_quiver","1889_quiver","1898_quiver","1908_quiver","1911_quiver","1933_quiver","2046_quiver","2054_quiver","2056_quiver","5713_quiver","2194_quiver","2226_quiver","2306_quiver","2357_quiver","2381_quiver","2392_quiver","2400_quiver","2435_quiver","2453_quiver","2513_quiver","2524_quiver","2569_quiver","2576_quiver","2580_quiver","2599_quiver","2693_quiver","2733_quiver","1210_quiver","2935_quiver","2958_quiver","2964_quiver","3034_quiver","3068_quiver","3080_quiver","3091_quiver"],
            "Y" : ["850_quiver","949_quiver","1088_quiver","1125_quiver","1159_quiver","1134_quiver","1224_quiver","1369_quiver","1410_quiver","1568_quiver","1577_quiver","1619_quiver","1634_quiver","1646_quiver","1652_quiver","1665_quiver","1681_quiver","1697_quiver","1722_quiver","1766_quiver","1783_quiver","1891_quiver","1937_quiver","1963_quiver","1790_quiver","1997_quiver","2073_quiver","2113_quiver","2163_quiver","2166_quiver","5705_quiver","2245_quiver","2259_quiver","2260_quiver","2334_quiver","2340_quiver","2382_quiver","2443_quiver","2511_quiver","2534_quiver","2573_quiver","2597_quiver","2651_quiver","2707_quiver","2766_quiver","2773_quiver","2791_quiver","2830_quiver","2875_quiver","3022_quiver","3070_quiver","3074_quiver","3075_quiver","3078_quiver"]
        },
        "C_maculatus" : {
            "X" : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
            "Y" : ['utg000322l_1','utg 000312c_1','utg 000610l_1','utg 001235l_1']},
        "D_carinulata" : {
            "X" : ['NC_079472.1'],
            "Y" : ['NC_079473.1']
        },
        "D_sublineata" : {
            "X" : ['NC_079485.1'],
            "Y" : ['NC_079486.1']
        },
    }

def associate_Aobt_contig_names(ENA_assembly, uppmax_assembly, minlen = 1e5):

    out_dict_ENA = {}
    for record in SeqIO.parse(ENA_assembly, "fasta"):
        contig_len = len(record.seq)
        if contig_len<minlen:
            continue
        else:
            header = record.id
            if "_" in header:
                contig = header.split("_Acanthoscelides_obtectus")[0].split("_")[-1]
            else:
                contig = header.split(" ")[0].split("|")[-1]
            out_dict_ENA[contig_len] = contig
    ENA_sorted_lengths = list(out_dict_ENA.keys())
    ENA_sorted_lengths.sort()
    print(f"{len(ENA_sorted_lengths)} contigs longer than {minlen} in the ENA assembly")
    
    out_dict_nat={}
    for record in SeqIO.parse(uppmax_assembly, "fasta"):
        contig_len = len(record.seq)
        if contig_len<minlen:
            continue
        else:
            out_dict_nat[contig_len] = record.id
            
    nat_sorted_lengths = list(out_dict_nat.keys())
    nat_sorted_lengths.sort()
    print(f"{len(nat_sorted_lengths)} contigs longer than {minlen} in the uppmax assembly")

    for i in range(len(nat_sorted_lengths)):
        print(f"{i}\n * {nat_sorted_lengths[i]}\t: {out_dict_nat[nat_sorted_lengths[i]]}\n * {ENA_sorted_lengths[i]}\t: {out_dict_ENA[ENA_sorted_lengths[i]]}")

    print(out_dict_ENA[43142905])
    print(out_dict_nat[43142905])


if __name__ == "__main__":

    tree, orthogroups, proteins_dict, nucleotides_dict, annotations_dict = filepaths()

    obtectus_dir="/Users/miltr339/work/a_obtectus/"
    associate_Aobt_contig_names(ENA_assembly=f"{obtectus_dir}A_obtectus_ENA_superscaffolded.fasta", uppmax_assembly=f"{obtectus_dir}A_obtectus_from_uppmax.fasta")