#### Parse the hierarchical orthogroups output from orthofinder to get the transcript IDs of orthologs of different categories

import parse_orthogroups as OGs
import parse_gff as gff
from Bio import SeqIO
from tqdm import tqdm

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
            "X" : ["CAVLJG010000002.1","CAVLJG010003236.1","CAVLJG010003544.1","CAVLJG010000099.1","CAVLJG010000155.1","CAVLJG010000244.1","CAVLJG010000377.1","CAVLJG010000488.1",],
            "Y" : ["CAVLJG010000343.1","CAVLJG010002896.1","CAVLJG010000233.1","CAVLJG010000566.1","CAVLJG010000588.1",]
        },
        "B_siliquastri" : {
            "X" : ['X'],
            "Y" : ['Y']
        },
        "C_chinensis" : {
            "X" : ["1092|quiver","1124|quiver","1080|quiver","1105|quiver","1148|quiver","1339|quiver","1435|quiver","1482|quiver","1501|quiver","1565|quiver","1694|quiver","1688|quiver","1758|quiver","1618|quiver","1786|quiver","1816|quiver","1815|quiver","1817|quiver","1826|quiver","1889|quiver","1898|quiver","1908|quiver","1911|quiver","1933|quiver","2046|quiver","2054|quiver","2056|quiver","5713|quiver","2194|quiver","2226|quiver","2306|quiver","2357|quiver","2381|quiver","2392|quiver","2400|quiver","2435|quiver","2453|quiver","2513|quiver","2524|quiver","2569|quiver","2576|quiver","2580|quiver","2599|quiver","2693|quiver","2733|quiver","1210|quiver","2935|quiver","2958|quiver","2964|quiver","3034|quiver","3068|quiver","3080|quiver","3091|quiver"],
            "Y" : ["850|quiver","949|quiver","1088|quiver","1125|quiver","1159|quiver","1134|quiver","1224|quiver","1369|quiver","1410|quiver","1568|quiver","1577|quiver","1619|quiver","1634|quiver","1646|quiver","1652|quiver","1665|quiver","1681|quiver","1697|quiver","1722|quiver","1766|quiver","1783|quiver","1891|quiver","1937|quiver","1963|quiver","1790|quiver","1997|quiver","2073|quiver","2113|quiver","2163|quiver","2166|quiver","5705|quiver","2245|quiver","2259|quiver","2260|quiver","2334|quiver","2340|quiver","2382|quiver","2443|quiver","2511|quiver","2534|quiver","2573|quiver","2597|quiver","2651|quiver","2707|quiver","2766|quiver","2773|quiver","2791|quiver","2830|quiver","2875|quiver","3022|quiver","3070|quiver","3074|quiver","3075|quiver","3078|quiver"]
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
    return contig_names_dict

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

    # print(out_dict_ENA[43142905])
    # print(out_dict_nat[43142905])


def print_contig_names_lengths(ENA_assembly, minlen = 1e5, xlist=[], ylist=[]):

    out_dict_ENA = {}
    for record in SeqIO.parse(ENA_assembly, "fasta"):
        contig_len = len(record.seq)
        header = record.description # description gives the whole fasta header even after spaces, which name and id leave out
        contig = header.split("contig: ")[-1]
        ENA_name = header.split(" ")[0].split("|")[-1]
        out_dict_ENA[contig] = [contig_len, ENA_name]
    
    if len(xlist)>0:
        print(f"X-contigs:")
        for x_contig in xlist:
            out_list = out_dict_ENA[x_contig]
            contig_len = out_list[0]
            if contig_len<minlen:
                continue
            ENA_name = out_list[1]
            print(f"* {x_contig} : {ENA_name}, {contig_len} bp")
    
    if len(ylist)>0:
        print(f"\n========================================================================\n")
        print(f"Y-contigs:")
        for y_contig in ylist:
            out_list = out_dict_ENA[y_contig]
            contig_len = out_list[0]
            if contig_len<minlen:
                continue
            ENA_name = out_list[1]
            print(f"* {y_contig} : {ENA_name}, {contig_len} bp")


def filter_orthogroups_dict(orthogroups_contigs_dict:dict):
    """
    filter the orthogroups dict for the fastX analysis. split into several dictionaries:
    * one_to_one : 1-to-1 orthologs present in all species, regardless of chromosomal position 
    * gametologs : orthogroups with two gene family members in each species
    """


def get_OG_member_contigs(orthogroups_dict:dict, annotations_dict:dict, max_GF_size:int = 30000):
    """
    Transform the orthogroups dict with the gene IDs into the same dict but with contig IDs
    the contig IDs come from the gff files in annotations_dict
    I also add a filtering step to get rid of gene families that have too many duplications already to make all future operations based on this quicker
    """

    annotations_class_dict = annotations_dict
    for species, annot_path in annotations_dict.items():
        print(f" * {species}")
        annotations_class_dict[species] = gff.parse_gff3_general(annot_path, keep_feature_category=gff.FeatureCategory.Transcript, verbose=False)

    orthogroups_contigs_dict = {}
    print(f"\n... replace all gene IDs with contigs")
    for OG_id, species_dict in tqdm(orthogroups_dict.items()):
        if any(len(GF_list) > max_GF_size for GF_list in species_dict.values()):
            # print(f"{OG_id} excluded because at least one GF has more than {max_GF_size} members")
            continue
        
        orthogroups_contigs_dict[OG_id] = {}
        for species, GF_transcript_list in species_dict.items():
            formatted_transcripts = [tr_ID[:-2] if tr_ID[-2:] == "_1" else tr_ID for tr_ID in GF_transcript_list]
            try:
                contigs = [annotations_class_dict[species][tr_ID].contig if tr_ID != "" else tr_ID for tr_ID in formatted_transcripts]
            except:
                print(f"{OG_id}, {species}: gene IDs not found {formatted_transcripts}")
            orthogroups_contigs_dict[OG_id][species] = contigs
        # raise RuntimeError
    return orthogroups_contigs_dict


if __name__ == "__main__":

    tree, orthogroups_path, proteins_dict, nucleotides_dict, annotations_dict = filepaths()

    if False:
        obtectus_dir="/Users/miltr339/work/a_obtectus/"
        aobt_xlist = ["chr_10","scaffold_49","scaffold_77","scaffold_108","scaffold_113","scaffold_121","scaffold_133","scaffold_143","scaffold_176","scaffold_186","scaffold_188","scaffold_192","scaffold_200","scaffold_207","scaffold_219","scaffold_227","scaffold_246","scaffold_276","scaffold_319","scaffold_327","scaffold_328","scaffold_341","scaffold_356","scaffold_363","scaffold_365","scaffold_370","scaffold_408","scaffold_411","scaffold_419","scaffold_420","scaffold_435","scaffold_482","scaffold_507","scaffold_524","scaffold_547","scaffold_563","scaffold_589","scaffold_602","scaffold_604","scaffold_621","scaffold_630","scaffold_633","scaffold_676","scaffold_697","scaffold_734","scaffold_768","scaffold_803","scaffold_838","scaffold_840","scaffold_855","scaffold_1045","scaffold_1086","scaffold_1100","scaffold_1154","scaffold_1176","scaffold_1195","scaffold_1209","scaffold_1267","scaffold_1338","scaffold_1339","scaffold_1356","scaffold_1498","scaffold_1564","scaffold_1663","scaffold_1704","scaffold_1759","scaffold_1786","scaffold_1796","scaffold_1822","scaffold_1875","scaffold_1902","scaffold_1913","scaffold_1914","scaffold_1922","scaffold_1949","scaffold_1956","scaffold_1988","scaffold_2012","scaffold_2027","scaffold_2033","scaffold_2041","scaffold_2045","scaffold_2061","scaffold_2071","scaffold_2101","scaffold_2107","scaffold_2124","scaffold_2144","scaffold_2194","scaffold_2225","scaffold_2265","scaffold_2289","scaffold_2371","scaffold_2372","scaffold_2403","scaffold_2469","scaffold_2509","scaffold_2524"]
        aobt_ylist = ["scaffold_13","scaffold_36","scaffold_120","scaffold_150","scaffold_152","scaffold_265","scaffold_284","scaffold_287","scaffold_303","scaffold_618","scaffold_1035","scaffold_1472","scaffold_1594","scaffold_1702","scaffold_2097","scaffold_2445"]
        print_contig_names_lengths(ENA_assembly=f"{obtectus_dir}A_obtectus_ENA_superscaffolded.fasta", minlen=100e3, xlist=aobt_xlist, ylist=aobt_ylist)

    ### TODO 
    species_list = gff.make_species_order_from_tree(tree)
    print(species_list)
    orthogroups = OGs.parse_orthogroups_dict(orthogroups_path)
    orthogroups_contigs_dict = get_OG_member_contigs(orthogroups, annotations_dict, max_GF_size = 2)
    print(f"{len(orthogroups)} orthogroups in original file, {len(orthogroups_contigs_dict)} in filtered file with contigs")
    OG_ex = orthogroups["N0.HOG0000000"]
    print(f"example N0.HOG0000000: {OG_ex}")
    cont_ex = orthogroups_contigs_dict["N0.HOG0000000"]
    print(f"example N0.HOG0000000: {cont_ex}")

