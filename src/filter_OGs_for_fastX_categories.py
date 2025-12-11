#### Parse the hierarchical orthogroups output from orthofinder to get the transcript IDs of orthologs of different categories

import parse_orthogroups as OGs
import parse_gff as gff
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
from collections import defaultdict
import upsetplot
import matplotlib.pyplot as plt
import numpy as np

def filepaths(username = "miltr339"):
    """
    get all filepaths with the correct username
    """
    tree = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder/species_tree.nw"
    orthogroups = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder/N0.tsv"

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
            "X" : ["1211|quiver","1844|quiver","854|quiver","5741|quiver","2866|quiver","658|quiver","1498|quiver","1455|quiver","2404|quiver","2935|quiver","1115|quiver","370|quiver","2273|quiver","1424|quiver","1865|quiver","767|quiver","2222|quiver","1525|quiver","5023|quiver","1925|quiver","1217|quiver","2328|quiver","2475|quiver","959|quiver","537|quiver","2776|quiver","325|quiver","2576|quiver","2336|quiver","988|quiver","2252|quiver","1388|quiver","1508|quiver","1712|quiver","1260|quiver","977|quiver","2202|quiver","2223|quiver","2397|quiver","693|quiver","1092|quiver","2189|quiver","1958|quiver","1355|quiver","2241|quiver","849|quiver","703|quiver","277|quiver","518|quiver","2589|quiver","1326|quiver","2962|quiver","2341|quiver","358|quiver","462|quiver","2786|quiver","1116|quiver","525|quiver","1358|quiver","5693|quiver","1429|quiver","1253|quiver","2372|quiver","326|quiver","474|quiver","777|quiver","955|quiver","1852|quiver","718|quiver","1024|quiver","1974|quiver","2295|quiver","2356|quiver","1484|quiver","1503|quiver","3076|quiver","2091|quiver","1262|quiver","1109|quiver","1475|quiver","1695|quiver","1168|quiver","1386|quiver","2201|quiver","2320|quiver","1117|quiver","769|quiver","2050|quiver","1805|quiver","2692|quiver","411|quiver","851|quiver","5703|quiver","1585|quiver","824|quiver","1816|quiver","1370|quiver","2416|quiver","1814|quiver","1277|quiver","619|quiver","1750|quiver","2709|quiver","2664|quiver","1250|quiver","971|quiver","3020|quiver","310|quiver","1176|quiver","2510|quiver","1699|quiver","1256|quiver","1420|quiver","5727|quiver","413|quiver","1124|quiver","682|quiver","1000|quiver","1313|quiver","5708|quiver","1556|quiver","274|quiver","1787|quiver","1137|quiver","360|quiver","1469|quiver","1853|quiver","2380|quiver","1239|quiver","993|quiver","791|quiver","2540|quiver","1510|quiver","868|quiver","505|quiver","1212|quiver","376|quiver","1564|quiver","1836|quiver","1670|quiver","500|quiver","2099|quiver","353|quiver","1042|quiver","419|quiver","1314|quiver","1339|quiver","1470|quiver","1576|quiver","1717|quiver","5692|quiver","2157|quiver","700|quiver","1284|quiver","1694|quiver","2306|quiver","2712|quiver","182|quiver","1973|quiver","882|quiver","2363|quiver","2482|quiver","1640|quiver","1913|quiver","2323|quiver","1240|quiver","161|quiver","1649|quiver","1164|quiver","1054|quiver","1096|quiver","313|quiver","1815|quiver","1831|quiver","1349|quiver","151|quiver","1478|quiver","1523|quiver","1888|quiver","739|quiver","1322|quiver","2338|quiver","1798|quiver","1391|quiver","1530|quiver","1519|quiver","1651|quiver","1105|quiver","509|quiver","1308|quiver","1833|quiver","1914|quiver","1741|quiver","1080|quiver","2292|quiver","2364|quiver","643|quiver","5745|quiver","1920|quiver","1725|quiver","125|quiver","1086|quiver","2552|quiver","5749|quiver","2120|quiver","2964|quiver","5722|quiver","2045|quiver","2422|quiver","593|quiver","1496|quiver","1772|quiver","799|quiver","2690|quiver","414|quiver","1531|quiver","1443|quiver","1408|quiver","1688|quiver","1371|quiver","1501|quiver","3090|quiver","1025|quiver","5698|quiver","347|quiver","1435|quiver","476|quiver","1883|quiver","2820|quiver","5728|quiver","342|quiver","1972|quiver","1826|quiver","968|quiver","2037|quiver","1723|quiver","252|quiver","1863|quiver","2983|quiver","1947|quiver","1430|quiver","1612|quiver","1701|quiver","839|quiver","613|quiver","1979|quiver","1584|quiver","2024|quiver","1486|quiver","1097|quiver"],
            # ["1092|quiver","1124|quiver","1080|quiver","1105|quiver","1148|quiver","1339|quiver","1435|quiver","1482|quiver","1501|quiver","1565|quiver","1694|quiver","1688|quiver","1758|quiver","1618|quiver","1786|quiver","1816|quiver","1815|quiver","1817|quiver","1826|quiver","1889|quiver","1898|quiver","1908|quiver","1911|quiver","1933|quiver","2046|quiver","2054|quiver","2056|quiver","5713|quiver","2194|quiver","2226|quiver","2306|quiver","2357|quiver","2381|quiver","2392|quiver","2400|quiver","2435|quiver","2453|quiver","2513|quiver","2524|quiver","2569|quiver","2576|quiver","2580|quiver","2599|quiver","2693|quiver","2733|quiver","1210|quiver","2935|quiver","2958|quiver","2964|quiver","3034|quiver","3068|quiver","3080|quiver","3091|quiver"],
            "Y" : ["850|quiver","949|quiver","1088|quiver","1125|quiver","1159|quiver","1134|quiver","1224|quiver","1369|quiver","1410|quiver","1568|quiver","1577|quiver","1619|quiver","1634|quiver","1646|quiver","1652|quiver","1665|quiver","1681|quiver","1697|quiver","1722|quiver","1766|quiver","1783|quiver","1891|quiver","1937|quiver","1963|quiver","1790|quiver","1997|quiver","2073|quiver","2113|quiver","2163|quiver","2166|quiver","5705|quiver","2245|quiver","2259|quiver","2260|quiver","2334|quiver","2340|quiver","2382|quiver","2443|quiver","2511|quiver","2534|quiver","2573|quiver","2597|quiver","2651|quiver","2707|quiver","2766|quiver","2773|quiver","2791|quiver","2830|quiver","2875|quiver","3022|quiver","3070|quiver","3074|quiver","3075|quiver","3078|quiver"]
        },
        # "C_maculatus_Kaufmann2023" : {
        #     "X" : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
        #     "Y" : ['utg000322l_1','utg000312c_1','utg000610l_1','utg001235l_1']
        # },
        "C_maculatus" : {
            "X" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'],
            "Y" : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']
        },
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




def get_OG_member_contigs(orthogroups_dict:dict, annotations_dict:dict, sex_chromosomes_dict:dict, max_GF_size:int = 30000, verbose = True):
    """
    Transform the orthogroups dict with the gene IDs into the same dict but with contig IDs
    the contig IDs come from the gff files in annotations_dict
    I also add a filtering step to get rid of gene families that have too many duplications already to make all future operations based on this quicker
    """

    annotations_class_dict = annotations_dict
    for species, annot_path in annotations_dict.items():
        print(f" * {species}")
        annotations_class_dict[species] = gff.parse_gff3_general(annot_path, keep_feature_category=gff.FeatureCategory.Transcript, verbose=False)

    ### TEST CMAC annotation
    # cmac_annot = annotations_class_dict["C_maculatus"]
    # print(cmac_annot["g14846.t1"])
    # raise RuntimeError

    orthogroups_contigs_dict = {}
    print(f"\n * assign chromosomes to transcript IDs")
    for OG_id, orthogroup_class in tqdm(orthogroups_dict.items()):

        species_dict = orthogroup_class.member_IDs_by_species
        
        if any(len(GF_list) > max_GF_size for GF_list in species_dict.values()):
            # print(f"{OG_id} excluded because at least one GF has more than {max_GF_size} members")
            continue
        
        genes_not_found_dict = {}
        for species, GF_transcript_list in species_dict.items():
            annot_dict = annotations_dict[species]
            sex_chr_dict = sex_chromosomes_dict[species]
            
            for tr_ID in GF_transcript_list:
                try:
                    transcript_contig = annot_dict[tr_ID].contig
                except:
                    raise RuntimeError(f"{species} transcript {tr_ID} (orthogroup {OG_id}) is not found in the given assembly for {species}: \n\t --> {annotations_dict[species]}")
                orthogroup_class.members[tr_ID].contig = transcript_contig
                autosome = True
                for chr_cat in ["X", "Y"]:
                    if transcript_contig in sex_chr_dict[chr_cat]:
                        orthogroup_class.members[tr_ID].chromosome_type = chr_cat
                        autosome = False
                        break
                # this is a bit sensitive because I don't actually have a dict that contains the autosomes.
                # it should be fine, since if transcript_contig can catch errors, but 
                if autosome:
                    orthogroup_class.members[tr_ID].chromosome_type = chr_cat = "A"

        orthogroups_contigs_dict[OG_id] = orthogroup_class

        if genes_not_found_dict !={} and verbose:
            print(f"* {OG_id} genes present in the orthofinder run not found in the annotation: {genes_not_found_dict}")
        # raise RuntimeError

    return orthogroups_contigs_dict
    # return orthogroups_dict


def write_orthogroups_contigs_to_file(orthogroups:dict, sex_chromosomes_contigs_file:str = "orthogroups_by_contig.tsv", species_list = [], binary_file = False, species_exclude = ""):
    """
    make a tsv file for all the 1-to-1 orthogroups that contains the contig positions instead of the transcript IDs
    You can make a file that contains the chromosome IDs, or a binary file that has 1 and 0 for presence/absence on the X chromosome
    """
    if species_list ==[]:
        ## get species list
        for orthogroup in orthogroups.values():
            member_species = orthogroup.member_species_list
            species_list.extend(member_species)
        species_list = list(set(member_species))
        print(f"species list undefined, complete species list inferred from orthogroups data: {species_list}")
    else:
        print(f"species list given: {species_list} ({len(species_list)} members)")
    
    if species_exclude != "" and species_exclude in species_list:
        species_list.remove(species_exclude)
    elif species_exclude != "" and species_exclude not in species_list:
        print(f"{species_exclude} was not excluded because it is not in {species_list}")
    elif species_exclude == "":
        print(f"all species included")

    header_t = "\t".join(species_list)
    header = f"orthogroup\t{header_t}\n"

    if binary_file:
        chr_dict = {
        "X" :"1",
        "A" :"0",
        "Y" :"0",
        "None": "0"}
    else:
        chr_dict = {
        "X" :"X",
        "A" :"A",
        "Y" :"Y",
        "None": "None"}

    with open(sex_chromosomes_contigs_file, "w") as outfile:
        outfile.write(header)
        for orthogroup in orthogroups.values():

            # only include 1-to-1 orthologs
            if orthogroup.is_one_to_one == False:
                continue

            OG_species_list = orthogroup.member_species_list
            # only include 1-to-1 orthologs that have a member in all species in the species list
            if len(set(species_list) - set(OG_species_list))>0:
                # print(f"{orthogroup.ID} does not contain all species")
                continue

            OG_dict = {orthogroup_member.species : orthogroup_member.chromosome_type for orthogroup_member in orthogroup.members.values()}
            # if len(OG_dict)<len(species_list):
            #     print(orthogroup)
            #     print(f"{orthogroup.ID} does not contain all species")
            OG_list = [chr_dict[OG_dict[species]] for species in species_list]
            OG_line = "\t".join(OG_list)
            whole_line = f"{orthogroup.ID}\t{OG_line}\n"
            outfile.write(whole_line)
    
    print(f"orthogroup contigs outfile written to {sex_chromosomes_contigs_file}")

    return sex_chromosomes_contigs_file

    
def plot_orthogroup_sex_chromosome_table(orthogroups_sex_chromosomes_file, filter_A_exclusive = True):
    # (!) the input file has to be generated with the function: write_orthogroups_contigs_to_file
    # (!) in binary mode, set with bin_file = True

    data = pd.read_csv(orthogroups_sex_chromosomes_file, sep = "\t")

    if filter_A_exclusive:
        # Filter out all orthogroups that are autosomal in all species
        data = data[data.iloc[:, 1:].sum(axis=1) > 0]

    # make a dictionary that can be used to plot the data
    presence_dict = defaultdict(int)
    for index, row in data.iterrows():
        conditions = tuple(row[1:])  # all columns except the first 
        presence_dict[conditions] += 1

    # prepare input for the upset plot library by transforming it to a pandas series
    # use MultiIndex to keep the column names
    multi_index = pd.MultiIndex.from_tuples(presence_dict.keys(), names=data.columns[1:])
    upset_data = pd.Series(list(presence_dict.values()), index=multi_index)

    # print(upset_data)

    upsetplot.plot(upset_data, show_counts=True)

    # Plotting
    plot_name = orthogroups_sex_chromosomes_file.split(".")[0]+"_upset_plot.png"
    plt.suptitle("X chromosome presence of single-copy orthologs")
    plt.savefig(plot_name, dpi = 300, transparent = False)
    print(f"plot saved here: {plot_name}") 

if __name__ == "__main__":

    username = "miltr339"
    data_dir = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder"
    tree, orthogroups_path, proteins_dict, nucleotides_dict, annotations_dict = filepaths(username=username)
    tree_species_list = gff.make_species_order_from_tree(tree)
    sex_chr_contigs_dict = sex_chromosome_names()
    binary_file_for_upset = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder/orthogroups_by_contig.tsv"

    if False:
        obtectus_dir="/Users/miltr339/work/a_obtectus/"
        aobt_xlist = ["chr_10","scaffold_49","scaffold_77","scaffold_108","scaffold_113","scaffold_121","scaffold_133","scaffold_143","scaffold_176","scaffold_186","scaffold_188","scaffold_192","scaffold_200","scaffold_207","scaffold_219","scaffold_227","scaffold_246","scaffold_276","scaffold_319","scaffold_327","scaffold_328","scaffold_341","scaffold_356","scaffold_363","scaffold_365","scaffold_370","scaffold_408","scaffold_411","scaffold_419","scaffold_420","scaffold_435","scaffold_482","scaffold_507","scaffold_524","scaffold_547","scaffold_563","scaffold_589","scaffold_602","scaffold_604","scaffold_621","scaffold_630","scaffold_633","scaffold_676","scaffold_697","scaffold_734","scaffold_768","scaffold_803","scaffold_838","scaffold_840","scaffold_855","scaffold_1045","scaffold_1086","scaffold_1100","scaffold_1154","scaffold_1176","scaffold_1195","scaffold_1209","scaffold_1267","scaffold_1338","scaffold_1339","scaffold_1356","scaffold_1498","scaffold_1564","scaffold_1663","scaffold_1704","scaffold_1759","scaffold_1786","scaffold_1796","scaffold_1822","scaffold_1875","scaffold_1902","scaffold_1913","scaffold_1914","scaffold_1922","scaffold_1949","scaffold_1956","scaffold_1988","scaffold_2012","scaffold_2027","scaffold_2033","scaffold_2041","scaffold_2045","scaffold_2061","scaffold_2071","scaffold_2101","scaffold_2107","scaffold_2124","scaffold_2144","scaffold_2194","scaffold_2225","scaffold_2265","scaffold_2289","scaffold_2371","scaffold_2372","scaffold_2403","scaffold_2469","scaffold_2509","scaffold_2524"]
        aobt_ylist = ["scaffold_13","scaffold_36","scaffold_120","scaffold_150","scaffold_152","scaffold_265","scaffold_284","scaffold_287","scaffold_303","scaffold_618","scaffold_1035","scaffold_1472","scaffold_1594","scaffold_1702","scaffold_2097","scaffold_2445"]
        print_contig_names_lengths(ENA_assembly=f"{obtectus_dir}A_obtectus_ENA_superscaffolded.fasta", minlen=100e3, xlist=aobt_xlist, ylist=aobt_ylist)

    
    ## make contigs file
    if True:
        orthogroups = OGs.parse_orthogroups_class(orthogroups_path, verbose=True)
        ## assign sex chromosomes
        orthogroups_contigs_dict = get_OG_member_contigs(orthogroups, annotations_dict, sex_chr_contigs_dict)# , max_GF_size = 2)
        if len(orthogroups) == len(orthogroups_contigs_dict):
            print(f"\n{len(orthogroups)} orthogroups in original file, {len(orthogroups_contigs_dict)} in filtered file with contigs\n")

        binary_file_for_upset = write_orthogroups_contigs_to_file(orthogroups=orthogroups_contigs_dict, sex_chromosomes_contigs_file = f"{data_dir}/orthogroups_by_contig.tsv", species_list=tree_species_list, binary_file = True)
        # binary_file_for_upset = write_orthogroups_contigs_to_file(orthogroups=orthogroups_contigs_dict, sex_chromosomes_contigs_file = f"{data_dir}/orthogroups_by_contig.tsv", species_list=['A_obtectus', 'B_siliquastri', 'C_maculatus', 'C_chinensis'], binary_file = True)
        # binary_file_for_upset = write_orthogroups_contigs_to_file(orthogroups=orthogroups_contigs_dict, sex_chromosomes_contigs_file = f"{data_dir}/orthogroups_by_contig.tsv", species_list=['D_carinulata', 'D_sublineata'], binary_file = True)
        
        if False:
            OG_example = "N0.HOG0005248"
            print(orthogroups_contigs_dict[OG_example])
            species_example = orthogroups_contigs_dict[OG_example].member_species_list
            print(species_example)

        if False:
            count_gametologs = 0
            count_onetoone = 0
            count_onetoone_X = 0
            count_onetoone_A = 0
            count_onetoone_mixed = 0
            for OG_id, orthogroup in orthogroups_contigs_dict.items():
                if orthogroup.has_gametolog:
                    # print(orthogroup)
                    count_gametologs += 1
                if orthogroup.is_one_to_one:
                    count_onetoone += 1
                    # try:
                    #     orthogroup_no_obt = orthogroup.exclude_species("A_obtectus")
                    # except:
                    #     orthogroup_no_obt = orthogroup

                    if orthogroup.is_on_chr_type("X", exclusive = True):
                        count_onetoone_X += 1
                    elif orthogroup.is_on_chr_type("A", exclusive = True):
                        count_onetoone_A += 1
                    else:
                        count_onetoone_mixed += 1

            print(f"{count_gametologs} orthogroups contain gametologs")
            print(f"""
{count_onetoone} 1-to-1 orthologs,
    * {count_onetoone_X} exclusively X-linked
    * {count_onetoone_A} exclusively A-linked
    * {count_onetoone_mixed} from different chromosome types
""")


    ## make upset plot
    if True:
        plot_orthogroup_sex_chromosome_table(binary_file_for_upset, filter_A_exclusive = True)


