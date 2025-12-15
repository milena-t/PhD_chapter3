# make barplot that shows the proportion of single-exon genes in species-specific genes compared to the proportion in the genome overall

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re
import numpy as np
from pathlib import Path
import parse_gff as gff
from collections import Counter
from Bio import SeqIO



def split_at_second_occurrence(s, char = "_"): # split the gene string at the second occurence of "_" to get only the species name
    if s.count(char)<2:
        return s
    else:
        first_occurrence = s.find(char)
        second_occurrence = s.find(char, first_occurrence + 1)
        species = s[:second_occurrence]
        return species

def make_species_order_from_tree(newick_tree_path):
    # Regular expression to extract leaf names
    # This matches strings between commas, parentheses, and before colons.
    leaf_pattern = r'(?<=\(|,)([a-zA-Z0-9_]+)(?=:)'
    with open(newick_tree_path, "r") as newick_tree_file:
        newick_tree_string = newick_tree_file.readlines()[0]
        # print(newick_tree_string)
        leaf_names = re.findall(leaf_pattern, newick_tree_string)
        # leaf names are like "A_obtectus_filtered_proteinfasta" but we only care about the species names in the beginning
        species_names = [split_at_second_occurrence(leaf, "_") for leaf in leaf_names]
    return species_names


def get_single_exon_genes(dir_path, species_list, write_to_file = False, outfile_name = "", include_total_gene_num = False):
    """
    get a dictionary with {species_name : [list, of, single, exon, transcripts]} from files generated in 04c_make_single_exon_gene_lists.sh
    """

    directory = Path(dir_path)
    # List all files in the directory
    all_files = [f for f in directory.iterdir() if f.is_file()]

    species_dict = {}
    num_single_exon_genes = {}
    num_transcripts = {}

    all_break = False
    for infile_path in all_files:

        print_statements = False
        infile_path = str(infile_path)
        species = infile_path.split("/")[-1]
        species = species.replace("native_", "")
        species = species.replace("orthoDB_", "")
        species = split_at_second_occurrence(species, "_")
        if len(species.split("_")[0])>1:
            species = species.split("_")[0][0].upper()+"_"+species.split("_")[-1]
        # replace with exact names from the tree if the species exists in the tree 
        species_mod = next((species_name for species_name in species_list if species in species_name), None) # replace the above parsed species name with the corresponding name from species_list
        if species_mod!=None:
            species = species_mod
        
        if len(species.split("_")[0])>1 and len(species.split("_"))>1:
            # if the species name is tribolium_castaneum instead of T_castaneum for example
            parts = species.split("_")
            species = f"{parts[0][0].upper()}_{parts[1]}"

        infile_dict = {}
        list_key = ""
        num_key = ""
        transcripts_num_key = ""
        with open(infile_path, "r") as infile:
            dict_elements = infile.readlines()
            for species_line in dict_elements:
                species_line = species_line.strip()
                try:
                    key = species_line.split(":")[0]
                    value = ":".join(species_line.split(":")[1:]) # some of the transcript names have ":" in them to spite me personally
                except:
                    print_statements = True
                    print(species_line[0:150])
                    continue
                value_list = [gene.split("|")[-1] for gene in value.strip().split(",")] # TODO # native annotations have some stuff going on before "|" and to actually make it work you need to get rid of that
                infile_dict[key]=value_list
                if "list" in key:
                    list_key = key
                if "single exon transcripts" in key:
                    num_key = key
                #if "number of transcripts" in key:
                if "gene features" in key:    
                    transcripts_num_key = key
                
 
            if all_break:
                break
            if print_statements:
                print(species)
                print(list_key)
                print(infile_dict.keys())
                print_statements = False
        # species_dict[species] = infile_dict["list of single-exon transcript IDs"]
        try:
            species_dict[species] = infile_dict[list_key]
        except:
            print(f"{species} single-exon list did not work with {list_key}")
        
        try:
            num_single_exon_genes[species] = int(infile_dict[num_key][0])
        except:
            infile_dict_small = {key: value[0][:150] for key, value in infile_dict.items()}
            print(infile_dict_small)
            print(f"{species} single-exon number did not work with {num_key}")
            raise RuntimeError

        if include_total_gene_num:
            try:
                num_transcripts[species] = int(infile_dict[transcripts_num_key][0])
            except:
                keys = list(infile_dict.keys())
                print(f"{species} transcript number did not work with {transcripts_num_key}, available keys are: {keys}")

    if write_to_file:
        gff.write_dict_to_file(species_dict, f"{dir_path}/../{outfile_name}")

    # num_single_exon_genes = {species : len(IDs) for species, IDs in species_dict.items()}
    if not include_total_gene_num:
        return(species_dict, num_single_exon_genes)
    elif include_total_gene_num:
        return(species_dict, num_single_exon_genes, num_transcripts)


def get_species_specific_genes(orthogroups_filepath, species_list, write_to_file = False, outfile_name = ""):
    # get a dictionary with {species_name : [list, of, species-specific, transcripts]}    

    list_spec_genes = {} # {species : [list, of, species-specific, transcripts]}
    num_spec_genes = {} # dictionary of species-specific genes {species : number_of_species-specific orthogroups}

    with open(orthogroups_filepath+"/Orthogroups.txt", "r") as orthogroups_file:
        orthogroups_lines = orthogroups_file.readlines()
        
        complete_gene_count_dict = {}
        complete_gene_count = 0

        # set up dictionary
        for species in species_list:
            list_spec_genes[species] = []
            complete_gene_count_dict[species] = 0
            
        for orthogroup in orthogroups_lines:
            orthogroup_members = orthogroup.split(": ")[1].strip().split(" ")
            complete_gene_count = complete_gene_count+len(orthogroup_members)

            # check if it's a single-species orthogroup
            # og_species_list = list(set([split_at_second_occurrence(gene, "_") for gene in orthogroup_members]))
            # og_species_list = [gene.rsplit("_",1)[0] for gene in orthogroup_members]
            og_species_list = [split_at_second_occurrence(gene, "_") for gene in orthogroup_members]
            og_species_list = [next((species_name for species_name in species_list if species_name in gene_species), gene_species) for gene_species in og_species_list]
            
            # count number of species in all genes
            all_species_counter = Counter(og_species_list)
            for species in species_list:
                complete_gene_count_dict[species] = complete_gene_count_dict[species]+all_species_counter[species]
            
            og_species_list = list(set(og_species_list)) 
            if len(og_species_list)<2 and og_species_list[0] in species_list:

                species = og_species_list[0]
                species = next((species_name for species_name in species_list if species in species_name), species) # replace the above parsed species name with the corresponding name from species_list
                # species_genes = [gene.replace("__", "_")[:-2].split(species+"_")[-1] for gene in orthogroup_members if species in gene] # get gene IDs for each species, replace the instances of "__" with "_" remove the "_1" at the end that orthofinder for some reason adds
                species_genes = []
                for gene in orthogroup_members:
                    # if species in gene:
                    species_genes.append(gene[:-2].replace(f"{species}_", "")) # remove the "_1" that orthofinder adds after each ID, and the species name
                list_spec_genes[species].extend(species_genes)

        # fill the number of species specific genes into a dictionary        
        for species in species_list:
            num_spec_genes[species] = len(list_spec_genes[species])

        if write_to_file:
            gff.write_dict_to_file(list_spec_genes, f"{orthogroups_filepath}/{outfile_name}")

    return(list_spec_genes, num_spec_genes, complete_gene_count_dict)


def check_proteinfasta_single_exon(proteinfasta_filepath, annotation_filepath = "", species = ""):
    fasta_headers_list = []
    transcript_IDs_single_exon = []
    for record in SeqIO.parse(proteinfasta_filepath, "fasta"):
        transcript_id = record.id[:-2]
        transcript_id = transcript_id.split(f"{species}_")[-1] # remove the species name and "_1" suffix from the orthfinder input files
        fasta_headers_list.append(transcript_id)
    # print(fasta_headers_list[:10])
    # print(f"{len(fasta_headers_list)} transcripts in {proteinfasta_filepath}")
    gff_parsed = gff.parse_gff3_general(annotation_filepath, verbose=False)
    # gff_parsed[fasta_headers_list[0]].show
    # print(gff_parsed[fasta_headers_list[0]].features[gff.FeatureCategory.Exon])
    for fasta_header in fasta_headers_list:
        if gff_parsed[fasta_header].category == gff.FeatureCategory.Gene:
            if len(gff_parsed[fasta_header].child_ids_list) == 1:
                transcript_IDs_single_exon.append(fasta_header)

    return(transcript_IDs_single_exon)

    


def get_transcript_number(proteinfasta_filepath):
    fasta_headers_list = []
    for record in SeqIO.parse(proteinfasta_filepath, "fasta"):
        fasta_headers_list.append(record.id)
    return(len(fasta_headers_list))


def read_gene_dicts(infile_path):
    species_dict = {}
    number_of_genes_dict = {}
    with open(infile_path, "r") as infile:
        dict_elements = infile.readlines()
        for species_line in dict_elements:
            species, genes = species_line.split(":")
            gene_id_list = [gene.split(".")[0] for gene in genes.split(",")] # remove the transcript ID at the end if it exists
            species_dict[species]=gene_id_list
            number_of_genes_dict[species]=len(gene_id_list)
    return(species_dict, number_of_genes_dict)


# get a dictionary with species specific genes that are also single-exon
def get_gene_overlap(single_exon_dict, species_specific_dict, verbose = True):

    overlapping_species = [single_exon_species for single_exon_species in list(single_exon_dict.keys()) if single_exon_species in list(species_specific_dict.keys())]
    single_exon_exclusive = [single_exon_species for single_exon_species in list(single_exon_dict.keys()) if single_exon_species not in list(species_specific_dict.keys())]
    species_specific_exclusive = [species_specifc_species for species_specifc_species in list(species_specific_dict.keys()) if species_specifc_species not in list(species_specific_dict.keys())]
    if verbose:
        print(f"overlapping species: {len(overlapping_species)}")
        print(f"single-exon genes species: {len(single_exon_dict.keys())}\t,\t species-specific genes species: {len(species_specific_dict.keys())}")
        if len(single_exon_exclusive)>0:
            print(f"\tspecies only present in the single-exon input: {single_exon_exclusive}")
        if len(species_specific_exclusive):
            print(f"\tspecies only present in the species specifc input: {species_specific_exclusive}")

    # initialize output dictionary that has values that are gene IDs that are species specific and only have one exon
    overlap_dict = {}
    number_of_genes_dict = {} 

    for species in overlapping_species:
        species = next((species_name for species_name in species_names if species in species_name), None) # replace the above parsed species name with the corresponding name from species_names
        overlap = [gene for gene in species_specific_dict[species] if gene in single_exon_dict[species]]
        overlap_dict[species] = overlap
        number_of_genes_dict[species] = len(overlap)
        # break
    
    return(overlap_dict, number_of_genes_dict)




### formatting of gene number y ticks
def format_func(value, tick_number):
    if value >= 1e6:
        return f'{value/1e6:.0f}M'
    elif value >= 1e3:
        return f'{value/1e3:.0f}k'
    else:
        return str(value)



def plot_single_exon_proportion(native_numbers, species_names = [], filename = "", transparent_bg = False):

    print(f" plotting for these {len(species_names)} species: {species_names}")

    # X coordinates for the groups
    x = np.arange(len(species_names))

    # figure proportions according to the data included (longer or shorter)

    # fontsize scales with the dpi somehow which i have to do extra because i change the aspect ratio manually below
    fs = 40 # 37 originally
    
    width = 0.6 # (this is a fraction of the standardized 1 unit of space between axis ticks)
    aspect_ratio = 15 / 10 # nice for a presentation
    ymax_factor = 1.25

    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    ax = fig.add_subplot(111)


    colors = {
        "native" : "#b82946", 
        "orthodb" : "#F2933A", 
        "orthoDB_unmasked" : "#4d7298",
    }
    hatch_color = '#ffffff' # '#E2D4CA' #kind of eggshell white
    plt.rcParams['hatch.color'] = hatch_color


    #### plot native annotation ####

    num_single_exon_genes = native_numbers["num_single_exon_genes"]
    num_total_genes = native_numbers["num_total_genes"]

    single_exon_genes = [num_single_exon_genes[species] for species in species_names]
    print(f"number of single exon genes: {single_exon_genes}")
    multi_exon_genes = [num_total_genes[species]-num_single_exon_genes[species] for species in species_names]
    print(f"number of multi exon genes: {multi_exon_genes}")

    color = [colors["native"], colors["native"]]
    hatching = ["//", "" ]

    category = ""
    x_subtr = 0
    
    print(f"bar width = {width}")

    # total number of genes (with single-exons hatched)
    native_rects1_base = ax.bar(x - x_subtr, single_exon_genes, width, label='proportion of which are single-exon', color= color[0], hatch=hatching[0])
    native_rects1_top = ax.bar(x - x_subtr, multi_exon_genes, width, bottom=single_exon_genes, label='all genes'+category, color= color[1], hatch=hatching[1])
    ymax = max(num_total_genes.values())*ymax_factor
    #### set up labels and stuff ####
    
    ax.set_ylabel('Number of genes', fontsize=fs+4)
    ax.set_title('Proportion of single-exon genes', fontsize=fs+4)
    ax.set_xticks(x)
    ax.set_xlabel('', fontsize=fs+4)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    xtick_labels = [species.replace("_", ". ") for species in species_names]
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=fs)
    # ax.set_yticklabels([f'{int(tick)/1e3:.0f}k' for tick in ax.get_yticks()], fontsize=fs)

    # make custom legend patch for the dashed bars
    plt.rcParams.update({'hatch.color': "#3f3832ff"})
    # dashed_handle = mpatches.Patch(hatch = "//", alpha = 0.0)
    # dashed_label = "proportion of genes that are single-exon"

    # Legend with custom order
    handles, labels = ax.get_legend_handles_labels()
    new_order = [1,0]
    handles = [handles[i] for i in new_order]
    labels = [labels[i] for i in new_order]

    ax.legend(handles, labels, fontsize=fs, ncol=2, loc='upper center')

    # add space at the top of the plot for the legend
    ax.set_ylim(0, int(ymax))
    ax.set_xlim(-0.5, len(xtick_labels)-0.5)

    plt.tight_layout()
    if transparent_bg:
        plt.savefig(filename, dpi = 300, transparent = True)
    else: 
        filename_white = ".".join(filename.split(".")[:-1])
        filename_white = f"{filename_white}_white_bg.png"
        plt.savefig(filename_white, dpi = 300, transparent = False)
    print("Figure saved in the current working directory directory as: "+filename)



if __name__ == '__main__':

    ## calculate all the numbers from input data
    try:
        username = "miltr339"
        tree_path = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder/species_tree.nw"
        species_names = make_species_order_from_tree(tree_path)
    except:
        username = "milena"
        tree_path = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/orthofinder/species_tree.nw"
        species_names = make_species_order_from_tree(tree_path)

    if True:
        ## get the data about single-exon genes by running PhD_chapter3/bash/calculate_single_exon_stats.sh
        ## the resulting output files are evaluated by get_single_exon_genes
        single_exon_stats_dir = f"/Users/{username}/work/chapter3/single_exon_stats"
        single_exon_dict, num_single_exon_dict, num_transcripts_dict = get_single_exon_genes(single_exon_stats_dir, species_names, write_to_file=False, outfile_name="native_single_exon_transcripts_list_14_species.txt", include_total_gene_num = True)

        SE_numbers = {
            "num_single_exon_genes" : num_single_exon_dict,
            "num_total_genes" : num_transcripts_dict
        }
        for key, species_num in SE_numbers.items():
            print(f"{key}:")
            for species, num in species_num.items():
                print(f"\t{species}: {num}")
        # print(f"\t --> gene_numbers = {SE_numbers}")
        print()
        

        data = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/annotation_evaluation/"
        
        ### plot all three annotations
        plot_single_exon_proportion(SE_numbers, species_names = species_names, filename = f"{data}single_exon_genes.png", transparent_bg=True)
        plot_single_exon_proportion(SE_numbers, species_names = species_names, filename = f"{data}single_exon_genes.png", transparent_bg=False)
        
        ### plot native and one other 
        # plot_single_exon_no_species_specific_three_annot(native_numbers, orthoDB_unmasked_numbers = orthodb_unmasked_numbers, species_names = species_names, filename = f"{data}single_exon_Genes_14_species_2_annotations_no_uniform.png", L50_values = L50_values)
        # plot_single_exon_no_species_specific_three_annot(native_numbers, orthodb_numbers = orthodb_numbers, species_names = species_names, filename = f"{data}single_exon_Genes_14_species_2_annotations.png", L50_values = L50_values)