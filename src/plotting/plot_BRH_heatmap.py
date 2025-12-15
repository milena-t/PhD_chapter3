"""
plot a heatmap from the results of identifying pairwise best reciprocal hits from blast searches
count either the number of X-linked or A-linked ortholog pairs, or if it's a self blast count the number of gametologs
"""

# import polars as pl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

def ortholog_tables(username = "miltr339"):
    dir_path = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/"
    dirs_list = [
        f"{dir_path}A_obtectus_A_obtectus_BRH.tsv",
        f"{dir_path}A_obtectus_B_siliquastri_BRH.tsv",
        f"{dir_path}A_obtectus_C_chinensis_BRH.tsv",
        f"{dir_path}A_obtectus_C_maculatus_BRH.tsv",
        f"{dir_path}A_obtectus_D_carinulata_BRH.tsv",
        f"{dir_path}A_obtectus_D_sublineata_BRH.tsv",
        f"{dir_path}B_siliquastri_B_siliquastri_BRH.tsv",
        f"{dir_path}B_siliquastri_C_chinensis_BRH.tsv",
        f"{dir_path}B_siliquastri_C_maculatus_BRH.tsv",
        f"{dir_path}B_siliquastri_D_carinulata_BRH.tsv",
        f"{dir_path}B_siliquastri_D_sublineata_BRH.tsv",
        f"{dir_path}C_chinensis_C_chinensis_BRH.tsv",
        f"{dir_path}C_chinensis_C_maculatus_BRH.tsv",
        f"{dir_path}C_chinensis_D_carinulata_BRH.tsv",
        f"{dir_path}C_chinensis_D_sublineata_BRH.tsv",
        f"{dir_path}C_maculatus_C_maculatus_BRH.tsv",
        f"{dir_path}C_maculatus_D_carinulata_BRH.tsv",
        f"{dir_path}C_maculatus_D_sublineata_BRH.tsv",
        f"{dir_path}D_carinulata_D_carinulata_BRH.tsv",
        f"{dir_path}D_carinulata_D_sublineata_BRH.tsv",
        f"{dir_path}D_sublineata_D_sublineata_BRH.tsv",
    ]
    return dirs_list, dir_path

def count_orthologs(brh_table_path:str, chr_type:str = "A", verbose = False) -> int:
    """
    count the number of brh pairs between two species in a brh table created with src/make_blastBRH_all_species.py
    You can choose if you count exclusively A-linked, X-linked or Y-linked pairs. pairs that are on different chromosome types cannot be counted
    """
    # brh_table = pl.read_csv(brh_table_path, separator="\t")
    brh_table = pd.read_csv(brh_table_path, sep="\t")
    brh_filtered = brh_table[brh_table["chromosome"] == chr_type]
    brh_filtered = brh_filtered[brh_filtered["chromosome.1"] == chr_type]
    rows_after = brh_filtered.shape[0]
    if verbose:
        rows_before = brh_table.shape[0]
        # print(brh_filtered.head())
        brh_filename = brh_table_path.split("/")[-1]
        print(f"--> {brh_filename}")
        print(f"all BRH hits: {rows_before}\nexclusively {chr_type}-linked: {rows_after}")
    return rows_after




def make_array_for_heatmap(brh_tables_list:list, chr_type:str = "A", verbose = False):
    """
    make an array to plot as a heatmap
    """

    ## initialize array of correct size, which is a square of the number of species
    species_index = {}
    species_count = 0
    for brh_table in brh_tables_list:
        brh_table = brh_table.split("/")[-1]
        brh_table = brh_table.replace("_BRH.tsv", "")
        try:
            gen1, spec1, gen2, spec2 =brh_table.split("_")
        except:
            raise RuntimeError(f"{brh_table} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        species2 = f"{gen2}_{spec2}"
        if species1 == species2:
            species_index[species1] = species_count
            species_count += 1
    species_count = len(species_index) 
    if verbose:
        print(f"{species_count} species: {species_index}")
    
    ortholog_counts = np.zeros((species_count,species_count), dtype=float)

    ## fill the initalized table with the counts
    for brh_table_path in brh_tables_list:
        brh_table = brh_table_path.split("/")[-1]
        brh_table = brh_table.replace("_BRH.tsv", "")
        try:
            gen1, spec1, gen2, spec2 =brh_table.split("_")
        except:
            raise RuntimeError(f"{brh_table} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        index1 = species_index[species1]
        species2 = f"{gen2}_{spec2}"
        index2 = species_index[species2]
        if species1 == species2:
            ortholog_counts[index1, index2] = np.NaN
        else:
            count_brh = count_orthologs(brh_table_path, chr_type=chr_type, verbose=False)
            ortholog_counts[index1, index2] = count_brh
            ortholog_counts[index2, index1] = count_brh

    if verbose:
        print(ortholog_counts)

    species_list = list(species_index.keys())
    return ortholog_counts, species_list


def plot_heatmap(counts_array, species_list, filename = "BRH_orthologs_heatmap.png", title = f"pairwise orthologs counts"):
    """
    plot the heatmap created in make_array_for_heatmap()
    """

    fs = 40
    aspect_ratio = 15 / 15
    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    ax = fig.add_subplot(111)

    # cmap=mpl.colormaps["rainbow"]
    cmap = LinearSegmentedColormap.from_list("red_to_orange",["#b82946", "#F2933A"])
    # cbarlabel="number of pairwise 1-to-1 ortholgs"
    im = ax.imshow(counts_array, cmap=cmap)

    # create text annotations
    for i in range(len(species_list)):
        for j in range(len(species_list)):
            try:
                count = int(counts_array[i, j])
                text = ax.text(j, i, count,ha="center", va="center", color="w", fontsize = fs)
            except:
                continue

    species_list = [species.replace("_", ". ") for species in species_list]
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(len(species_list)), labels=species_list,rotation=45, ha="right", rotation_mode="anchor", fontsize = fs)
    ax.set_yticks(range(len(species_list)), labels=species_list, fontsize = fs)
    plt.title(label=title, fontsize = fs*1.3)
    
    plt.tight_layout()
    # plt.show()
    plt.savefig(filename, dpi = 300, transparent = False)
    print(f"figure saved here: {filename}")



def count_gametologs(brh_table_path:str, verbose = False) -> int:
    """
    count the number of brh pairs between X and Y chromosome in a brh table created with src/make_blastBRH_all_species.py
    the BRH table should be from a species blasted against itself! everything else doesn't make sense
    """
    # brh_table = pl.read_csv(brh_table_path, separator="\t")
    brh_table = pd.read_csv(brh_table_path, sep="\t")
    # I only need to check one direction, because of the self blast, the other direction would be identical.
    # so no need to make a table that filters in column 1 for x and in column 3 for y, and then the inverse
    # where you filter in column 1 for y and in 3 or x, because those would be the same
    brh_filtered = brh_table[brh_table["chromosome"] == "X"]
    brh_filtered = brh_filtered[brh_filtered["chromosome.1"] == "Y"]
    rows_after = brh_filtered.shape[0]
    if verbose:
        rows_before = brh_table.shape[0]
        # print(brh_filtered.head())
        brh_filename = brh_table_path.split("/")[-1]
        print(f"--> {brh_filename}")
        print(f"all BRH hits: {rows_before}\nX-Y gametologs: {rows_after}")
    return rows_after


def count_all_gametologs(brh_tables_list:list, verbose = False)->dict:
    """
    return a dictionary with the number of X-Y gametologs in all species
    """
    gam_counts = {}
    for brh_table_path in brh_tables_list:
        brh_table = brh_table_path.split("/")[-1]
        brh_table = brh_table.replace("_BRH.tsv", "")
        try:
            gen1, spec1, gen2, spec2 =brh_table.split("_")
        except:
            raise RuntimeError(f"{brh_table} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        species2 = f"{gen2}_{spec2}"
        if species1 == species2:
            gam_counts[species1] = count_gametologs(brh_table_path)
            if verbose:
                print(f"{species1} : {gam_counts[species1]}")
    return gam_counts

if __name__ == "__main__":
    username = "miltr339"
    all_tables_list, dir_path = ortholog_tables(username=username)

    # count_A = count_orthologs(all_tables_list[1], chr_type="A", verbose=True)
    # count_X = count_orthologs(all_tables_list[1], chr_type="X", verbose=True)
    # count_g = count_gametologs(all_tables_list[0], verbose=True)

    if True:
        for chromosome in ["A", "X"]:
            ortholog_counts_array,species_list = make_array_for_heatmap(all_tables_list, chr_type=chromosome, verbose=False)
            plot_heatmap(ortholog_counts_array,species_list, filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/BRH_{chromosome}_linked_counts_heatmap.png", title = f"{chromosome}-linked pairwise orthologs counts")
    
    ## count gametologs
    print(f"--> count within-species gametologs")
    count_all_gametologs(all_tables_list, verbose=True)