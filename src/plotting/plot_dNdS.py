"""
plot results from dNdS analysis summaries
"""

import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

username = "miltr339"
summary_paths = {
    "A" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/dNdS_summary_A-linked.txt",
    "X" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/dNdS_summary_X-linked.txt",
}


def read_dNdS_summary_file(summary_path):
    out_dict = {}
    with open(summary_path, "r") as summary:
        for line in summary.readlines():
            try:
                pair_ident, dNdS_vals = line.strip().split(" : ")
            except:
                raise RuntimeError(f"could not parse line: \n{line[:500]} ...")
            pair = pair_ident.split("_pairwise")[0]
            dNdS_list = [float(dNdS) for dNdS in dNdS_vals.split(",")]
            out_dict[pair] = dNdS_list

    return out_dict


def calculate_num_species(dNdS_dict):
    """
    rearrange the binomial coefficient (multiplicative formula, k = 2)
    to get the original number of species from the number of unique unordered pairwise comparisons (excl. self comparison)
    
    """
    pairs = len(dNdS_dict)
    num_ind = (1+sqrt(1+8*pairs))/2
    return int(num_ind)

def get_species_list(dNdS_dict):
    """
    get species list from dNdS dict
    """
    species = []
    for pair_names in dNdS_dict.keys():
        split_names = pair_names.split("_")
        sp1 = f"{split_names[0]}_{split_names[1]}"
        species.append(sp1)
        sp2 = f"{split_names[2]}_{split_names[3]}"
        species.append(sp2)
    species = list(set(species))
    assert len(species) == calculate_num_species(dNdS_dict)
    return(species)


def make_means_array_from_dict(dNdS_dict, verbose = True):
    """
    make an array to plot as a heatmap
    """
    species_list = get_species_list(dNdS_dict)
    species_count = len(species_list)
    species_index = {species : i for i, species in enumerate(species_list)}


    pairwise_dNdS = np.full((species_count,species_count), np.NaN)

    ## fill the initalized table with the counts
    for pair, dNdS_list in dNdS_dict.items():
        try:
            gen1, spec1, gen2, spec2 =pair.split("_")
        except:
            raise RuntimeError(f"{pair} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        index1 = species_index[species1]
        species2 = f"{gen2}_{spec2}"
        index2 = species_index[species2]
        if species1 == species2:
            pairwise_dNdS[index1, index2] = np.NaN
        else:
            if not np.isnan(dNdS_list).all():
                mean_dNdS = np.nanmean(dNdS_list)
            else:
                mean_dNdS = np.NaN
            pairwise_dNdS[index1, index2] = mean_dNdS
            pairwise_dNdS[index2, index1] = mean_dNdS

    if verbose:
        print(pairwise_dNdS)
        print(type(pairwise_dNdS))

    return pairwise_dNdS, species_list


def plot_heatmap(counts_array, species_list, filename = "mean_dNdS.png", title = f"mean dNdS"):
    """
    plot the heatmap created in make_means_array_from_dict()
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
                # count = counts_array[i, j]
                count = f"{counts_array[i, j]:.3}"
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


if __name__ == "__main__":
    username = "miltr339"
    chromosome = "A"

    dNdS_dict = read_dNdS_summary_file(summary_paths[chromosome])
    species = get_species_list(dNdS_dict)
    dNdS_array, species_list = make_means_array_from_dict(dNdS_dict)
    plot_heatmap(counts_array= dNdS_array, species_list=species, filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/mean_dNdS_{chromosome}_linked_heatmap.png", title = f"{chromosome}-linked orthologs mean pairwise dNdS")