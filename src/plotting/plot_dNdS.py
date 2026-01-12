"""
plot results from dNdS analysis summaries
"""

import numpy as np
import pandas as pd
from math import sqrt, isnan
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def get_summary_paths(username = "miltr339"):
    summary_paths = {
        "A" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/dNdS_summary_A-linked.txt",
        "X" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/dNdS_summary_X-linked.txt",
    }
    return summary_paths


def read_dNdS_summary_file(summary_path):
    out_dict = {}
    with open(summary_path, "r") as summary:
        for line in summary.readlines():
            try:
                pair_ident, dNdS_vals = line.strip().split(" : ")
            except:
                raise RuntimeError(f"could not parse line: \n{line[:500]} ...")
            pair = pair_ident.split("_pairwise")[0]
            dNdS_list = [float(dNdS) if dNdS != 0.0 else np.NaN for dNdS in dNdS_vals.split(",")]
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

    # initialize array of np.NaN
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


def violinplot_pair(data_A_X, row, col, n_A, n_X, mean_A, mean_X, axes, colors_dict,fs, xticks = ["A", "X"], xlab = ""):
    ## make general function so i can repeat it easily for the "mirror" species where row and col are switched
    violins = axes[row,col].violinplot(data_A_X, showmeans = False, showextrema = False)
    colors = [colors_dict["A"], colors_dict["X"]]
    for body, color in zip(violins['bodies'], colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.7)
    
    max_dNdS_add = 0.3

    axes[row, col].set_xlabel('')
    if col-row == 1:
        axes[row, col].set_ylabel('dN/dS', fontsize = fs*0.8)
    elif xlab != "" and col == 0:
        axes[row, col].set_ylabel(xlab, fontsize = fs)
    else:
        axes[row, col].set_ylabel('')
    axes[row, col].tick_params(axis='x', labelsize=fs*0.8)
    axes[row, col].tick_params(axis='y', labelsize=fs*0.8)
    axes[row, col].set_ylim([0,1+max_dNdS_add])
    axes[row, col].set_xticks([1,2])
    axes[row, col].set_xticklabels(xticks)
    
    axes[row, col].text(1-0.2, 0.78+max_dNdS_add, f"n={n_A}", fontsize = fs, color = colors_dict["A"])
    axes[row, col].text(2-0.2, 0.78+max_dNdS_add, f"n={n_X}", fontsize = fs, color = colors_dict["X"])
    axes[row, col].hlines(y=mean_A, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["A"])
    axes[row, col].hlines(y=mean_X, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["X"])
    axes[row, col].hlines(y=1, xmin=0.5, xmax=2.5, linewidth=2, linestyle = ":", color="#818181")

    return violins

def plot_dNdS_violins(A_dict:dict, X_dict:dict, filename = "dNdS_ratios_A_X.png", legend_in_last = True, dark_mode=False):
    """
    plot a grid of violin plots for all pairwise comparisons
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(dNdS_dict_A)
    
    ### check that A and X are about the same species set
    assert species_list == get_species_list(dNdS_dict_X)
    assert list(A_dict.keys()) == list(X_dict.keys())

    species_count = len(species_list)
    species_index = {species : i for i, species in enumerate(species_list)}

    cols = species_count
    rows = cols
    if rows>2:
        fig, axes = plt.subplots(rows, cols, figsize=(25, 25)) # for more than three rows
    else:
        fig, axes = plt.subplots(rows, cols, figsize=(15, 10)) # for more than three rows
    
    fs = 25

    colors_dict = {
        # "A" : "#4d7298", # uniform_unfiltered blue
        "A" : "#F2933A", # uniform_filtered orange
        "X" : "#b82946", # native red
    }

    diagonals_done = []

    for pair in A_dict.keys():
        ### get pair indices for species pair
        try:
            gen1, spec1, gen2, spec2 =pair.split("_")
        except:
            raise RuntimeError(f"{pair} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        row = species_index[species1]
        species2 = f"{gen2}_{spec2}"
        col = species_index[species2]

        # only do top right matrix
        if row>col:
            col_temp = col
            col = row
            row = col_temp
            species2_temp = species2
            species2 = species1
            species1 = species2_temp

        ## plot species name on diagonals
        if row not in diagonals_done:
            axes[row,row].text(0.8,0.2,f"{species1}", rotation = 90, fontsize = fs)
            diagonals_done.append(row)
            
        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(A_dict[pair], dtype=float)
        data_X_nan = np.array(X_dict[pair], dtype=float)
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            axes[row,col].axis('off')
            axes[row, col].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            # axes[col,row].axis('off')
            # axes[col, row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]
        # plot mirror
        violinplot_pair(data_A_X=data_AX, row=row, col=col, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs=fs)
        # axes[row, col].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
        axes[row, col].set_title(f'{species2}', fontsize = fs)
        # violinplot_pair(data_A_X=data_AX, row=col, col=row, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X)
        # axes[col, row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        axes[col,row].axis('off')
        axes[row,row].axis('off')
        axes[col,col].axis('off')

        print(f"{row}, {col} : {species1} vs. {species2} --> mean dNdS A: {mean_A:.3f}, mean dNdS X: {mean_X:.3f}")
        # if species1 == "D_carinulata" or species2 == "D_carinulata":
        #     print(f"sample sizes, n_A = {n_A}, n_X = {n_X}, data X  = {data_X}")
    
    # fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=fs)
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")
    plt.savefig(filename, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename}")



if __name__ == "__main__":
    username = "miltr339"
    chromosome = "A"
    summary_paths = get_summary_paths
    dNdS_dict_A = read_dNdS_summary_file(summary_paths["A"])
    dNdS_dict_X = read_dNdS_summary_file(summary_paths["X"])
    species = get_species_list(dNdS_dict_A)
    plot_dNdS_violins(A_dict=dNdS_dict_A, X_dict=dNdS_dict_X,filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/dNdS_violin_plot.png")
    

    ### HEATMAP
    if False:
        dNdS_array, species_list = make_means_array_from_dict(dNdS_dict)
        plot_heatmap(counts_array= dNdS_array, species_list=species, filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/mean_dNdS_{chromosome}_linked_heatmap.png", title = f"{chromosome}-linked orthologs mean pairwise dNdS")