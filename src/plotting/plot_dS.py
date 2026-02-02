"""
plot results from dNdS analysis summaries, use dS also as a proxy for baseline mutation rate
"""

import numpy as np
import pandas as pd
from math import sqrt, isnan
import matplotlib.pyplot as plt
import scipy.stats
from plot_dNdS import get_summary_paths,violinplot_pair_single


def read_dNdS_dS_summary_file(summary_path, only_dS, exclude_list = []):
    """
    read dNdS and dS summary outfiles into a dict by species pair
    out_dict = { pair : {
                    "dS" : [],
                    "dNdS" : []
                    }
                }
    """
    out_dict = {}
    pairs_done = []
    with open(summary_path, "r") as summary:
        for line in summary.readlines():
            try:
                pair_ident, dNdS_vals = line.strip().split(" : ")
            except:
                raise RuntimeError(f"could not parse line: \n{line[:500]} ...")
            
            d_sp = pair_ident.split("_")
            try:
                species1 = f"{d_sp[0]}_{d_sp[1]}"
                species2 = f"{d_sp[2]}_{d_sp[3]}"
            except:
                raise RuntimeError(f"{pair_ident} cannot be parsed as species")
            if d_sp[-1] != "dS" and d_sp[-1] != "dNdS":
                raise RuntimeError(f"{pair_ident} cannot be parsed as dNdS or dS values")
            
            # skip when species pair contains one excluded species
            if species1 in exclude_list or species2 in exclude_list:
                continue

            pair = f"{species1}_{species2}"
            if pair not in pairs_done and only_dS == False:
                out_dict[pair] = {
                    "dS" : [],
                    "dNdS" : []
                }
            pairs_done.append(pair)

            values_list = [float(dNdS) if dNdS != 0.0 else np.NaN for dNdS in dNdS_vals.split(",")]
            if d_sp[-1] == "dS":
                if only_dS:
                    out_dict[pair] = values_list
                else:
                    out_dict[pair]["dS"] = values_list
            if d_sp[-1] == "dNdS" and only_dS == False:
                out_dict[pair]["dNdS"] = values_list

    return out_dict

def read_filtered_dNdS_summary(summary_path, excl_list=[], max_dS=2):
    """
    read the dNdS and dS values, filter for min_dS, and return only dNdS values that meet the criteria
    """
    summary_dict = read_dNdS_dS_summary_file(summary_path=summary_path, exclude_list=excl_list, only_dS=False)
    dNdS_dict = {pair : [] for pair in summary_dict.keys()}
    for pair, lists_dict in summary_dict.items():
        dNdS_unfiltered = lists_dict["dNdS"]
        dS_list = lists_dict["dS"]
        assert len(dNdS_unfiltered) == len(dS_list)

        dNdS_filtered = [dNdS for i, dNdS in enumerate(dNdS_unfiltered) if dS_list[i]<max_dS]
        print(f"{pair} : {len(dNdS_unfiltered)} dNdS values, {len(dS_list)} dS values; \t {len(dNdS_filtered)} have dS < {max_dS}")
        dNdS_dict[pair] = dNdS_filtered

    return dNdS_dict


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
        axes[row, col].set_ylabel('dS', fontsize = fs)
    elif xlab != "" and col == 0:
        axes[row, col].set_ylabel(xlab, fontsize = fs)
    else:
        axes[row, col].set_ylabel('')
    axes[row, col].tick_params(axis='x', labelsize=fs)
    axes[row, col].tick_params(axis='y', labelsize=fs)
    # axes[row, col].set_ylim([0,1+max_dNdS_add])
    axes[row, col].set_xticks([1,2])
    axes[row, col].set_xticklabels(xticks)
    
    axes[row, col].text(1-0.2, 1.5+max_dNdS_add, f"n={n_A}", fontsize = fs, color = colors_dict["A"])
    axes[row, col].text(2-0.2, 1.5+max_dNdS_add, f"n={n_X}", fontsize = fs, color = colors_dict["X"])
    axes[row, col].hlines(y=mean_A, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["A"])
    axes[row, col].hlines(y=mean_X, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["X"])
    axes[row, col].hlines(y=1, xmin=0.5, xmax=2.5, linewidth=2, linestyle = ":", color="#818181")

    return violins

def plot_dS_violins(A_dict:dict, X_dict:dict, filename = "dNdS_ratios_A_X.png", legend_in_last = True, dark_mode=False):
    """
    plot a grid of violin plots for all pairwise comparisons
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(A_dict)
    
    ### check that A and X are about the same species set
    assert species_list == get_species_list(A_dict)
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

        # put this here before otherwise the last row/col label never gets reached
        if row == len(species_list)-1:
            species1_lab = species1.replace("_", ". ")
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.2)

        # separate top right matrix from bottom left matrix
        if row>col:
            col_temp = col
            col = row
            row = col_temp
            species2_temp = species2
            species2 = species1
            species1 = species2_temp

        ## plot species name on diagonals
        if row not in diagonals_done:
            species1_lab = species1.replace("_", ". ")
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.2)
            diagonals_done.append(row)
            
        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(A_dict[pair], dtype=float)
        data_X_nan = np.array(X_dict[pair], dtype=float)
        data_A = [dS_A for dS_A in data_A_nan if not np.isnan(dS_A) ]
        data_X = [dS_X for dS_X in data_X_nan if not np.isnan(dS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            axes[row,col].axis('off')
            axes[row,col].text(0.1,0.4,f"{species1}\n{species2}:\nmissing data", fontsize = fs*0.75)
            axes[col,row].axis('off')
            # axes[col, row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        median_A = np.nanmedian(data_A)
        mean_A = np.nanmean(data_A)
        median_X = np.nanmedian(data_X)
        mean_X = np.nanmean(data_X)

        data_AX = [data_A, data_X]
        # plot mirror
        violinplot_pair(data_A_X=data_AX, row=row, col=col, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs=fs)
        # axes[row, col].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
        # axes[row, col].set_title(f'{species2}', fontsize = fs)
        # violinplot_pair(data_A_X=data_AX, row=col, col=row, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X)
        # axes[col, row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        axes[col,row].axis('off')
        axes[row,row].axis('off')
        axes[col,col].axis('off')

        print(f"{row}, {col} : {species1} vs. {species2} --> mean/median dS A: {mean_A:.3f}/{median_A:.3f}, mean/median dS X: {mean_X:.3f}/{median_X:.3f}")
        # if species1 == "D_carinulata" or species2 == "D_carinulata":
        #     print(f"sample sizes, n_A = {n_A}, n_X = {n_X}, data X  = {data_X}")
    
    # fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=fs)
    
    # Adjust layout to prevent overlap  (left, bottom, right, top)
    plt.tight_layout(rect=[0.01, 0, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")
    plt.savefig(filename, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename}")




def calculate_dS_dNdS_lin_reg(dS_list:list, dNdS_list:list, species_pair:str):
    """
    Calculate the linear regression of dS vs. dNdS values of a set of genes
    """
    
    result = scipy.stats.linregress(dS_list, dNdS_list)
    
    ## test normality of residuals
    def predict(x):
        pred_PIC = x*result.slope + result.intercept
        return(pred_PIC)

    residuals = [dNdS_list[i] - predict(dS_list[i]) for i in range(len(dS_list))]
    stat, p_value = scipy.stats.shapiro(residuals)
    
    if p_value < 0.05:
        print(f"\t * non-normal residuals in species pair: {species_pair}")
        # return np.nan, np.nan 
        return result.slope, result.intercept, False
    else:
        return result.slope, result.intercept, True


def make_line_vectors(slope, intercept, x_data, y_data):
    """
    return two vectors to plot slope and intercept of x_data and y_data. the data is needed to correctly set the limits
    """
    low = min([min(x_data), min(y_data)])
    high = max([max(x_data), max(y_data)])
    x = np.linspace(low,high,100)
    y = slope*x + intercept

    return x,y


def plot_dS_vs_dNdS(A_dict:dict, X_dict:dict, filename = "dS_vs_dNdS.png", dark_mode=False, max_dS=2):
    """
    plot a grid of scatterplots of dS vs. dNdS in X and autosomes for each pair
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(A_dict)
    
    ### check that A and X are about the same species set
    assert species_list == get_species_list(A_dict)
    assert list(A_dict.keys()) == list(X_dict.keys())

    species_count = len(species_list)
    species_index = {species : i for i, species in enumerate(species_list)}

    cols = species_count
    rows = cols
    if species_count>4:
        fig, axes = plt.subplots(rows, cols, figsize=(27, 25)) # for more than three rows
        fs = 25
    else:
        fig, axes = plt.subplots(rows, cols, figsize=(20, 18)) # for more than three rows
        fs = 25
    

    colors_dict = {
        # "A" : "#4d7298", # uniform_unfiltered blue
        "A" : "#F2933A", # uniform_filtered orange
        "A_line" : "#D36D0D", # darker orange
        "X" : "#b82946", # native red
        "X_line" : "#861D32" #dark red
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

        # put this here before otherwise the last row/col label never gets reached
        if row == len(species_list)-1:
            species1_lab = species1.replace("_", ". ")
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.2)

        # separate top right matrix from bottom left matrix
        if row>col:
            col_temp = col
            col = row
            row = col_temp
            species2_temp = species2
            species2 = species1
            species1 = species2_temp

        ## plot species name on diagonals
        if row not in diagonals_done:
            species1_lab = species1.replace("_", ". ")
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.2)
            diagonals_done.append(row)
            
        # Extract dS and dNdS
        ## exclude all the NaNs because violinplot can't handle them
        dS_A_nan = np.array(A_dict[pair]["dS"], dtype=float)
        dS_X_nan = np.array(X_dict[pair]["dS"], dtype=float)
        dNdS_A_nan = np.array(A_dict[pair]["dNdS"], dtype=float)
        dNdS_X_nan = np.array(X_dict[pair]["dNdS"], dtype=float)

        # remove nan's from both (remove both elements if one of them is nan)
        dS_A_all = [dS_A for i, dS_A in enumerate(dS_A_nan) if not np.isnan(dS_A_nan[i]) and not np.isnan(dNdS_A_nan[i]) ]
        dS_X_all = [dS_X for i, dS_X in enumerate(dS_X_nan) if not np.isnan(dS_X_nan[i]) and not np.isnan(dNdS_X_nan[i]) ]
        dNdS_A_all = [dNdS_A for i, dNdS_A in enumerate(dNdS_A_nan) if not np.isnan(dS_A_nan[i]) and not np.isnan(dNdS_A_nan[i]) ]
        dNdS_X_all = [dNdS_X for i, dNdS_X in enumerate(dNdS_X_nan) if not np.isnan(dS_X_nan[i]) and not np.isnan(dNdS_X_nan[i]) ]
        assert len(dS_A_all) == len(dNdS_A_all)
        assert len(dS_X_all) == len(dNdS_X_all)
        # remove dS>mad_dS from both (remove dNdS also if dS is removed)
        dS_A = [dS for dS in dS_A_all if dS < max_dS]
        dS_X = [dS for dS in dS_X_all if dS < max_dS]
        dNdS_A = [dNdS for i, dNdS in enumerate(dNdS_A_all) if dS_A_all[i] < max_dS]
        dNdS_X = [dNdS for i, dNdS in enumerate(dNdS_X_all) if dS_X_all[i] < max_dS]
        assert len(dS_A) == len(dNdS_A)
        assert len(dS_X) == len(dNdS_X)

        print(f"\tremoved dS>{max_dS}:\n\t{len(dNdS_A_all)} A dNdS values (no np.nan) before filtering, {len(dNdS_A)} after filtering, \n\t{len(dNdS_X_all)} X dNdS values (no np.nan) before filtering , {len(dNdS_X)} after filtering")

        if len(dS_A)==0 or len(dS_X)==0 or len(dNdS_A)==0 or len(dNdS_X)==0:
            axes[row,col].axis('off')
            axes[row,col].text(0.1,0.4,f"{species1}\n{species2}:\nmissing data", fontsize = fs*0.75)
            axes[col,row].axis('off')
            axes[col,row].text(0.1,0.4,f"{species1}\n{species2}:\nmissing data", fontsize = fs*0.75)
            continue

        n_A = len(dS_A)
        n_X = len(dS_X)
        mean_A = np.nanmedian(dS_A)
        mean_X = np.nanmedian(dS_X)

        dS_AX = [dS_A, dS_X]

        # plot dS violins
        violinplot_pair(data_A_X=dS_AX, row=row, col=col, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs=fs)

        # plot dNdS scatters
        axes[col,row].scatter(dS_A, dNdS_A, color = colors_dict["A"], s=35)
        axes[col,row].scatter(dS_X, dNdS_X, color = colors_dict["X"], s=35)
        axes[col,row].tick_params(axis='x', labelsize=fs)
        axes[col,row].tick_params(axis='y', labelsize=fs)
        axes[species_count-1,row].set_xlabel("dS", fontsize = fs)
        axes[col,0].set_ylabel("dNdS", fontsize = fs)
        axes[col,row].set_xlim(-0.08,2.08)
        ## linear regression
        slope_A, intercept_A, normal_residuals_A = calculate_dS_dNdS_lin_reg(dS_list = dS_A, dNdS_list = dNdS_A, species_pair= pair)
        slope_X, intercept_X, normal_residuals_X = calculate_dS_dNdS_lin_reg(dS_list = dS_X, dNdS_list = dNdS_X, species_pair= pair)
        linreg_x_A, linreg_y_A = make_line_vectors(slope=slope_A, intercept=intercept_A, x_data=dS_A, y_data=dNdS_A)
        linreg_x_X, linreg_y_X = make_line_vectors(slope=slope_X, intercept=intercept_X, x_data=dS_X, y_data=dNdS_X)

        if normal_residuals_A:
            linestyle_A = ":"
        else:
            linestyle_A = "-"
        if normal_residuals_X:
            linestyle_X = ":"
        else:
            linestyle_X = "-"

        axes[col,row].plot(linreg_x_A, linreg_y_A, color = colors_dict["A_line"], linewidth=2, label=f"slope: {slope_A:.3f}", linestyle=linestyle_A)    
        axes[col,row].plot(linreg_x_X, linreg_y_X, color = colors_dict["X_line"], linewidth=2, label=f"slope: {slope_X:.3f}", linestyle=linestyle_X)
        axes[col,row].legend(fontsize=fs*0.75)

        axes[row,row].axis('off')
        axes[col,col].axis('off')

        print(f"{row}, {col} : {species1} vs. {species2}, slope A: {slope_A:.2f}, slope X: {slope_X:.2f}")
        # if species1 == "D_carinulata" or species2 == "D_carinulata":
        #     print(f"sample sizes, n_A = {n_A}, n_X = {n_X}, data X  = {data_X}")
    
    # fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=fs)
    # Adjust layout to prevent overlap  (left, bottom, right, top)
    plt.tight_layout(rect=[0.01, 0, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")
    # transparent background
    plt.savefig(filename, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = filename.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename} and {filename_tr}")



if __name__ == "__main__":
    username = "miltr339"
    chromosome = "A"
    data_files = {"A" : ["A_dNdS", "A_LRT"],
                  "X" : ["X_dNdS", "X_LRT"]}
    summary_paths = get_summary_paths(username=username)

    # bruchini
    if False:
        species_excl = ["D_carinulata", "D_sublineata", "T_castaneum", "T_freemani", "C_septempunctata", "C_magnifica"]
        filename =f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/dS_vs_dNdS_scatterplot_bruchini.png"
    # coccinella
    elif False:
        species_excl = ["D_carinulata", "D_sublineata", "T_castaneum", "T_freemani", "B_siliquastri", "A_obtectus", "C_maculatus", "C_chinensis"]
        filename =f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/dS_vs_dNdS_scatterplot_coccinella.png"
    # tribolium
    elif True:
        species_excl = ["D_carinulata", "D_sublineata", "C_septempunctata", "C_magnifica", "B_siliquastri", "A_obtectus", "C_maculatus", "C_chinensis"]
        filename =f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/dS_vs_dNdS_scatterplot_tribolium.png"
    
    
    ## make the violinplots of only dS for each pair
    if False:
        print(f"/////////////// A ///////////////")
        dS_dict_A = read_dNdS_summary_file(summary_paths["A"], only_dS = True)
        # print(dS_dict_A)

        print(f"/////////////// X ///////////////")
        dS_dict_X = read_dNdS_summary_file(summary_paths["X"], only_dS = True)
        # print(dS_dict_X)
        
        species = get_species_list(dS_dict_A)
        plot_dS_violins(A_dict=dS_dict_A, X_dict=dS_dict_X,filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/dS_violin_plot.png")
    
    ## make the scatterplots
    if True:
        print(f"reading A ...")
        dS_dict_A = read_dNdS_dS_summary_file(summary_paths[data_files["A"][0]], only_dS = False, exclude_list=species_excl)
        # print(dS_dict_A)

        print(f"reading X ...")
        dS_dict_X = read_dNdS_dS_summary_file(summary_paths[data_files["X"][0]], only_dS = False, exclude_list=species_excl)
        # print(dS_dict_X)
        
        species = get_species_list(dS_dict_A)
        plot_dS_vs_dNdS(A_dict=dS_dict_A, X_dict=dS_dict_X,filename=filename, max_dS=2)
    
