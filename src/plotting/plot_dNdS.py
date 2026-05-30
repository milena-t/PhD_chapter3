"""
plot results from dNdS analysis summaries
"""

import numpy as np
import scipy.stats as sts
from math import sqrt, isnan
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def get_summary_paths(username = "miltr339"):
    summary_paths = {
        "A_dNdS_old" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/dNdS_summary_A-linked.txt",
        "A_dNdS" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/dNdS_by_A_ortholog_pairwise_revisions.txt",
        "A_LRT" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/site_classes_summary_A-linked.txt",
        "X_dNdS_old" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/dNdS_summary_X-linked.txt",
        "X_dNdS" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/dNdS_by_X_ortholog_pairwise_revisions.txt",
        "X_LRT" : f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/site_classes_summary_X-linked.txt",
        "A_LRT_BH_corr" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_A_BH_corrected.txt",
        "X_LRT_BH_corr" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_X_BH_corrected.txt",
    }
    return summary_paths


def read_dNdS_dS_summary_file(summary_path, only_metric="", exclude_list = [], max_dS=10, by_ortholog_file=True):
    """
    read dNdS and dS summary outfiles into a dict by species pair
    out_dict = { pair : {
                    "dS" : [],
                    "dNdS" : []
                    }
                }
    by_ortholog_file=False is the setting for when each line is a species pair and all dN or dNdS values
    if True, then it's the output file from extract dNdS rresults

    max_dS=10 is a very high default that effectively filters nothing
    """
    out_dict = {}
    pairs_done = []
    if by_ortholog_file:
        separator = ":"
        max_dN = max_dS
    else:
        separator = " : "
    with open(summary_path, "r") as summary:
        for line in summary.readlines():
            try:
                pair_ident, dNdS_vals = line.strip().split(separator)
            except:
                raise RuntimeError(f"could not parse line: \n{line[:500]} ...")
            
            if ".log" in pair_ident or len(pair_ident)<5:
                continue
            d_sp = pair_ident.split("_")
            try:
                species1 = f"{d_sp[0]}_{d_sp[1]}"
                species2 = f"{d_sp[2]}_{d_sp[3]}"
            except:
                raise RuntimeError(f"{pair_ident} cannot be parsed as species from '{line[:100]}'")
            if d_sp[-1] != "dS" and d_sp[-1] != "dNdS":
                raise RuntimeError(f"{pair_ident} cannot be parsed as dNdS or dS values")
            
            # skip when species pair contains one excluded species
            if species1 in exclude_list or species2 in exclude_list:
                continue

            pair = f"{species1}_{species2}"
            if pair not in out_dict and only_metric == "":
                out_dict[pair] = {
                    "dS" : [],
                    "dNdS" : []
                }
            elif pair not in out_dict and by_ortholog_file and only_metric!="":
                out_dict[pair] = []

            if by_ortholog_file:
                vals_dict = {"dN" : np.nan, "dS" : np.nan, "dNdS" : np.nan}
                for value in dNdS_vals.split(","):
                    name,number_ = value.split("=")
                    if number_ == "not_found":
                        number = np.nan
                    else:
                        number = float(number_)
                        ## this is set up so that i can also pick different maximum values for dN, dS or dNdS
                        if name == "dN" and number > max_dN:
                            number = np.nan
                        elif name == "dS" and number > max_dS:
                            number = np.nan
                        elif name == "dNdS" and number > max_dS:
                            number = np.nan
                    vals_dict[name] = number

                if np.isnan(vals_dict["dS"]) or np.isnan(vals_dict["dN"]) or np.isnan(vals_dict["dNdS"]):
                    continue
                elif only_metric != "":
                    out_dict[pair].append(vals_dict[only_metric])
                else:
                    out_dict[pair]["dS"].append(vals_dict["dS"])
                    out_dict[pair]["dNdS"].append(vals_dict["dNdS"])
                    
                
            else:
                pairs_done.append(pair)
                values_list = [float(dNdS) if dNdS != 0.0 else np.NaN for dNdS in dNdS_vals.split(",")]
                if d_sp[-1] == "dS":
                    if only_metric=="dS" and max_dS==0:
                        out_dict[pair] = values_list
                    elif only_metric=="dS" and max_dS>0:
                        out_dict[pair] = [val for val in values_list if val < max_dS]
                    else:
                        out_dict[pair]["dS"] = values_list
                if d_sp[-1] == "dNdS" and only_metric == "":
                    out_dict[pair]["dNdS"] = values_list
    
    print(f"read summary file with {len(out_dict)} species pair(s)")
    if True and only_metric!="":
        for key,value in out_dict.items():
            print(f"\t * {key}: {len(value)}")
    elif True:
        for key,value in out_dict.items():
            print(f"\t * {key}:")
            for metric,list in value.items():
                print(f"\t\t- {metric}: {len(list)}")
    return out_dict


def calculate_num_species(dNdS_dict):
    """
    rearrange the binomial coefficient (multiplicative formula, k = 2)
    to get the original number of species from the number of unique unordered pairwise comparisons (excl. self comparison)
    !! this does only work when there's a complete set of comparisons of every species against every other species
    """
    pairs = len(dNdS_dict)
    num_ind = (1+sqrt(1+8*pairs))/2
    return int(num_ind)

def get_species_list(dNdS_dict, exclude_list = []):
    """
    get species list from dNdS dict
    if you add exclude_list, then species in this list will not be added to the total species list
    """
    species = []
    for pair_names in dNdS_dict.keys():
        split_names = pair_names.split("_")
        sp1 = f"{split_names[0]}_{split_names[1]}"
        if sp1 not in exclude_list:
            species.append(sp1)
        sp2 = f"{split_names[2]}_{split_names[3]}"
        if sp2 not in exclude_list:
            species.append(sp2)
    species = sorted(list(set(species)))
    # assert len(species) == calculate_num_species(dNdS_dict)
    # the assertion doesn't hold any more with new species selection because I don't do complete comparisons of every species vs. every other species any more
    return(species)



def permutate_dNdS(dNdS_A, dNdS_X, num_permut = 1000, mean=False):
    """
    permutate n times and calculate the differences of median between all pairs
    """
    medians_diff_list = [np.NaN] * num_permut
    print(f"... running {num_permut} permutations ...")

    def permute_dNdS(dNdS_A, dNdS_X):
        """
        return resampled A and X
        """
        n_A = len(dNdS_A)
        dNdS_all = dNdS_A + dNdS_X
        permut = np.random.permutation(dNdS_all) ## permutation: sampling without replacement
        new_A = permut[0:n_A]
        new_X = permut[n_A:]
        return new_A, new_X

    if mean:
        for i in range(num_permut):
            new_A, new_X = permute_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X)
            m_A = np.nanmedian(new_A)
            m_X = np.nanmedian(new_X)
            medians_diff_list[i] = m_A - m_X
    else: # median
        for i in range(num_permut):
            new_A, new_X = permute_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X)
            m_A = np.nanmean(new_A)
            m_X = np.nanmean(new_X)
            medians_diff_list[i] = m_A - m_X
    
    return medians_diff_list


def calculate_list_CI(values_list:list, cl = 0.95, verbose = False):
    """
    calculate 95% confidence interval of a list of float values
    """
    mean_coeff = np.mean(values_list)
    std_coeff = np.std(values_list)
            
    # Sample statistics
    lower, upper = sts.norm.interval(cl, loc = mean_coeff, scale = std_coeff) 
    norm_coeffs = [mean_coeff,std_coeff, lower, upper]
    if verbose:
        print(f"\t\tmean correlation coefficient: {mean_coeff:.3f}, standard deviation {std_coeff:.3f}, 95% confidence interval: [{lower:.3f}, {upper:.3f}]")
    # ci = sts.t.interval(cl, df=len(values_list)-1, loc=np.mean(values_list), scale=np.std(values_list, ddof=1) / np.sqrt(len(values_list)))
    return(norm_coeffs)



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


def violinplot_pair(data_A_X, row, col, n_A, n_X, mean_A, mean_X, axes, colors_dict,fs, xticks = ["A", "X"], xlab = "", ymax = 0, text_x_offset = 0.4):
    ## make general function so i can repeat it easily for the "mirror" species where row and col are switched
    violins = axes[row,col].violinplot(data_A_X, showmeans = False, showextrema = False)
    colors = [colors_dict["A"], colors_dict["X"]]
    for body, color in zip(violins['bodies'], colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.7)
    
    axes[row, col].set_xlabel('')
    if col-row == 1:
        axes[row, col].set_ylabel('dN/dS', fontsize = fs*0.8)
    elif xlab != "" and col == 0:
        axes[row, col].set_ylabel(xlab, fontsize = fs)
    else:
        axes[row, col].set_ylabel('')
    axes[row, col].tick_params(axis='x', labelsize=fs)
    axes[row, col].tick_params(axis='y', labelsize=fs)
    if ymax == 0:
        max_dNdS_add = 0.3
        axes[row, col].set_ylim([0,1+max_dNdS_add])
        axes[row, col].text(1-text_x_offset, 0.78+max_dNdS_add, f"n={n_A}\nmedian={mean_A:.2f}", fontsize = fs, color = colors_dict["A"])
        axes[row, col].text(2-text_x_offset, 0.78+max_dNdS_add, f"n={n_X}\nmedian={mean_X:.2f}", fontsize = fs, color = colors_dict["X"])
        axes[row, col].hlines(y=mean_A, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["A"])
        axes[row, col].hlines(y=mean_X, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["X"])
        axes[row, col].hlines(y=1, xmin=0.5, xmax=2.5, linewidth=2, linestyle = ":", color="#818181")
    else:
        axes[row, col].set_ylim([0,ymax])
        axes[row, col].text(1-text_x_offset, 0.8*ymax, f"n={n_A}\nmedian={mean_A:.2f}", fontsize = fs, color = colors_dict["A"])
        axes[row, col].text(2-text_x_offset, 0.8*ymax, f"n={n_X}\nmedian={mean_X:.2f}", fontsize = fs, color = colors_dict["X"])
    axes[row, col].set_xticks([1,2])
    axes[row, col].set_xticklabels(xticks)

    return violins


def violinplot_pair_single(data_A_X, col, n_A, n_X, mean_A, mean_X, axes, colors_dict,fs, xticks = ["A", "X"], xlab = "",ylab='dN/dS', ymax = 0,  text_x_offset = 0.4):
    ## make general function so i can repeat it easily for the "mirror" species where row and col are switched
    violins = axes[col].violinplot(data_A_X, showmeans = False, showextrema = False)
    colors = [colors_dict["A"], colors_dict["X"]]
    for body, color in zip(violins['bodies'], colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.7)
    
    axes[col].set_xlabel('')
    axes[col].set_ylabel(ylab, fontsize = fs)
    axes[col].set_xlabel(xlab, fontsize = fs)
    axes[col].tick_params(axis='x', labelsize=fs)
    axes[col].tick_params(axis='y', labelsize=fs*0.85)
    if ymax == 0:
        max_dNdS_add = 0.3
        axes[col].set_ylim([0,1+max_dNdS_add])
        axes[col].text(1- text_x_offset, 0.8+max_dNdS_add, f"n={n_A}\nmedian={mean_A:.2f}", fontsize = fs, color = colors_dict["A"])
        axes[col].text(2- text_x_offset, 0.8+max_dNdS_add, f"n={n_X}\nmedian={mean_X:.2f}", fontsize = fs, color = colors_dict["X"])
        axes[col].hlines(y=1, xmin=0.5, xmax=2.5, linewidth=2, linestyle = ":", color="#818181")
    else:
        axes[col].set_ylim([0,ymax])
        axes[col].text(1-0.5, 2, f"n={n_A}\nmedian={mean_A:.3f}", fontsize = fs*0.9, color = colors_dict["A"])
        axes[col].text(2-0.45, 2, f"n={n_X}\nmedian={mean_X:.3f}", fontsize = fs*0.9, color = colors_dict["X"])
    
    axes[col].hlines(y=mean_A, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["A"])
    axes[col].hlines(y=mean_X, xmin=0.5, xmax=2.5, linewidth=2, color=colors_dict["X"])
    axes[col].set_xticks([1,2])
    axes[col].set_xticklabels(xticks)

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
    # username = "miltr339"
    username = "milena" 
    chromosome = "A"
    data_files = {"A" : ["A_dNdS", "A_LRT"],
                  "X" : ["X_dNdS", "X_LRT"]}
    summary_paths = get_summary_paths(username=username)
    dNdS_dict_A = read_dNdS_dS_summary_file(summary_paths[data_files["A"][0]], max_dS=2,only_metric="dNdS", exclude_list=["D_carinulata", "D_sublineata"])
    dNdS_dict_X = read_dNdS_dS_summary_file(summary_paths[data_files["X"][0]], max_dS=2,only_metric="dNdS", exclude_list=["D_carinulata", "D_sublineata"])
    species = get_species_list(dNdS_dict_A)

    # plot_dNdS_violins(A_dict=dNdS_dict_A, X_dict=dNdS_dict_X,filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/dNdS_violin_plot.png")
    
    if True:
        ## bootstrap significance test to see if dS is sig. diff. on X vs. A
        ## this is not plotted, only command line output
        num_permutations = 10000

        for pair in dNdS_dict_A.keys():
            dNdS_A = dNdS_dict_A[pair]
            dNdS_X = dNdS_dict_X[pair]
            median_diffs = np.nanmedian(dNdS_A) - np.nanmedian(dNdS_X)
            bootstraps = permutate_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X, num_permut=num_permutations)
            mean_cor,std_cor,lower_CI,upper_CI = calculate_list_CI(bootstraps)
            mean_boot = np.mean(bootstraps)
            if median_diffs<lower_CI or median_diffs>upper_CI:
                print(f" *  {pair} --> \t median(dNdS_A)-median(dNdS_X) = {median_diffs:.3f}, mean bootstrap diff = {mean_boot:.5f} with CI [{lower_CI:.5f},{upper_CI:.5f}] --> SIGNIFICANT")
            else:
                print(f" *  {pair} --> \t median(dNdS_A)-median(dNdS_X) = {median_diffs:.3f}, mean bootstrap diff = {mean_boot:.5f} with CI [{lower_CI:.5f},{upper_CI:.5f}] --> (nonsignificant)")


    ### HEATMAP
    if False:
        dNdS_array, species_list = make_means_array_from_dict(dNdS_dict)
        plot_heatmap(counts_array= dNdS_array, species_list=species, filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/mean_dNdS_{chromosome}_linked_heatmap.png", title = f"{chromosome}-linked orthologs mean pairwise dNdS")

    # rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_A/site_classes_summary_A-linked.txt /Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/site_classes_summary_A-linked.txt
    # rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_A_branch_model/dNdS_dS_summary_A-linked_updated_species.txt /Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_A/dNdS_summary_A-linked.txt
    # rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_X/site_classes_summary_X-linked.txt /Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/site_classes_summary_X-linked.txt
    # rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_X_branch_model/dNdS_dS_summary_X-linked_updated_species.txt /Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_X/dNdS_summary_X-linked.txt