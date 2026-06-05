import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from bootstrap_dNdS import permutate_dNdS,calculate_list_CI
import numpy as np
import pandas as pd

def in_paths(username="miltr339"):
    files_dir=f"/Users/{username}/work/chapter3/revision_tests"
    out_dict = {
        "default" : f"{files_dir}/dNdS_by_ortholog_default_revisions.txt",
        "pairwise" : f"{files_dir}/dNdS_by_ortholog_pairwise_revisions.txt",
    }
    return out_dict


def read_dNdS_dS_summary_file(summary_path, exclude_list = [], max_dS=2, max_dN=2):
    """
    read dNdS and dS summary outfiles into a dict by species pair
    out_dict = { orthologID : {
                    "dS" : float,
                    "dN" : float,
                    "dNdS" : float,
                    }
                }
    """
    out_dict = {}
    pairs_done = []
    with open(summary_path, "r") as summary:
        for line in summary.readlines():
            try:
                ortholog_ID, dNdS_vals = line.strip().split(":")
            except:
                raise RuntimeError(f"could not parse line: \n{line[:500]} ...")

            out_dict[ortholog_ID] = {}
            for value in dNdS_vals.split(","):
                name,number_ = value.split("=")

                if number_ == "not_found":
                    number = np.nan
                else:
                    number = float(number_)
                    if name == "dN" and number > max_dN:
                        number = np.nan
                    elif name == "dS" and number > max_dS:
                        number = np.nan
 
                out_dict[ortholog_ID][name] = number

    return out_dict


def plot_correlation(def_file, pair_file, filename = ""):
    """
    plot correlation between two summary files
    """
    def_dict = read_dNdS_dS_summary_file(def_file)
    pair_dict = read_dNdS_dS_summary_file(pair_file)
    # print(pair_dict)

    orthologs_order = list(def_dict.keys())
    def_dN_list = []
    def_dS_list = []
    def_dNdS_list = []
    pair_dN_list = []
    pair_dS_list = []
    pair_dNdS_list = []

        
    for ortholog in orthologs_order:
        def_vals = def_dict[ortholog]
        pair_vals = pair_dict[ortholog]
        def_dN_list.append(def_vals["dN"])
        pair_dN_list.append(pair_vals["dN"])
        def_dS_list.append(def_vals["dS"])
        pair_dS_list.append(pair_vals["dS"])
        def_dNdS_list.append(def_vals["dN"]/def_vals["dS"])
        pair_dNdS_list.append(pair_vals["dN"]/pair_vals["dS"])
    
    pair_lists = {
        "dN" : pair_dN_list,
        "dS" : pair_dS_list,
        "dNdS" : pair_dNdS_list,
    }
    def_lists = {
        "dN" : def_dN_list,
        "dS" : def_dS_list,
        "dNdS" : def_dNdS_list,
    }

    fs = 30
    aspect_ratio = 15 / 15
    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig, ax = plt.subplots(1, 3, figsize=(30, 10))
    
    for i,metric in enumerate(["dN", "dS", "dNdS"]):
        ax[i].plot([0, 1], [0, 1], transform=ax[i].transAxes, color = "black", linestyle="dashed", linewidth=4)
        ax[i].scatter(pair_lists[metric], def_lists[metric], color="#b82946", s=fs*1.5)
        ax[i].tick_params(axis='x', labelsize=fs)
        ax[i].tick_params(axis='y', labelsize=fs)
        ax[i].set_xlabel("pairwise", fontsize = fs)
        ax[i].set_ylabel("default", fontsize = fs)

        # pandas automatically removes the orthologs where at least one value is nan
        corr = pd.Series(pair_lists[metric]).corr(pd.Series(def_lists[metric]), method="pearson")
        print(f"{metric}")
        print(corr)
        print()
        ax[i].set_title(f"{metric} values (pearson R = {corr:.4f})", fontsize = fs)

    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {filename}")



def bootstrap_pos_sel(X_data:dict, A_data:dict, num_permutations=1000):
    """
    calcualte bootstrap confidence intervals of enrichment of pos.sel. genes on X or A.
    Input data are two dictionaries that contain count numbers of positively selected genes and not pos. sel. genes
    { pos_sel : int , non_pos_sel : int }
    """
    
    X_vec = [1]*X_data["pos_sel"] + [0]*X_data["non_pos_sel"]
    A_vec = [1]*A_data["pos_sel"] + [0]*A_data["non_pos_sel"]
    X_mean = np.nanmean(X_vec)
    A_mean = np.nanmean(A_vec)
    mean_diff = A_mean - X_mean

    print(f"X_mean prop: {X_mean:.5} ({len(X_vec)} samples)")
    print(f"A_mean prop: {A_mean:.5} ({len(A_vec)} samples)")
    
    bootstraps = permutate_dNdS(dNdS_A=A_vec, dNdS_X=X_vec, num_permut=num_permutations)
    mean_cor,std_cor,lower_CI,upper_CI = calculate_list_CI(bootstraps)
    mean_boot = np.mean(bootstraps)
    if mean_diff<lower_CI or mean_diff>upper_CI:
        print(f" * mean(pos_sel_A)-mean(pos_sel_X) --> \t{mean_diff:.5f}, mean bootstrap diff = {mean_boot:.5f} with CI [{lower_CI:.5f},{upper_CI:.5f}] --> SIGNIFICANT")
    else:
        print(f" * mean(pos_sel_A)-mean(pos_sel_X) --> \t{mean_diff:.5f}, mean bootstrap diff = {mean_boot:.5f} with CI [{lower_CI:.5f},{upper_CI:.5f}] --> (nonsignificant)")


if __name__ == "__main__":

    username = "miltr339"
    paths = in_paths(username=username)

    if False:
        # plot the correlation of codeml dN and dS values for pairwise comparisons of runmode=1 and runmode=-2
        plot_correlation(def_file=paths["default"], pair_file=paths["pairwise"], filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/revision_tests/codeml_dNdS_correlation.png")
    if False:
        # do permutation test on the 4-way 1-to-1 orthologs in Bruchini
        bootstrap_pos_sel(X_data={"pos_sel" : 6 , "non_pos_sel" : 268}, A_data={"pos_sel" : 60 , "non_pos_sel" : 1881}, num_permutations=10000)
    
    if True:
        # Compare positively selected orthologs in the pairwise tests to if they are also positively selected in the 4-way comparison
        # TODO
        pass