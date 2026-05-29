import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
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




if __name__ == "__main__":

    username = "miltr339"
    paths = in_paths(username=username)

    plot_correlation(def_file=paths["default"], pair_file=paths["pairwise"], filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/revision_tests/codeml_dNdS_correlation.png")