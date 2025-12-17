"""
extract the dNdS values from the results
example filepath:
A_obtectus_B_siliquastri_pairwise_dNdS/A_obtectus_B_siliquastri_A-linked_ortholog_1122_dNdS/2NG.dNdS
make violinplot or boxplot or something
"""

import os
from statistics import mean
import numpy as np


def extract_dNdS_file(dNdS_path):
    if not os.path.exists(dNdS_path):
        # print(f"{dNdS_path} does not exist")
        return(np.NaN)
    if os.path.getsize(dNdS_path) == 0:
        print(f"{dNdS_path} has size 0")
        return(np.NaN)
    with open(dNdS_path, "r") as dNdS_file:
        lines = dNdS_file.readlines()
        # only look at second line
        value= lines[1].strip().split("\t")[-1]
        return(float(value))


def dNdS_list_of_pair(pair_dir, results_dir):
    """
    give the dir that contains all the results for the dNdS of the orthologs generated with calculate_pairwise_dNdS.py
    """
    try:
        subdirectories = [f"{results_dir}{pair_dir}/{d}/2NG.dNdS" for d in os.listdir(f"{results_dir}{pair_dir}")]
    except:
        ### TODO the error gi
        raise RuntimeError(f"pair directory not found! {results_dir}{pair_dir}")
    dNdS_list = [extract_dNdS_file(f) for f in subdirectories]
    return dNdS_list

def get_dNdS_pairs_dict(results_dir, outfile_name = ""):
    pair_dirs = []
    for d in os.listdir(results_dir):
        if os.path.isfile(d):
            continue
        d_sp = d.split("_")
        try:
            species1 = f"{d_sp[0]}_{d_sp[1]}"
            species2 = f"{d_sp[2]}_{d_sp[3]}"
        except:
            print(f"{d} cannot be parsed as species")
            continue
        if species1 != species2:
            pair_dirs.append(d)

    pair_lists = {f"{d}":[] for d in pair_dirs}
    if outfile_name == "":
        for pair_dir in pair_lists.keys():
            if not os.path.isdir(f"{results_dir}{pair_dir}"):
                raise RuntimeError(f"parsed dir {results_dir}{pair_dir} does not exist!")
            pair_lists[pair_dir] = dNdS_list_of_pair(pair_dir, results_dir)
            print(f"{pair_dir} : {pair_lists[pair_dir]}")
        return pair_lists

    else:
        outfile_name = f"{results_dir}{outfile_name}"
        with open(outfile_name, "w") as outfile:
            for pair_dir in pair_lists.keys():
                pair_list = dNdS_list_of_pair(pair_dir, results_dir)
                pair_list = ",".join([str(dNdS) for dNdS in pair_list])
                outfile.write(f"{pair_dir} : {pair_list}\n")
        print(f"outfile saved to: {outfile_name}\nin {results_dir}")


if __name__ == "__main__":
    
    chr_types = ["A","X"]
    for chr_type in chr_types:
        # results_path = f"/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_{chr_type}/"
        results_path = f"/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_{chr_type}/"
        print(f"\n//////////////////// {chr_type} ////////////////////\n")
        get_dNdS_pairs_dict(results_path, f"dNdS_summary_{chr_type}-linked.txt")

#     [f"{dirpath}{d}/2NG.dNdS" for d in os.listdir(results_path)]

