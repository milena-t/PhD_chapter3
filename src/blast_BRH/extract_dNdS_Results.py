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
        return(np.NaN)
    with open(dNdS_path, "r") as dNdS_file:
        lines = dNdS_file.readlines()
        # only look at second line
        value= lines[1].strip().split("\t")[-1]
        return(float(value))


def mean_dNdS_of_pair(results_dir):
    """
    give the dir that contains all the results for the dNdS of the orthologs generated with calculate_pairwise_dNdS.py
    """


if __name__ == "__main__":
    
    chr_type = "A"
    AobtBsil_comp_path = f"A_obtectus_D_sublineata_pairwise_dNdS/"
    results_path = f"/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_{chr_type}/"
    dirpath = f"{results_path}{AobtBsil_comp_path}"
    print(dirpath)
    subdirectories = [f"{dirpath}{d}/2NG.dNdS" for d in os.listdir(dirpath)]
    dNdS_list = [extract_dNdS_file(f) for f in subdirectories]
    print(chr_type)
    print(mean(dNdS_list))

    chr_type = "X"
    results_path = f"/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_{chr_type}/"
    dirpath = f"{results_path}{AobtBsil_comp_path}"
    print(dirpath)
    subdirectories = [f"{dirpath}{d}/2NG.dNdS" for d in os.listdir(dirpath)]
    dNdS_list = [extract_dNdS_file(f) for f in subdirectories]
    print(chr_type)
    print(mean(dNdS_list))