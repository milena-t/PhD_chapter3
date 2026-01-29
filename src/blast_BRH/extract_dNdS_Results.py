"""
extract the dNdS values from the results
example filepath:
A_obtectus_B_siliquastri_pairwise_dNdS/A_obtectus_B_siliquastri_A-linked_ortholog_1122_dNdS/2NG.dNdS
make violinplot or boxplot or something
"""

import os
from statistics import mean
import numpy as np


def extract_dNdS_file(dNdS_path, dS_file = False):
    if not os.path.exists(dNdS_path):
        # print(f"{dNdS_path} does not exist")
        # return(np.NaN)
        return(np.nan)
    if os.path.getsize(dNdS_path) == 0:
        print(f"{dNdS_path} has size 0")
        # return(np.NaN)
        return(np.nan)
    with open(dNdS_path, "r") as dNdS_file:
        lines = dNdS_file.readlines()

        if dS_file:
            # only look at third line for dS file
            value= lines[2].strip().split()[-1]
        else:
            # only look at second line for dNdS file
            value= lines[1].strip().split("\t")[-1]

    return(float(value))


def dNdS_list_of_pair(pair_dir, results_dir, only_dNdS = True):
    """
    give the dir that contains all the results for the dNdS of the orthologs generated with calculate_pairwise_dNdS.py
    """
    try:
        os.listdir(f"{results_dir}{pair_dir}")
    except Exception as e:
        # raise RuntimeError(f"{results_dir}{pair_dir} ---> something didn't work! \n{e}")
        print(f"{results_dir}{pair_dir} ---> something didn't work! \n{e}")

    try:
        subdirectories = [f"{results_dir}{pair_dir}/{d}" for d in os.listdir(f"{results_dir}{pair_dir}")]
        subdirectories_dNdS = [f"{d}/2NG.dNdS" for d in subdirectories]
        subdirectories_dS = [f"{d}/2NG.dS" for d in subdirectories]
    except:
        raise RuntimeError(f"pair directory not found! {results_dir}{pair_dir}")
    
    assert len(subdirectories_dS) == len(subdirectories)
    assert len(subdirectories_dNdS) == len(subdirectories)

    if only_dNdS:
        dNdS_list = [extract_dNdS_file(f) for f in subdirectories_dNdS]
        return dNdS_list
    else:
        dNdS_dS_pairs = { subdirectories[i] : [extract_dNdS_file(subdirectories_dNdS[i]), extract_dNdS_file(subdirectories_dS[i], dS_file=True)] for i in range(len(subdirectories)) }
        return dNdS_dS_pairs


def get_dNdS_pairs_dict(results_dir, outfile_name = "", only_dNdS = True):
    """
    Extracts either only dNdS values as a list, or corresponding dS values of the same ortholog as well, in a list of the same order
    """
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
            if ".out" in pair_dir or ".log" in pair_dir or ".txt" in pair_dir:
                continue
            if not os.path.isdir(f"{results_dir}{pair_dir}"):
                raise RuntimeError(f"parsed dir {results_dir}{pair_dir} does not exist!")
            pair_lists[pair_dir] = dNdS_list_of_pair(pair_dir, results_dir, only_dNdS=only_dNdS)
            print(f"{pair_dir} : {pair_lists[pair_dir]}")
        return pair_lists

    else:
        outfile_name = f"{results_dir}{outfile_name}"
        with open(outfile_name, "w") as outfile:
            for pair_dir in pair_lists.keys():
                print(f" --> {pair_dir}")
                if ".out" in pair_dir or ".log" in pair_dir or ".txt" in pair_dir:
                    print(f"\t! log file ignired")
                    continue
                if not os.path.isdir(f"{results_dir}{pair_dir}"):
                    raise RuntimeError(f"parsed dir {results_dir}{pair_dir} does not exist!")
                pair_list = dNdS_list_of_pair(pair_dir, results_dir, only_dNdS=only_dNdS)
                if only_dNdS:
                    pair_list = ",".join([str(dNdS) for dNdS in pair_list])
                    outfile.write(f"{pair_dir} : {pair_list}\n")
                else:
                    pairs_names = list(pair_list.keys())
                    pair_list_dNdS = ",".join([str(pair_list[pair][0]) for pair in pairs_names])
                    pair_list_dS = ",".join([str(pair_list[pair][1]) for pair in pairs_names])
                    pair_name = pair_dir.replace("_pairwise_dNdS", "")
                    outfile.write(f"{pair_name}_dNdS : {pair_list_dNdS}\n")
                    outfile.write(f"{pair_name}_dS : {pair_list_dS}\n")
        print(f"outfile saved to: {outfile_name}\nin {results_dir}")


def extract_site_classes(site_classes_path):
    if not os.path.exists(site_classes_path):
        # print(f"{dNdS_path} does not exist")
        # return(np.NaN)
        return("no_file")
    if os.path.getsize(site_classes_path) == 0:
        print(f"{site_classes_path} has size 0")
        # return(np.NaN)
        return("empty_file")
    with open(site_classes_path, "r") as dNdS_file:
        try:
            p_line, w_line = dNdS_file.readlines()
        except:
            lines = [line.strip() for line in dNdS_file.readlines()]
            lines_str = ";".join(lines)
            return f"parsing_error:{lines_str}"
        
        lines_str = f"{p_line.strip()};{w_line.strip()}"

    return(lines_str)


def site_classes_list_of_pair(pair_dir, results_dir):
    """
    give the dir that contains all the results for the site class tables of the orthologs generated with calculate_pairwise_dNdS.py
    """
    try:
        os.listdir(f"{results_dir}{pair_dir}")
    except Exception as e:
        # raise RuntimeError(f"{results_dir}{pair_dir} ---> something didn't work! \n{e}")
        print(f"{results_dir}{pair_dir} ---> something didn't work! \n{e}")

    try:
        subdirectories = [f"{results_dir}{pair_dir}/{d}" for d in os.listdir(f"{results_dir}{pair_dir}")]
        subdirectories_site_classes = [f"{d}/codeml_M1a_site_class_table.out"  if os.path.exists(f"{d}/codeml_M1a_site_class_table.out") else f"{d}/codeml_M2a_site_class_table.out" for d in subdirectories]
    except:
        raise RuntimeError(f"pair directory not found! {results_dir}{pair_dir}")
    
    # assert len(subdirectories_site_classes) == len(subdirectories)

    site_classes = [f"{f} : {extract_site_classes(f)}" for f in subdirectories_site_classes]
    return site_classes

def get_site_classes(results_dir, outfile_name = ""):
    """
    Extracts either only dNdS values as a list, or corresponding dS values of the same ortholog as well, in a list of the same order
    """
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
            if ".out" in pair_dir or ".log" in pair_dir or ".txt" in pair_dir:
                continue
            if not os.path.isdir(f"{results_dir}{pair_dir}"):
                raise RuntimeError(f"parsed dir {results_dir}{pair_dir} does not exist!")
            pair_lists[pair_dir] = site_classes_list_of_pair(pair_dir, results_dir)
            print(f"{pair_dir} : {pair_lists[pair_dir]}")
        return pair_lists

    else:
        outfile_name = f"{results_dir}{outfile_name}"
        with open(outfile_name, "w") as outfile:
            for pair_dir in pair_lists.keys():
                print(f" --> {pair_dir}")
                if ".out" in pair_dir or ".log" in pair_dir or ".txt" in pair_dir:
                    print(f"\t! log file ignired")
                    continue
                if not os.path.isdir(f"{results_dir}{pair_dir}"):
                    raise RuntimeError(f"parsed dir {results_dir}{pair_dir} does not exist!")
                pair_list = site_classes_list_of_pair(pair_dir, results_dir)
                
                pair_list = "\n".join([str(dNdS) for dNdS in pair_list])
                outfile.write(f"{pair_list}\n") ## !! do not forget fucking tailing newline !!

        print(f"outfile saved to: {outfile_name}\nin {results_dir}")

if __name__ == "__main__":
    
    chr_types = ["X","A"]
    # chr_types = ["A"]
    for chr_type in chr_types:

        # results_path = f"/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_results_{chr_type}/"
        results_path = f"/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_{chr_type}_branch_model/"
        print(chr_type)
        print(f"\n//////////////////// {chr_type} ////////////////////\n")

        get_dNdS_pairs_dict(results_path, outfile_name= f"dNdS_dS_summary_{chr_type}-linked_updated_species.txt", only_dNdS=False)
        # get_site_classes(results_path, outfile_name= f"site_classes_summary_{chr_type}-linked.txt")

#     [f"{dirpath}{d}/2NG.dNdS" for d in os.listdir(results_path)]

# interactive -A uppmax2026-1-8 -t 5:00:00
# module load Biopython/1.86-gfbf-2025b