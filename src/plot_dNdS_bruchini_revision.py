"""
Plot all the branch-wise dNdS estimates from the branch model of the Bruchini phylogeny for the revision
"""
import pandas as pd
import numpy as np
from tqdm import tqdm
from plot_dNdS import permutate_dNdS, calculate_list_CI

def filepaths(username="miltr339"):
    data_dir_X = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_seqs_revision_X/"
    data_dir_A = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_seqs_revision_A/"
    files_dict = {
        "branch_model_X" : f"{data_dir_X}run_branch_model.log",
        "site_model_X" : f"{data_dir_X}run_site_model.log",
        "site_model_beta_X" : f"{data_dir_X}run_site_model_beta.log",
        "branch_model_A" : f"{data_dir_A}branch_model_bruchini_res_A.log",
        "site_model_beta_A" : f"{data_dir_A}site_model_beta_bruchini_res_A.log",
        "site_model_A" : f"{data_dir_A}site_model_bruchini_res_A.log",
    }
    return files_dict


def read_branch_log_outfile(filename:str, species_association_dict:dict):
    """
    The orthologs log file comes with a matrix that shows all the dNdS values for all pairwise comparisons 
    I the names are the gene IDs and I hardcode how the naming sceme associates back to the species names 
    """

    out_dict = {}
    curr_OG_dNdS=False
    ortholog = ""
    with open(filename, "r") as log_file:
        all_log_lines = log_file.readlines()
        for i, log_line in tqdm(enumerate(all_log_lines)):
            line = log_line.strip()
            if line == "":
                continue
            if "====================== Bruchini" in line:
                curr_ortholog = line.split(" ")[1]
                if curr_ortholog != ortholog:
                    ortholog=curr_ortholog
                    curr_OG_dNdS=True

            if "dN/dS values for " in line and curr_OG_dNdS:
                vals_df = pd.DataFrame(columns=range(1,5), index=range(1,5))
                species_order = [] # use like row/col labels in a symmetrical matrix
                for line_down in range(1,5):
                    val_line = all_log_lines[i+line_down].strip()
                    gID_10char = val_line[:11] # first 10 chars of geneIDs from leaf names

                    for gID_patt, sp in species_association_dict.items():
                        if gID_patt in gID_10char:
                            species_order.append(sp)
                    
                        if line_down==1:
                            continue # symmetrical matrix has no entry at (0,0)
                        else:
                            dNdS_vals = val_line.split()
                            for col_count, dNdS_val in enumerate(dNdS_vals[1:]): # exclude geneID
                                row_ind = col_count+1
                                vals_df.loc[line_down,row_ind] = float(dNdS_val)
                vals_df.rename(columns={ i+1: species_order[i] for i in range(4) }, index={ i+1: species_order[i] for i in range(4) }, inplace=True)

                out_dict[ortholog] = vals_df
                curr_OG_dNdS=False

    return out_dict


def get_pair_dNdS_list(orthologs_dict, species_list, dNdS_max=2):
    """
    returns a list of dNdS values from each pair based on an input dict from read_branch_log_outfile()
    """

    print(f"only dNdS values less or equal to {dNdS_max} are included! ")
    species_list=sorted(species_list)
    pairs_list = []
    for species1 in species_list:
        for species2 in species_list:
            if species1==species2:
                continue
            pairs_list.append(tuple(sorted([species1,species2])))
    pairs_list = sorted([f"{sp1}_{sp2}" for sp1,sp2 in list(set(pairs_list))])
    out_lists = { pair : [] for pair in pairs_list}

    for ortholog, vals_df in orthologs_dict.items():
        for pair in pairs_list:
            f1,s1,f2,s2 = pair.split("_")
            species1 = f"{f1}_{s1}"
            species2 = f"{f2}_{s2}"
            if pd.isna(vals_df.loc[species1,species2]):
                species2 = f"{f1}_{s1}"
                species1 = f"{f2}_{s2}"
            dNdS_val = vals_df.loc[species1,species2]
            if dNdS_val <= dNdS_max:
                out_lists[pair].append(dNdS_val)

    return out_lists 
        
    


if __name__ == "__main__":

    username = "miltr339"
    codeml_outfiles = filepaths(username=username)
    
    geneID_species_association = {
        "AOBTE" : "A_obtectus",      
        "BRAKERILHT" : "B_siliquastri", 
        "CALCHI" : "C_chinensis", 
        "g" : "C_maculatus",
    }
    species_list=sorted(["C_maculatus", "C_chinensis", "B_siliquastri", "A_obtectus"])
    dNdS_dict={}

    for chromosome in ["A", "X"]:
        print(f"\n////////////////////////// {chromosome} //////////////////////////")

        orthologs_vals = read_branch_log_outfile(filename=codeml_outfiles[f"branch_model_{chromosome}"], species_association_dict=geneID_species_association)

        # for ortholog, vals_df in orthologs_vals.items():
        #     print(f"\n------------- {ortholog} -------------")
        #     print(vals_df)
        #     break
        
        pair_vals = get_pair_dNdS_list(orthologs_dict=orthologs_vals, species_list=species_list)

        for pair, vals_list in pair_vals.items():
            pair_med = np.median(vals_list)
            print(f"{pair} : \t{pair_med:.4f}")

        dNdS_dict[chromosome]=pair_vals

    ### test with 100, takes a bit of time otherwise
    ### actual analysis with 10000
    num_permutations = 10000
    print(f"\n\n run permutation test with {num_permutations} permutations:\n")
    for pair in dNdS_dict["A"].keys():

        dNdS_A = dNdS_dict["A"][pair]
        dNdS_X = dNdS_dict["X"][pair]
        median_diffs = np.nanmedian(dNdS_A) - np.nanmedian(dNdS_X)
        bootstraps = permutate_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X, num_permut=num_permutations)
        mean_cor,std_cor,lower_CI,upper_CI = calculate_list_CI(bootstraps)
        mean_boot = np.mean(bootstraps)
        if median_diffs<lower_CI or median_diffs>upper_CI:
            print(f" *  {pair} --> \t median(dNdS_A)-median(dNdS_X) = {median_diffs:.5f}, mean bootstrap diff = {mean_boot:.5f} with CI [{lower_CI:.5f},{upper_CI:.5f}] --> SIGNIFICANT")
        else:
            print(f" *  {pair} --> \t median(dNdS_A)-median(dNdS_X) = {median_diffs:.5f}, mean bootstrap diff = {mean_boot:.5f} with CI [{lower_CI:.5f},{upper_CI:.5f}] --> (nonsignificant)")