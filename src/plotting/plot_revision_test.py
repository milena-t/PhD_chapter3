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


def make_4way_bin_dict(four_way_orthologs_IDs:str, four_way_pos_sel:str, species_list:list):
    """
    make dict for every geneID in the four-way orthologs that assigns positive selection, separated by species
    returns { species : {geneID:True, geneID:False, ...}, ...} with True and False showing positive selection
    """
    sel_bin_dir = {"codeml_M1a_site_class_table.out" : False, "codeml_M2a_site_class_table.out" : True}
    # read pos_sel_ortholog_IDs into dict like {ortholog_ID : True} for pos. sel.=True or false
    pos_sel_IDs_dict = {}
    with open(four_way_pos_sel, "r") as pos_sel_IDs_file:
        for line in pos_sel_IDs_file.readlines():
            line = line.strip().split(" : ")[0].split("linked_ortholog_")[1]
            OG_id,paml_outfile = line.split("_dNdS/")
            pos_sel_IDs_dict[OG_id] = sel_bin_dir[paml_outfile]

    # get dict with pos sel IDs
    pos_sel_genes = {species : {} for species in species_list}
    count = 0
    with open(four_way_orthologs_IDs, "r") as orthologs_IDs_file:
        for line in orthologs_IDs_file.readlines():
            OG_id,chromosome,Aobt,Bsil,Cchi,Cmac = line.strip().split(",")
            try:
                pos_sel = pos_sel_IDs_dict[OG_id]
            except: 
                count +=1
                continue
            pos_sel_genes["A_obtectus"][Aobt] = pos_sel
            pos_sel_genes["B_siliquastri"][Bsil] = pos_sel
            pos_sel_genes["C_maculatus"][Cchi] = pos_sel
            pos_sel_genes["C_chinensis"][Cmac] = pos_sel
    print(f"{count} 4-way orthogroups failed the paml run or i did not run them because i was impatient")
    return pos_sel_genes


def make_pairwise_bin_dict(pairwise_orthologs_IDs:str, pairwise_pos_sel:str, species_list:list):
    """
    make dict for every geneID in the pairwise orthologs that assigns positive selection
    returns { pair : {gene_id:True, gene_id:False, ...}, ...}
    """
    # read pos_sel_ortholog_IDs into dict like {pair : {ortholog_ID:True, ...}, ...} for pos.sel.=True or False
    sel_bin_dir = {"codeml_M1a_site_class_table.out" : False, "codeml_M2a_site_class_table.out" : True}
    pos_sel_IDs_dict = {}
    with open(pairwise_pos_sel, "r") as pos_sel_IDs_file:
        for line in pos_sel_IDs_file:
            line = line.strip()
            # get pair
            pair = line.split("_pairwise_dNdS")[0].split("/")[-1]
            # get ID
            OG_id,paml_outfile = line.split(" : ")[0].split("linked_ortholog_")[1].split("_dNdS/")
            if pair not in pos_sel_IDs_dict.keys():
                pos_sel_IDs_dict[pair] = {OG_id : sel_bin_dir[paml_outfile]}
            else:
                pos_sel_IDs_dict[pair][OG_id] = sel_bin_dir[paml_outfile]
    
    pos_sel_genes = {pair : {} for pair in pos_sel_IDs_dict.keys()} # make { pair : {gene_id:True, gene_id:False, ...}, ...}
    count = 0            
    with open(pairwise_orthologs_IDs, "r") as orthologs_IDs_file:
        for line in orthologs_IDs_file.readlines():
            line=line.strip()
            if "not_found" in line:
                continue
            # get species_pair
            sp = line.split("_")
            sp1 = f"{sp[0]}_{sp[1]}"
            sp2 = f"{sp[2]}_{sp[3]}"
            pair = f"{sp1}_{sp2}"
            OG_id = line.split(":")[0].split("linked_ortholog_")[-1]
            try:
                pos_sel = pos_sel_IDs_dict[pair][OG_id]
            except:
                count +=1
            try:
                gene_id1,gene_id2=line.split(":")[1].split(",")
            except:
                print(f"-------- {line} --------")
                raise RuntimeError
            pos_sel_genes[pair][gene_id1] = pos_sel
            pos_sel_genes[pair][gene_id2] = pos_sel
    print(f"{count} pairwise orthogroups failed the paml run")
    return pos_sel_IDs_dict



def plot_pos_sel_overlap_bruchini(four_way_orthologs_IDs:str, four_way_pos_sel:str, pairwise_orthologs_IDs:str, pairwise_pos_sel:str, species_list:list):
    """
    plot the number of pairwise positive selection for every 4-way ortholog
    """
    four_way_pos_sel_genes = make_4way_bin_dict(four_way_orthologs_IDs=four_way_orthologs_IDs, four_way_pos_sel=four_way_pos_sel, species_list=species_list)
    pairwise_pos_sel_genes = make_pairwise_bin_dict(pairwise_orthologs_IDs=pairwise_orthologs_IDs, pairwise_pos_sel=pairwise_pos_sel, species_list=species_list)

    print(f"FOUR-WAY ANALYSIS")
    for species, pos_sel_dict in four_way_pos_sel_genes.items():
        pos_sel_list = list(pos_sel_dict.values())
        print(f"\t- {species}: {sum(pos_sel_list)} out of {len(pos_sel_list)} positively selected")

    print(f"\nPAIRWISE ANALYSIS")
    for pair, pos_sel_dict in pairwise_pos_sel_genes.items():
        pos_sel_list = list(pos_sel_dict.values())
        print(f"\t- {pair}: {sum(pos_sel_list)} out of {len(pos_sel_list)} positively selected")


    ### TODO something is very wrong from here down and i can't work it out
    
    print("\n\n")
    pairs = list(pairwise_pos_sel_genes.keys())
    four_way_pos_stats = {species : {0:0, 1:0, 2:0, 3:0, "miss":0} for species in four_way_pos_sel_genes.keys()} # count the number of positively selected pairs when positive selection is detected in the four-way comparison
    four_way_neg_stats = {species : {0:0, 1:0, 2:0, 3:0, "miss":0} for species in four_way_pos_sel_genes.keys()}
    for species, four_way_pos_sel_dict in four_way_pos_sel_genes.items():
        species_pairs = [pair for pair in pairs if species in pair]
        print(f"{species}: {species_pairs}")
    
        count_pos_miss = 0
        count_neg_miss = 0
        for four_way_species_geneID,four_way_pos_sel in four_way_pos_sel_dict.items():
            count_pos_pair = 0
            count_neg_pair = 0
            if four_way_pos_sel==False:
                continue
                # for pair in species_pairs:
                #     try:
                #         pairwise_pos_sel_genes[pair][four_way_species_geneID]
                #     except:
                #         count_neg_miss+=1
                #         continue
                #     if pairwise_pos_sel_genes[pair][four_way_species_geneID]:
                #         count_pos_pair+=1

            else:
                for pair in species_pairs:
                    try:
                        pairwise_pos_sel_genes[pair][four_way_species_geneID]
                    except:
                        count_pos_miss+=1
                        # print(four_way_species_geneID)
                        continue
                    if pairwise_pos_sel_genes[pair][four_way_species_geneID]==True:
                        count_pos_pair+=1
                    elif pairwise_pos_sel_genes[pair][four_way_species_geneID]==False:
                        count_neg_pair+=1
                    else:
                        print(four_way_species_geneID)
                        raise RuntimeError
                four_way_neg_stats[species][count_neg_pair]+=1
                four_way_pos_stats[species][count_pos_pair]+=1

        four_way_pos_stats[species]["miss"] = count_pos_miss
        four_way_neg_stats[species]["miss"] = count_neg_miss
        
        print(f"\t - 4-way pos. sel genes stats: {four_way_pos_stats[species]}")
        print(f"\t - 4-way neg. sel genes stats: {four_way_neg_stats[species]}")
                



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
        four_way_ortholog_IDs=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_orthologs.csv"
        four_way_pos_sel=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_all_orthologs_pos_sel_no_err.txt"
        pairwise_ortholog_IDs=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/pairwise_ortholog_IDs_association.txt"
        pairwise_pos_sel=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/pairwise_site_classes_summary.txt"
        
        plot_pos_sel_overlap_bruchini(four_way_orthologs_IDs=four_way_ortholog_IDs, four_way_pos_sel=four_way_pos_sel, 
        pairwise_orthologs_IDs=pairwise_ortholog_IDs, pairwise_pos_sel=pairwise_pos_sel, 
        species_list=["C_maculatus", "C_chinensis", "B_siliquastri", "A_obtectus"])