import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from bootstrap_dNdS import permutate_dNdS,calculate_list_CI
import numpy as np
import pandas as pd
from tqdm import tqdm
import analyze_site_classes

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


def make_4way_bin_tuples(four_way_orthologs_IDs:str, four_way_pos_sel:str, species_list:list):
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
    pos_sel_genes_sp_split = {species : {} for species in species_list}
    pos_sel_genes = []
    count = 0
    with open(four_way_orthologs_IDs, "r") as orthologs_IDs_file:
        for line in orthologs_IDs_file.readlines():
            OG_id,chromosome,Aobt,Bsil,Cchi,Cmac = line.strip().split(",")
            try:
                pos_sel = pos_sel_IDs_dict[OG_id]
            except: 
                count +=1
                continue
            pos_sel_genes.append(([Aobt,Bsil,Cchi,Cmac],pos_sel))
            pos_sel_genes_sp_split["A_obtectus"][Aobt] = pos_sel
            pos_sel_genes_sp_split["B_siliquastri"][Bsil] = pos_sel
            pos_sel_genes_sp_split["C_maculatus"][Cchi] = pos_sel
            pos_sel_genes_sp_split["C_chinensis"][Cmac] = pos_sel
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
            if gene_id1 == "rna-AOBTE_LOCUS493" or gene_id2 == "rna-AOBTE_LOCUS493":
                print(line, pos_sel)
    print(f"{count} pairwise orthogroups failed the paml run")
    return pos_sel_genes



def plot_pos_sel_overlap_bruchini(four_way_orthologs_IDs:str, four_way_pos_sel:str, pairwise_orthologs_IDs:str, pairwise_pos_sel:str, species_list:list, pairs_list:list, filename = ""):
    """
    plot the number of pairwise positive selection for every 4-way ortholog
    """
    four_way_pos_sel_genes = make_4way_bin_tuples(four_way_orthologs_IDs=four_way_orthologs_IDs, four_way_pos_sel=four_way_pos_sel, species_list=species_list)
    pairwise_pos_sel_genes = make_pairwise_bin_dict(pairwise_orthologs_IDs=pairwise_orthologs_IDs, pairwise_pos_sel=pairwise_pos_sel, species_list=species_list)

    print(f"\nFOUR-WAY ANALYSIS")
    count = 0
    for geneID_list, pos_sel in four_way_pos_sel_genes:
        if pos_sel:
            count += 1
    print(f"\t- {count} out of {len(four_way_pos_sel_genes)} positively selected")

    print(f"\nPAIRWISE ANALYSIS")
    for pair, pos_sel_dict in pairwise_pos_sel_genes.items():
        pos_sel_list = list(pos_sel_dict.values())
        print(f"\t- {pair}: {sum(pos_sel_list)} out of {len(pos_sel_list)} positively selected")


    ### TODO something is very wrong from here down and i can't work it out

    print("\n")
    print("""
check for every four-way ortholog how many of the pairs are under positive selection.
Check separately for four-way orthologs under positive selection and four-way orthologs under negative selection:
    """)
    pairs = list(pairwise_pos_sel_genes.keys())
    #print(pairwise_pos_sel_genes[pairs[0]])
    four_way_pos_stats = {i:0 for i in range(len(pairs_list)+1)}
    four_way_neg_stats = {i:0 for i in range(len(pairs_list)+1)}
    
    for geneID_list, four_way_pos_sel in four_way_pos_sel_genes:
        Aobt,Bsil,Cchi,Cmac=geneID_list
        geneID_dict = {"C_maculatus":Cmac, "C_chinensis":Cchi, "B_siliquastri":Bsil, "A_obtectus":Aobt}
        count_pos_miss = 0
        count_pos_pair = 0
        
        pairs_done = []
        for species,geneID in geneID_dict.items():
            for pair in pairs_list:
                if species not in pair or pair in pairs_done:
                    continue
                
                try:
                    pairwise_pos_sel_genes[pair][geneID]
                except:
                    count_pos_miss+=1
                    continue
                if pairwise_pos_sel_genes[pair][geneID]==True:
                    count_pos_pair+=1
                pairs_done.append(pair)
        
        if four_way_pos_sel==False:
            four_way_neg_stats[count_pos_pair]+=1
        else:
            four_way_pos_stats[count_pos_pair]+=1
        

    four_way_pos_stats_str = "\n".join([f"\t\tpositively selected in {num_pairs} of 3 pairs: {pos_sel} orthologs" for num_pairs,pos_sel in four_way_pos_stats.items() if num_pairs != "miss"])
    four_way_neg_stats_str = "\n".join([f"\t\tpositively selected in {num_pairs} of 3 pairs: {pos_sel} orthologs" for num_pairs,pos_sel in four_way_neg_stats.items() if num_pairs != "miss"])
    print(f"\t - 4-way pos. sel genes: \n{four_way_pos_stats_str}")
    print(f"\t - 4-way neg. sel genes: \n{four_way_neg_stats_str}")

    fig, ax = plt.subplots(1,2, figsize=(30, 15))
    fs = 35
    plt.rcParams.update({'font.size': fs})

    ylab = f"num 4-way orthologs that have\nthis number of pos. sel. pairs"
    xlab = f"num of pairs in each 4-way ortholog"
    ax[0].bar(x=four_way_pos_stats.keys(), height=four_way_pos_stats.values(), color="#BD351E")
    ax[0].set_title(f"4-way positively sel. genes", fontsize=fs)
    ax[0].set_ylabel(ylab, fontsize = fs)
    ax[0].set_xlabel(xlab, fontsize = fs)
    ax[0].set_xticks(list(four_way_pos_stats.keys()))
    ax[0].set_xticklabels(list(four_way_pos_stats.keys()))
    ax[0].tick_params(axis='x', labelsize=fs)
    ax[0].tick_params(axis='y', labelsize=fs)

    ax[1].bar(x=four_way_neg_stats.keys(), height=four_way_neg_stats.values(), color="#BD351E")
    ax[1].set_title(f"4-way not pos. sel. genes", fontsize=fs)
    ax[1].set_xlabel(xlab, fontsize = fs)
    ax[1].set_xticks(list(four_way_neg_stats.keys()))
    ax[1].set_xticklabels(list(four_way_neg_stats.keys()))
    ax[1].tick_params(axis='x', labelsize=fs)
    ax[1].tick_params(axis='y', labelsize=fs)

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)


    


def make_filtered_pairs_list(pairs_lookup_table:str, four_way_lookup_table:str):
    """
    Filter the pairs lookup table to a new file that only includes pairs where both members are also part of a 4-way ortholog
    """
    #read 4-way lookup table
    four_way_df = pd.read_csv(four_way_lookup_table)
    outfile_filt_pairs_name = pairs_lookup_table.replace(".txt", "_only_4way_orthologs.txt")
    four_way_dict=four_way_df.to_dict(orient="list")
    
    with open(pairs_lookup_table, "r") as pairs_infile, open(outfile_filt_pairs_name, "w") as pairs_outfile_filt:
        count_filt = 0
        count_all = 0
        for line in tqdm(pairs_infile.readlines()):
            if "not_found" in line:
                continue
            count_all += 1
            OGs,IDs = line.strip().split(":")
            lines_list = OGs.split("_")
            sp1 = f"{lines_list[0]}_{lines_list[1]}"
            sp2 = f"{lines_list[2]}_{lines_list[3]}"
            gID1,gID2 = IDs.split(",")
        
            if gID1 in four_way_dict[sp1] and gID2 in four_way_dict[sp2]:
                count_filt +=1
                pairs_outfile_filt.write(line)

    print(f"out of {count_all} pairwise orthologs, {count_filt} are in 4-way Bruchini orthologs")
    print(f"outfile written to :\n{outfile_filt_pairs_name}")



def filt_site_classes(filt_lookup_table_file, site_classes_table):
    site_classes_table_filt = site_classes_table.replace(".txt", "_only_4way_orthologs.txt")
    
    # get orthologIDs list from lookup table
    with open(filt_lookup_table_file, "r") as filt_lookup_table, open(site_classes_table, "r") as site_classes, open(site_classes_table_filt, "w") as site_classes_filt:
        four_way_ortholog_IDs =[line.strip().split(":")[0] for line in filt_lookup_table.readlines()] 
        count_filt = 0
        count_all = 0
        for line in tqdm(site_classes.readlines()):
            orthologID = line.strip().split("_dNdS")[1].split("/")[-1]
            count_all +=1
            if orthologID in four_way_ortholog_IDs:
                count_filt += 1
                site_classes_filt.write(line)

    
    print(f"out of {count_all} pairwise orthologs, {count_filt} are in 4-way Bruchini orthologs")
    print(f"outfile written to :\n{site_classes_table_filt}")


def make_orthologs_association_file(four_way_lookup_table, pairs_lookup_table, outfile):
    # associate ortholog groups from lookup table
        #read 4-way lookup table
    four_way_df = pd.read_csv(four_way_lookup_table)
    four_way_df["four_way_OG_ID"] = four_way_df["orthogroup_ID"]
    four_way_df = four_way_df.set_index("orthogroup_ID")
    four_way_df["pairs_OG_IDs_lists"] = [ [] for i in range(len(four_way_df))]
    print(four_way_df)
    
    pairs_df = pd.read_csv(pairs_lookup_table, sep=":|,", names=["pair_ortholog_ID", "species1_gID", "species2_gID"])
    pairs_df = pairs_df.set_index("pair_ortholog_ID")
    print(pairs_df)

    for pair_ortholog, (species1_gID, species2_gID) in tqdm(pairs_df.iterrows()):
        sp_list = pair_ortholog.split("_")
        sp1 = f"{sp_list[0]}_{sp_list[1]}"
        sp2 = f"{sp_list[2]}_{sp_list[3]}"

        # get_orthogroup_index
        sp1_4way = four_way_df.index[four_way_df[sp1] == species1_gID][0]
        sp2_4way = four_way_df.index[four_way_df[sp2] == species2_gID][0]
        if sp1_4way != sp2_4way:
            print(f"!!!! {pair_ortholog} --> {sp1}:{species1_gID}:4way_ID={sp1_4way} + {sp2}:{species2_gID}:4way_ID={sp2_4way} ")
            raise RuntimeError(f"somehow the above pairwise geneIDs aren't part of the same four-way orthogroup")
        four_way_df.loc[sp1_4way, "pairs_OG_IDs_lists"].append(pair_ortholog)
    
    print(four_way_df)
    four_way_df.to_csv(outfile, index=False)
    print(f"outfile saved to: {outfile}")
                


def plot_omega_association(ortholog_IDs_association, site_classes_pairs, site_classes_4way):
    orthologIDs_df = pd.read_csv(ortholog_IDs_association)
    print(orthologIDs_df)

    ## TODO read site classes with analyze_site_classes.read_site_classes(summary_path)



if __name__ == "__main__":

    username = "miltr339"
    paths = in_paths(username=username)

    if False:
        # plot the correlation of codeml dN and dS values for pairwise comparisons of runmode=1 and runmode=-2
        plot_correlation(def_file=paths["default"], pair_file=paths["pairwise"], filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/revision_tests/codeml_dNdS_correlation.png")
    if False:
        # do permutation test on the 4-way 1-to-1 orthologs in Bruchini
        bootstrap_pos_sel(X_data={"pos_sel" : 6 , "non_pos_sel" : 268}, A_data={"pos_sel" : 60 , "non_pos_sel" : 1881}, num_permutations=10000)
    
    four_way_ortholog_IDs=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_orthologs.csv"
    four_way_pos_sel=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_bruchini_all_orthologs_pos_sel_no_err.txt"
    pairwise_ortholog_IDs=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/pairwise_ortholog_IDs_association.txt"
    pairwise_ortholog_IDs_4way_filt = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/pairwise_ortholog_IDs_association_only_4way_orthologs.txt"
    pairwise_pos_sel=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/pairwise_site_classes_summary.txt"
    site_classes_files = {
            "A_LRT" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_A_BH_corrected_only_4way_orthologs.txt",
            "X_LRT" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_X_BH_corrected_only_4way_orthologs.txt",
        }
    species_list=sorted(["C_maculatus", "C_chinensis", "B_siliquastri", "A_obtectus"])
    pairs_list = []
    for species1 in species_list:
        for species2 in species_list:
            if species1==species2:
                continue
            pairs_list.append(tuple(sorted([species1,species2])))
    pairs_list = sorted([f"{sp1}_{sp2}" for sp1,sp2 in list(set(pairs_list))])
    
    if False:
        # Compare positively selected orthologs in the pairwise tests to if they are also positively selected in the 4-way comparison
        plot_pos_sel_overlap_bruchini(four_way_orthologs_IDs=four_way_ortholog_IDs, four_way_pos_sel=four_way_pos_sel, 
            pairwise_orthologs_IDs=pairwise_ortholog_IDs, pairwise_pos_sel=pairwise_pos_sel, 
            species_list=species_list, pairs_list=pairs_list, 
            filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/revision_tests/codeml_pos_sel_ortholog_overlap.png")

    if False:
        # filter pair IDs of pairwise orthologs to include only ones that are present in 4-way orthologs
        
        # make_filtered_pairs_list(pairs_lookup_table=pairwise_ortholog_IDs, four_way_lookup_table=four_way_ortholog_IDs)
        
        # fitler site classes to only include orthologs that are in the above list:
        site_classes_files = {
            "A_LRT_BH_corr" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_A_BH_corrected.txt",
            "X_LRT_BH_corr" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_X_BH_corrected.txt",
        }
        filt_site_classes(filt_lookup_table_file = pairwise_ortholog_IDs_4way_filt, site_classes_table = site_classes_files["X_LRT_BH_corr"])
        filt_site_classes(filt_lookup_table_file = pairwise_ortholog_IDs_4way_filt, site_classes_table = site_classes_files["A_LRT_BH_corr"])

    if True:
        # correlate w_0 for the positively selected genes
        pairs_fourWay_association = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/ortholog_IDs_geneIDs_4_way_and_pairs_Bruchini.txt"
        # make_orthologs_association_file(four_way_lookup_table=four_way_ortholog_IDs, pairs_lookup_table=pairwise_ortholog_IDs_4way_filt, outfile = pairs_fourWay_association)
        plot_omega_association(ortholog_IDs_association = pairs_fourWay_association, site_classes_pairs = site_classes_files["X_LRT"], site_classes_4way=four_way_pos_sel)
