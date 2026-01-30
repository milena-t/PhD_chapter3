"""
plot the results extracted by src/blast_BRH/extract_dNdS_Results.get_site_classes()
The results file contains the path to the ortholog, with pair and number ID, M1a or M2a model (depending on if M2a is sig. better than M1a in the LRT), 
and then p and w from the site classes table. Either two or three site classes, depending on if M2a (3 classes) was significantly better than M1a (2 classes)
"""

from plot_dNdS import get_summary_paths,read_dNdS_summary_file,get_species_list
import os
import numpy as np
import bootstrap_dNdS 
from matplotlib.ticker import FuncFormatter

site_classes_files = {
    "A" : "/Users/miltr339/work/chapter3/dNdS/site_classes_summary_A-linked.txt",
    "X" : "/Users/miltr339/work/chapter3/dNdS/site_classes_summary_X-linked.txt",
}


class SiteClassesTable:
    """
    General class for the site classes tables that are calculated for each pairwise comparison
    contains several utilities for basic operations
    """
    def __init__(self, p_neg:float, p_neutr:float, p_pos:float, w_neg:float, w_neutr:float, w_pos:float) -> None:
        """
        if site classes table from M1a, pass 0.0 to p_pos and w_pos! 
        """
        self.p_neg = p_neg
        self.p_neutr = p_neutr
        self.p_pos = p_pos
        self.w_neg = w_neg
        self.w_neutr = w_neutr
        self.w_pos = w_pos

    @property
    def sig_pos_selection(self):
        if self.p_pos==0.0 and self.w_pos==0.0:
            return False
        else:
            return True

    @property
    def mean_w(self):
        """
        weighted average of w using p as weights
        """
        return np.average([self.w_neg,self.w_neutr,self.w_pos], weights=[self.p_neg,self.p_neutr,self.p_pos])


    def __str__(self) -> str:
        if self.sig_pos_selection:
            return f"significant positive selection detected!\np:\t{self.p_neg}\t{self.p_neutr}\t{self.p_pos}\nw:\t{self.w_neg}\t{self.w_neutr}\t{self.w_pos}\nweighted avg. w = {self.mean_w}"
        else:
            return f"no positive selection detected.\np:\t{self.p_neg}\t{self.p_neutr}\t{self.p_pos}\nw:\t{self.w_neg}\t{self.w_neutr}\t{self.w_pos}\nweighted avg. w = {self.mean_w}"





class Ortholog:
    """
    class that contains all information about an ortholog and the site classes table
    """
    def __init__(self, pair:list, ortholog_number:int, site_classes:SiteClassesTable) -> None:
        assert len(pair) ==2
        self.pair = pair
        self.ortholog_number = ortholog_number
        self.site_classes = site_classes

    def __str__(self) -> str:
        site_classes_string = f"\tp:\t{self.site_classes.p_neg}\t{self.site_classes.p_neutr}\t{self.site_classes.p_pos}\n\tw:\t{self.site_classes.w_neg}\t{self.site_classes.w_neutr}\t{self.site_classes.w_pos}"
        return f"1-to-1 ortholog nr. {self.ortholog_number} between {self.pair[0]} and {self.pair[1]}, site class table:\n{site_classes_string}"



def species_names_from_pair(pair_string):
    """
    'G_species1_G_species2_...' -> 'G_species1','G_species2'
    """
    try:
        pair_list = pair_string.split("_")
        g1, s1, g2, s2 = pair_list[:4]
    except:
        raise RuntimeError(f"invalid species pair! {pair_string}")
    
    return f"{g1}_{s1}", f"{g2}_{s2}"


def get_pairs_from_summary(summary_path):
    """
    get a list of species pair sets from the summary path
    """
    if not os.path.isfile(summary_path):
        raise RuntimeError(f"FILE: {summary_path} does not exist!")
    

    with open(summary_path, "r") as summary_file:
        lines = summary_file.readlines()
        pairs_list = ['']*len(lines)
        for i, line in enumerate(lines):
            filepath = line.strip().split(" : ")[0]
            filepath = filepath.split("/")
            pair_name = filepath[-3]
            species1,species2=species_names_from_pair(pair_name)
            if "-" in species1 or "-" in species2:
                raise RuntimeError(f"wrong parsing!\n {line} \n-> {filepath} \n-> {pair_name}")
            pairs_list[i] = f"{species1}_{species2}"
        pairs_list_unique = list(set(pairs_list))
    
    return pairs_list_unique


def read_site_classes(summary_path):
    """
    read the site classes summary file into a data structure
    """
    if not os.path.isfile(summary_path):
        raise RuntimeError(f"FILE: {summary_path} does not exist!")

    out_dict = {pair : [] for pair in get_pairs_from_summary(summary_path)}
    no_dNdS = {pair : [] for pair in get_pairs_from_summary(summary_path)}
    print(f" *  {len(out_dict)} unique pairs")
    
    with open(summary_path, "r") as summary_file:
        for line in summary_file.readlines():
            try:
                filepath_str, site_classes_string = line.strip().split(" : ")
            except:
                raise RuntimeError(f"line cannot be parsed! \n{line}")
                # continue
            filepath = filepath_str.split("/")
            
            pair_name,site_class_name = filepath[-2:]
            species1,species2=species_names_from_pair(pair_name)
            pair = f"{species1}_{species2}"
            model = site_class_name.split("_")[1]
            try:
                pair_number = int(pair_name.split("_")[-2])
            except:
                if pair_name[-4:] == ".out":
                    continue
                else:
                    raise RuntimeError(f"cannot parse pair number in {pair_name}, from filepath: {filepath}")

            try:
                p_string,w_string = site_classes_string.split(";")
            except:
                no_dNdS[pair].append(filepath_str)
                continue

            p_list = [float(val) for val in p_string.split()[1:]]
            w_list = [float(val) for val in w_string.split()[1:]]
            assert len(p_list) == len(w_list)
            if model == "M1a" and len(p_list) == 2:
                p_list.append(0.0)
                w_list.append(0.0)
            
            stc = SiteClassesTable(p_neg=p_list[0], p_neutr=p_list[1], p_pos=p_list[2], w_neg=w_list[0], w_neutr=w_list[1], w_pos=w_list[2])
            og = Ortholog(pair=[species1,species2], ortholog_number=pair_number, site_classes=stc)

            out_dict[pair].append(og)

    return out_dict, no_dNdS


def count_pos_sel_genes(summary_dict, full_list = False):
    """
    takes a summary dict created by read_site_classes and returns the total number of genes investigated 
    and the number of genes with positively selected sites according to the LRT
    if full_list, it doesn't return counts but a binary list of 1 and 0 depending on if a gene is positively selected or not
    """
    out_dict = {pair : [np.nan, np.nan] for pair in summary_dict.keys()}

    for pair, orthologs_list in summary_dict.items():
        binary_list = [np.nan]*len(orthologs_list)
        count_all = 0
        count_pos = 0
        for i, ortholog in enumerate(orthologs_list):
            count_all+=1
            if ortholog.site_classes.sig_pos_selection:
                count_pos+=1
                binary_list[i]=1
            else:
                binary_list[i]=0

        if full_list:
            out_dict[pair] = binary_list
        else:    
            out_dict[pair][0] = count_all
            out_dict[pair][1] = count_pos
    
    return out_dict


def avg_prop_pos_sel_sites(summary_dict, full_list = False):
    """
    takes a summary dict created by read_site_classes and returns the average proportion of positively selected sites 
    """
    out_dict = {pair : [np.nan] for pair in summary_dict.keys()}

    for pair, orthologs_list in summary_dict.items():
        proportion_list = [np.nan]*len(orthologs_list)

        for i, ortholog in enumerate(orthologs_list):
            proportion_list[i] = ortholog.site_classes.p_pos

        if full_list:
            out_dict[pair] = proportion_list
        else:
            out_dict[pair] = np.nanmean(proportion_list)
    
    return out_dict



def binary_barplot_pair(data_A, data_X, row, col, n_A, n_X, axes, colors_dict, fs, ylab ="percent", xticks = ["A", "X"]):
    """
    plot a bar chart to show proportion of genes that have a positive site-class model LRT
    """
    A_nonsig = len([val for val in data_A if val == 0.0])*100.0
    X_nonsig = len([val for val in data_X if val == 0.0])*100.0
    nonsig = [A_nonsig/len(data_A), X_nonsig/len(data_X)]
    A_sig = len([val for val in data_A if val == 1.0])*100.0
    X_sig = len([val for val in data_X if val == 1.0])*100.0
    sig = [A_sig/len(data_A), X_sig/len(data_X)]

    axes[row,col].bar([1,2], nonsig, width = 0.75, label='M1a', color=[colors_dict["A"], colors_dict["X"]], alpha=0.7)
    axes[row,col].bar([1,2], sig, width = 0.75, bottom=nonsig, label='M2a', color= [colors_dict["A"], colors_dict["X"]])
    print(f"\t - row: {row} col: {col}\t nonsignificant: A= {nonsig[0]:.2f}% X={nonsig[1]:.2f}%")
    print(f"\t - row: {row} col: {col}\t significant: A= {sig[0]:.2f}% X={sig[1]:.2f}%")

    axes[row, col].set_ylim([0,118])
    axes[row, col].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<1 else f'{int(x)}%'))

    axes[row, col].text(1-0.2, 105, f"n={n_A}", fontsize = fs, color = colors_dict["A"])
    axes[row, col].text(2-0.2, 105, f"n={n_X}", fontsize = fs, color = colors_dict["X"])
    
    axes[row, col].set_xticks([1,2])
    axes[row, col].set_xticklabels(xticks)

    axes[row, col].set_xlabel('')
    if col-row == 1:
        axes[row, col].set_ylabel(ylab, fontsize = fs*0.8)
    else:
        axes[row, col].set_ylabel('')
    axes[row,col].tick_params(axis='y', labelsize=fs)
    axes[row,col].tick_params(axis='x', labelsize=fs) 


if __name__ == "__main__":

    chr_types = ["X", "A"]
    username = f"miltr339" #"milena"

    if False:
        ## compute statistics to terminal
        for chr_type in chr_types:
            print(f"\n//////////////////// {chr_type} ////////////////////\n")

            summary_dict_X, no_dNdS_X = read_site_classes(site_classes_files[chr_type])
            pos_counts = count_pos_sel_genes(summary_dict_X)
            mean_pos_prop = avg_prop_pos_sel_sites(summary_dict_X)
            print(f" *  number of positively selected genes")
            for pair, counts_list in pos_counts.items():
                avg_prop = mean_pos_prop[pair]
                try:
                    ratio = counts_list[1]*100.0/counts_list[0]
                    print(f"\t{pair} : \t {counts_list[1]} of {counts_list[0]} ({ratio:.2f} %) genes have positively selected sites, avg. {avg_prop:.3}")
                except:
                    print(f"\t{pair} : \t {counts_list[1]} of {counts_list[0]} genes have positively selected sites")

            print(f"\n *  missing data:")
            for pair, absent_list in no_dNdS_X.items():
                print(f"\t{pair} : \t {len(absent_list)} genes could not be computed, {pos_counts[pair][0]} are computed")

    ## plot 


    username = "milena"# "miltr339"
    chromosome = "A"
    data_files = {"A" : ["A_dNdS", "A_LRT"],
                  "X" : ["X_dNdS", "X_LRT"]}
    summary_paths = get_summary_paths(username=username)

    summary_dict_A, no_dNdS_A = read_site_classes(summary_paths[data_files["A"][1]])
    summary_dict_X, no_dNdS_X = read_site_classes(summary_paths[data_files["X"][1]])


    type_plot = "bin" # bin for binary or prop for proportion of sites

    if type_plot == "bin":
        pos_list_A = count_pos_sel_genes(summary_dict_A, full_list=True)
        pos_list_X = count_pos_sel_genes(summary_dict_X, full_list=True)
    elif type_plot == "prop":
        pos_list_A = avg_prop_pos_sel_sites(summary_dict_A, full_list=True)
        pos_list_X = avg_prop_pos_sel_sites(summary_dict_X, full_list=True)

    pairs_list = list(summary_dict_X.keys())

### specify species subsets
    # pairs_list = [pair for pair in pairs_list if "T_castaneum" in pair]
    # plot_name = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/fastX_{type_plot}_pos_sites_permutation_Tcas_Tfre.png"

    # pairs_list = [pair for pair in pairs_list if "C_magnifica" in pair]
    # plot_name = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/fastX_{type_plot}_pos_sites_permutation_Csep_Cmag.png"

    pairs_list = [pair for pair in pairs_list if "C_maculatus" in pair or "C_chinensis" in pair or "A_obtectus" in pair or "B_siliquastri" in pair]
    plot_name = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/fastX_{type_plot}_pos_sites_permutation_Bruchini.png"
###

    bootstraps = {pair : [] for pair in pairs_list}
    mean_num_pos_sel = {pair : np.NaN for pair in pairs_list}
    
    ####
    #  test with 100, takes a bit of time otherwise
    num_permutations = 10000
    ####

    for pair in pairs_list:
        pos_props_A = pos_list_A[pair]
        pos_props_X = pos_list_X[pair]
        mean_num_pos_sel[pair] = np.nanmean(pos_props_A) - np.nanmean(pos_props_X)
        bootstraps[pair] = bootstrap_dNdS.permutate_dNdS(dNdS_A=pos_props_A, dNdS_X=pos_props_X, num_permut=num_permutations)
        mean_boot = np.mean(bootstraps[pair])
        print(f" *  {pair} mean({type_plot}_A)-mean({type_plot}_X)  --> \t{mean_num_pos_sel[pair]:.4f}, mean bootstrap diff {mean_boot:.6f}")
        
    if type_plot == "prop":
        violin_ymax=0.2
        binary=False
    else:
        violin_ymax=0
        binary=True

    bootstrap_dNdS.plot_dNdS_permutations(boot_diff=bootstraps,measure_diff=mean_num_pos_sel, A_dict=pos_list_A, X_dict=pos_list_X, 
                                            filename=plot_name, hist_label = f"mean({type_plot}_A)-mean({type_plot}_X)", violin_label=f"{type_plot}. pos. sel.", 
                                            violin_ymax=violin_ymax, transparent=False, binary=binary)
    