"""
Make a table for the phylogenetic rank analyses that has all the cmac genes in a row with the following columns
* geneID
* most distant conservation rank
* significant sex bias (1: female biased, -1: male biased, 0: unbiased) for
    * abdomen
    * head+thorax tissue 
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel
from plot_basic_DE_stuff import get_summary_paths

def get_full_table_path(username="miltr339"):
    out_dict = {
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_A_chr.tsv",
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_X_chr.tsv",
    }
    return out_dict


def make_rank_summary_table(full_table_path:str, outfile_path:str, min_LFC=1,p_threshold=0.05):
    """
    make summary table for conservation rank analyses
    """

    full_summary_df = pd.read_csv(full_table_path, sep="\t")
    full_summary_df.drop_duplicates(subset="focal_transcript", inplace=True) # I only need one instance since I don't include the dNdS values

    def classify_expression(log2FC, pvalue_FDR, min_log2FC=min_LFC, p_threshold=p_threshold):
        if pvalue_FDR > p_threshold:
            return 0
        elif log2FC > min_log2FC: # positive -> female bias
            return 1
        elif log2FC < min_log2FC*-1 : # negative -> male bias
            return -1
        else:
            return 0 ## some DE values are NaN (because expression is so low that they are filtered in edgeR and not included in the output)
    
    header_list = [
        "Cmac_transcript_ID",
        "conservation_rank",
        "abdomen_sex_bias_category",
        "head_thorax_sex_bias_category"
    ]
    header_string = "\t".join(header_list)

    with open(outfile_path, "w") as outfile:
        outfile.write(f"{header_string}\n")
        counts_dict_a = {1 : 0, -1 : 0, 0 : 0}
        counts_dict_ht = {1 : 0, -1 : 0, 0 : 0}
        for i,row in full_summary_df.iterrows():
            
            # classify expression
            abd_sex_bias = classify_expression(log2FC = row["LFC_abdomen"], pvalue_FDR=row["FDR_pval_abdomen"])
            ht_sex_bias = classify_expression(log2FC = row["LFC_head+thorax"], pvalue_FDR=row["FDR_pval_head+thorax"])
            counts_dict_a[abd_sex_bias]+=1
            counts_dict_ht[ht_sex_bias]+=1
            # other vals
            trans_ID = row["focal_transcript"]
            cons_rank = row["level_most_dist_ortholog"]

            out_row = f"{trans_ID}\t{cons_rank}\t{abd_sex_bias}\t{ht_sex_bias}\n"
            outfile.write(out_row)

    print(f"done writing to: {outfile_path}!")
    print(f"counts Abdomen: {counts_dict_a}")
    print(f"counts head+thorax: {counts_dict_ht}")

        

def categorize_SB_int(SB_int:int):
    if SB_int==0:
        return "unbiased"
    elif SB_int==-1:
        return "male_biased"
    elif SB_int==1:
        return "female_biased"
    else:
        raise RuntimeError(f"{SB_int} could not be categorized")

        
def make_category_proportion_plot(rank_summary_path_A:str, rank_summary_path_X:str, outfile:str):
    """
    plot a stacked bar to show the proportion of male/female/unbiased for X and A, and for abdominal and head+thorax
    """
    
    X_df = pd.read_csv(rank_summary_path_X, sep="\t")
    A_df = pd.read_csv(rank_summary_path_A, sep="\t")
    
    def modify_counts_dict(df):
        unique_ranks = df["conservation_rank"].unique()
        rank_sums_dict = { rank : {"abdomen" : {"female_biased" : 0,"male_biased" : 0,"unbiased" : 0},"head_thorax" : {"female_biased" : 0,"male_biased" : 0,"unbiased" : 0}} for rank in unique_ranks}

        for i, row in df.iterrows():
            row_rank = row["conservation_rank"]
            SB_int_a = row["abdomen_sex_bias_category"]
            SB_int_h = row["head_thorax_sex_bias_category"]
            cat_a = categorize_SB_int(SB_int_a)
            cat_h = categorize_SB_int(SB_int_h)
            # print(f"{cmac_ID}\t : abdomen sex bias category: {SB_int_a}={cat_a}, head+thorax: {SB_int_h}={cat_h}")
            # print(f"\t * {row_rank} --> {rank_sums_dict[row_rank]}")
            rank_sums_dict[row_rank]["abdomen"][cat_a] += 1
            rank_sums_dict[row_rank]["head_thorax"][cat_h] += 1
        return rank_sums_dict


    def make_percentage_dict(counts_dict, df):
        unique_ranks = sorted(list(df["conservation_rank"].unique()))
        perc_dict = { rank : {"abdomen" : {"female_biased" : 0,"male_biased" : 0,"unbiased" : 0},"head_thorax" : {"female_biased" : 0,"male_biased" : 0,"unbiased" : 0}} for rank in range(1,max(unique_ranks)+1)}
        gene_numbers_dict = { rank : 0 for rank in range(1,max(unique_ranks)+1)}
        
        for row_rank in perc_dict.keys():
            num_genes = df[df["conservation_rank"]==row_rank].shape[0]
            print(f"{row_rank} : {num_genes} genes")
            gene_numbers_dict[row_rank] = num_genes
            perc_dict[row_rank]["abdomen"] = { sex_bias : 100*float(val)/float(num_genes) for sex_bias,val in counts_dict[row_rank]["abdomen"].items()}
            perc_dict[row_rank]["head_thorax"] = { sex_bias : 100*float(val)/float(num_genes) for sex_bias,val in counts_dict[row_rank]["head_thorax"].items()}
        
        print(perc_dict)
        return perc_dict, gene_numbers_dict


    X_sums_dict = modify_counts_dict(X_df)
    print(f"\n ---> X rank percentages dict")
    X_perc_dict,X_num_genes_rank = make_percentage_dict(X_sums_dict, X_df)
    A_sums_dict = modify_counts_dict(A_df)
    print(f"\n ---> A rank percentages dict")
    A_perc_dict,A_num_genes_rank = make_percentage_dict(A_sums_dict, A_df)
    assert len(X_perc_dict) == len(A_perc_dict)

    ### plot stacked bar chart
    fs = 25
    width=0.2
    x_gap = 0.025
    # x_subtr=width*1.1/2.0
    x_coords = list(range(1,len(X_sums_dict.keys())+1))

    x_coords_bars = [
        value
        for x in x_coords
        for value in (
            x - width*1.5 - x_gap,
            x - width*0.5 - x_gap,
            x + width*0.5 + x_gap,
            x + width*1.5 + x_gap,
        )
    ]
    fig, ax = plt.subplots(1, 1, figsize=(22, 11)) 

    colors_dict ={
        "abdomen" : {
            "female_biased" : "#DC4141", # scarlet rush
            "male_biased" : "#7B8CE0", # wisteria blue
            "unbiased" : "#44455F", # vintage grape
        },
        "head_thorax" : { # colors in a lighter shade than abdomen
            "female_biased" : "#E86D6D", # coral rush
            "male_biased" : "#9DAAE9", # wisteria blue
            "unbiased" : "#5E5F7A", # dusty grape
        }
    }
    labels_dict ={
        "female_biased" : "female biased",
        "male_biased" : "male biased", 
        "unbiased" : "unbiased",
        "X" : "X-linked",
        "A" : "autosomal",
        "abdomen" : "reproductive",
        "head_thorax" : "somatic"
    }

    bottom = [0,0,0,0]*5
    subcat_tissue = ["abdomen","head_thorax"]
    chr_labels = [f"{labels_dict[tissue]}:{chr}" for chr in ["A","X"] for tissue in subcat_tissue]
    bar_labels = [
        (
            f"{sublabel} ({A_num_genes_rank[conserved_rank]})"
            if "A" in sublabel
            else f"{sublabel} ({X_num_genes_rank[conserved_rank]})"
        )
        for conserved_rank in x_coords
        for sublabel in chr_labels
    ]


    # stack bars from bottom to top
    for sex_bias in reversed(colors_dict["abdomen"].keys()):
        print(f"\t * {sex_bias}")
        # alternate the colors for the abdomen/ht bars
        colors_list = [colors_dict["abdomen"][sex_bias] if i%2==0 else colors_dict["head_thorax"][sex_bias] for i, _ in enumerate(range(len(bottom))) ]
        # make bar height for the four bars in all the conservation rank categories
        y_coords = []
        for conserved_rank in x_coords:
            y_coords.extend([A_perc_dict[conserved_rank]["abdomen"][sex_bias]]) 
            y_coords.extend([A_perc_dict[conserved_rank]["head_thorax"][sex_bias]]) 
            y_coords.extend([X_perc_dict[conserved_rank]["abdomen"][sex_bias]])
            y_coords.extend([X_perc_dict[conserved_rank]["head_thorax"][sex_bias]])

        assert len(x_coords_bars)==len(y_coords)==len(colors_list)
        ax.bar(x_coords_bars, y_coords, width = width, bottom=bottom, color= colors_list)
        
        bottom = [bottom[i]+y_coords[i] for i in range(len(y_coords))]

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<1 else f'{int(x)}%'))

    ax.tick_params(axis='y', labelsize=fs)
    ax.tick_params(axis='x', labelsize=fs*0.75, rotation=90)
    ax.set_xticks(ticks=x_coords_bars, labels=bar_labels)

    ax2 = ax.secondary_xaxis('bottom')
    ax2.set_xticks(x_coords)
    ax2.set_xticklabels(x_coords, fontsize=fs*1.2)
    ax2.spines['bottom'].set_position(('outward', 180))   
    ax2.xaxis.set_ticks_position('none')
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(axis='x', labelsize=fs*1.2)
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    # transparent background
    plt.savefig(outfile, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = outfile.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile} and {filename_tr}")



def make_log_reg_table(rank_summary_path_A, rank_summary_path_X):
    """
    make one table of A and X with the columns:
    * conservation rank
    * female biased (bin)
    * unbiased (bin)
    * male based (bin)
    * chromosome (A/X)
    """
    X_df = pd.read_csv(rank_summary_path_X, sep="\t")
    A_df = pd.read_csv(rank_summary_path_A, sep="\t")
    # add columns for chromosome cat
    X_df["chromosome"] = [1]*X_df.shape[0]
    A_df["chromosome"] = [0]*A_df.shape[0]

    df = pd.concat([A_df,X_df], ignore_index=True)
    print(df)

    ### ordinal logistic regression
    df['interaction'] = df['conservation_rank'] * df['chromosome']
    explanatory = df[['interaction','conservation_rank','chromosome']]
    response_abdomen = df['abdomen_sex_bias_category']
    response_head_thorax = df['head_thorax_sex_bias_category']

    print(f"\n////////////////// ABDOMEN //////////////////")
    model = OrderedModel(response_abdomen, explanatory, distr='logit')
    result = model.fit(method='bfgs', disp=False)
    print(result.summary())


    print(f"\n////////////////// HEAD+THORAX //////////////////")
    model = OrderedModel(response_head_thorax, explanatory, distr='logit')
    result = model.fit(method='bfgs', disp=False)
    print(result.summary())

############################


def make_phylogeny_rank_dict(summary_file_df, min_p):
    """
    make a dict with rank orders like { 1 : [list, of, Log2FC, numbers] ,  2 [more, log2FC, numbers] , ... }
    """
    filtered_df = summary_file_df.drop_duplicates(subset=['focal_transcript']) # remove duplicate data points
    
    LFC_dict_abdomen = { i : [] for i in range(1,6)}
    LFC_dict_head_thorax = { i : [] for i in range(1,6)}

    if min_p == 0:
        for p_rank, log2FC_val  in zip(filtered_df["level_most_dist_ortholog"], filtered_df["LFC_abdomen"]):
            LFC_dict_abdomen[p_rank].append(log2FC_val)
        for p_rank, log2FC_val  in zip(filtered_df["level_most_dist_ortholog"], filtered_df["LFC_head+thorax"]):
            LFC_dict_head_thorax[p_rank].append(log2FC_val)
    else:
        for p_rank, log2FC_val,pval  in zip(filtered_df["level_most_dist_ortholog"], filtered_df["LFC_abdomen"], filtered_df["FDR_pval_abdomen"]):
            if pval<min_p:
                LFC_dict_abdomen[p_rank].append(log2FC_val)
        for p_rank, log2FC_val,pval  in zip(filtered_df["level_most_dist_ortholog"], filtered_df["LFC_head+thorax"], filtered_df["FDR_pval_head+thorax"]):
            if pval<min_p:
                LFC_dict_head_thorax[p_rank].append(log2FC_val)

    return LFC_dict_abdomen, LFC_dict_head_thorax



def check_DE_phylogeny_rank_conserved(summary_paths_AX_list:dict, outfile = "", abs_LFC=False, sig_p_threshold = 0, sep_MF=True):
    """
    check if genes with a higher phylogeny conservarion rank have higher log2FC values 
    if abs_LFC=True then do abs() around LFC to assess general sex bias and don't differentiate male-female contrast
    """

    summary_data_A = pd.read_csv(summary_paths_AX_list["A"], sep = "\t", index_col=False)
    summary_data_X = pd.read_csv(summary_paths_AX_list["X"], sep = "\t", index_col=False)

    LFC_dict_abdomen_A, LFC_dict_head_thorax_A = make_phylogeny_rank_dict(summary_data_A, min_p = sig_p_threshold)
    LFC_dict_abdomen_X, LFC_dict_head_thorax_X = make_phylogeny_rank_dict(summary_data_X, min_p = sig_p_threshold)

    if abs_LFC:
        LFC_lists_abdomen_A = [[abs(val) for val in vals_list] for vals_list in LFC_dict_abdomen_A.values()]
        LFC_lists_head_thorax_A = [[abs(val) for val in vals_list] for vals_list in LFC_dict_head_thorax_A.values()]
        LFC_lists_abdomen_X = [[abs(val) for val in vals_list] for vals_list in LFC_dict_abdomen_X.values()]
        LFC_lists_head_thorax_X = [[abs(val) for val in vals_list] for vals_list in LFC_dict_head_thorax_X.values()]
        y_label = f"|log2FC|"
    else:
        LFC_lists_abdomen_A = [vals_list for vals_list in LFC_dict_abdomen_A.values()]
        LFC_lists_head_thorax_A = [vals_list for vals_list in LFC_dict_head_thorax_A.values()]
        LFC_lists_abdomen_X = [vals_list for vals_list in LFC_dict_abdomen_X.values()]
        LFC_lists_head_thorax_X = [vals_list for vals_list in LFC_dict_head_thorax_X.values()]
        y_label = f"log2FC for female-male"
    
    if sig_p_threshold>0:
        y_label = f"{y_label}, p<{sig_p_threshold}"


    #### plot 2x2 boxplots
    
    fs = 25 # font size

    # set figure aspect ratio
    aspect_ratio = 18 / 12
    height_pixels = 1200  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(2,2,figsize=(width_pixels / 100, height_pixels / 100), dpi=100)

    def plot_DE_subplot(lists, row,col, fs, title, colors_dict, abs_logFC, lw=2):
        tick_labels = [f"{i+1} ({len(vals)})" for i,vals in enumerate(lists)]
        bp = ax[row,col].boxplot(lists, patch_artist=True)   

        # set axis labels
        ax[row,col].tick_params(axis='x', labelsize=fs) 
        ax[row,col].set_xticks(ticks = range(1,6), labels = tick_labels, fontsize=fs*0.8)
        ax[row,col].tick_params(axis='y', labelsize=fs)
        ax[row,col].set_title(title, fontsize=fs)

        ## modify boxplot colors
        for box in bp['boxes']:
            box.set(facecolor=colors_dict["fill"], edgecolor=colors_dict["edge"], linewidth=2)
        for median in bp['medians']:
            median.set(color=colors_dict['medians'], linewidth=lw)
        for whisker in bp['whiskers']:
            whisker.set(color=colors_dict['edge'], linestyle='-',linewidth=lw)
        for cap in bp['caps']:
            cap.set(color=colors_dict['edge'],linewidth=lw)
        for flier in bp['fliers']:
            flier.set(marker='.', markerfacecolor=colors_dict['edge'], markeredgecolor=colors_dict['edge'])
    
        if abs_logFC == False:
            ax[row,col].axhline(y=0, linestyle = "--", color = "black")


    def plot_DE_subplot_sep_MF(lists, row,col, fs, title, colors_dict, abs_logFC, lw=2):
        
        ## split into male- and female biased for separate boxplots
        sex_biased_lists = []
        for lst in lists:
            male_biased = [x for x in lst if x < 0]
            female_biased = [x for x in lst if x > 0]
            sex_biased_lists.extend([male_biased, female_biased])

        # tick_labels = [f"{i // 2 + 1}\n({len(vals)})" for i,vals in enumerate(sex_biased_lists)]
        tick_labels = [f"({len(vals)})" for i,vals in enumerate(sex_biased_lists)] # plot conservation rank label on separate axis
        pos_adjust=0.2
        tick_pos = [i+pos_adjust if i%2==1 else i-pos_adjust for i in range(1,len(tick_labels)+1)]
        bp = ax[row,col].boxplot(sex_biased_lists, positions=tick_pos, patch_artist=True)

        # set axis labels
        tick_fs_factor = 0.6
        ax[row,col].set_xticks(ticks = tick_pos, labels = tick_labels, fontsize=fs*tick_fs_factor)
        ax[row,col].tick_params(axis='x', labelsize=fs*tick_fs_factor, rotation = 90)
        ax[row,col].tick_params(axis='y', labelsize=fs*0.9)
        ax[row,col].set_title(title, fontsize=fs)

        ax2 = ax[row,col].secondary_xaxis('bottom')
        ax2.set_xticks([i+0.5 for i in range(1,10,2)])
        ax2.set_xticklabels([i for i in range(1,6)], fontsize=fs)
        ax2.spines['bottom'].set_position(('outward', 40))   
        ax2.xaxis.set_ticks_position('none')
        ax2.spines['bottom'].set_visible(False)
        ax2.tick_params(axis='x', labelsize=fs)

        ## modify boxplot colors
        for i, box in enumerate(bp['boxes']):
            if i%2==0:
                box.set(facecolor=colors_dict["fill"], edgecolor=colors_dict["edge"], linewidth=2)
            else:
                box.set(facecolor=colors_dict["FB_fill"], edgecolor=colors_dict["FB_edge"], linewidth=2)
        for i, median in enumerate(bp['medians']):
            if i%2==0:
                median.set(color=colors_dict['medians'], linewidth=lw)
            else:
                median.set(color=colors_dict['FB_medians'], linewidth=lw)
        for i, whisker in enumerate(bp['whiskers']):
            # print(f"whisker: {i}")
            if i//2 % 2==0:
                whisker.set(color=colors_dict['edge'], linestyle='-',linewidth=lw)
            else:
                whisker.set(color=colors_dict['FB_edge'], linestyle='-',linewidth=lw)
        for i, cap in enumerate(bp['caps']):
            if i//2 % 2==0:
                cap.set(color=colors_dict['edge'],linewidth=lw)
            else:
                cap.set(color=colors_dict['FB_edge'],linewidth=lw)
        for i, flier in enumerate(bp['fliers']):
            if i%2==0:
                flier.set(marker='.', markerfacecolor=colors_dict['edge'], markeredgecolor=colors_dict['edge'])
            else:
                flier.set(marker='.', markerfacecolor=colors_dict['FB_edge'], markeredgecolor=colors_dict['FB_edge'])
    
        if abs_logFC == False:
            ax[row,col].axhline(y=0, linestyle = "--", color = "black")

    colors = {
        "fill" : "#495E83", # dusk blue
        "edge" : "#374C6E", # dusk blue darker
        "medians" : "#A7CCED", # icy blue
        "lines" : "#7D93B5", # lavender grey
        "FB_fill" : "#AB354A", # cherry rose
        "FB_edge" : "#771C2C", # dark amaranth
        "FB_medians" : "#EA9AA9", # cotton candy
        # "FB_lines" : "#7D93B5" # lavender grey
    }

    if sep_MF:
        print(f" separate male and female sex biased genes")
        plot_DE_subplot_sep_MF(LFC_lists_abdomen_A, row=0, col=0, fs=fs, title=f"Autosomes: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot_sep_MF(LFC_lists_head_thorax_A, row=0, col=1, fs=fs, title=f"Autosomes: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot_sep_MF(LFC_lists_abdomen_X, row=1, col=0, fs=fs, title=f"X-chromosome: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot_sep_MF(LFC_lists_head_thorax_X, row=1, col=1, fs=fs, title=f"X-chromosome: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        fig.supxlabel(f"(number of genes)\nconservation rank of C. maculatus ortholog", fontsize = fs)
    else:
        plot_DE_subplot(LFC_lists_abdomen_A, row=0, col=0, fs=fs, title=f"Autosomes: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot(LFC_lists_head_thorax_A, row=0, col=1, fs=fs, title=f"Autosomes: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot(LFC_lists_abdomen_X, row=1, col=0, fs=fs, title=f"X-chromosome: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot(LFC_lists_head_thorax_X, row=1, col=1, fs=fs, title=f"X-chromosome: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        fig.supxlabel(f"conservation rank, and (number of genes)", fontsize = fs)
    fig.supylabel(y_label, fontsize = fs)

    # layout (left, bottom, right, top)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = outfile.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile} and {filename_tr}")




if __name__ == "__main__":
    username = "miltr339"

    table_paths_dict = get_full_table_path(username=username)
    summary_table_paths = {}

    if False:
        for chromosome, path in table_paths_dict.items():
            # if chromosome=="X":
            #     continue
            print(f" --> {chromosome}")
            outfile_path = path.replace(".tsv", "_conservation_rank_analysis.tsv")
            # make_rank_summary_table(path, outfile_path=outfile_path, min_LFC=1, p_threshold=0.05)
            summary_table_paths[chromosome] = outfile_path


        plot_outfile_name=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DE_conservation_rank_proportions.png"
        make_category_proportion_plot(rank_summary_path_A=summary_table_paths["A"], 
                                  rank_summary_path_X=summary_table_paths["X"],
                                  outfile=plot_outfile_name)

        ### statistical analysis with logistic regression for binary categories
        make_log_reg_table(rank_summary_path_A=summary_table_paths["A"], 
                        rank_summary_path_X=summary_table_paths["X"])
    
    if True:
        summary_paths = get_summary_paths(username=username)
        abs_logFC = False
        check_DE_phylogeny_rank_conserved(summary_paths_AX_list=summary_paths,
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/conservation_rank_all_sex_bias_proportion.png",
            abs_LFC=abs_logFC, sig_p_threshold=0)
        check_DE_phylogeny_rank_conserved(summary_paths_AX_list=summary_paths,
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/conservation_rank_sig_sex_bias_proportion.png",
            abs_LFC=abs_logFC, sig_p_threshold=0.05)