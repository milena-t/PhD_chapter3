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

        
def make_category_proportion_plot(rank_summary_path_A:str, rank_summary_path_X:str, outfile:str):
    """
    plot a stacked bar to show the proportion of male/female/unbiased for X and A, and for abdominal and head+thorax
    """
    sum_dict = {"abdomen" : {"female_biased" : 0,"male_biased" : 0,"unbiased" : 0},"head_thorax" : {"female_biased" : 0,"male_biased" : 0,"unbiased" : 0}}
    X_df = pd.read_csv(rank_summary_path_X, sep="\t")
    A_df = pd.read_csv(rank_summary_path_A, sep="\t")

    def categorize_SB_int(SB_int:int):
        if SB_int==0:
            return "unbiased"
        elif SB_int==-1:
            return "male_biased"
        elif SB_int==1:
            return "female_biased"
        else:
            raise RuntimeError(f"{SB_int} could not be categorized")
    
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
        
        for row_rank in perc_dict.keys():
            num_genes = df[df["conservation_rank"]==row_rank].shape[0]
            print(f"{row_rank} : {num_genes} genes")
            perc_dict[row_rank]["abdomen"] = { sex_bias : 100*float(val)/float(num_genes) for sex_bias,val in counts_dict[row_rank]["abdomen"].items()}
            perc_dict[row_rank]["head_thorax"] = { sex_bias : 100*float(val)/float(num_genes) for sex_bias,val in counts_dict[row_rank]["head_thorax"].items()}
        
        print(perc_dict)
        return perc_dict


    X_sums_dict = modify_counts_dict(X_df)
    print(f"\n ---> X rank percentages dict")
    X_perc_dict = make_percentage_dict(X_sums_dict, X_df)
    A_sums_dict = modify_counts_dict(A_df)
    print(f"\n ---> A rank percentages dict")
    A_perc_dict = make_percentage_dict(A_sums_dict, A_df)
    assert len(X_perc_dict) == len(A_perc_dict)

    ### plot stacked bar chart
    fs = 25
    width=0.2
    x_gap = 0.025
    # x_subtr=width*1.1/2.0
    x_coords = (list(range(1,len(X_sums_dict.keys())+1)))
    x_coords_bars = [
        value
        for x in x_coords
        for value in (
            x - width*1.5 - x_gap,
            x - width*0.5 - x_gap,
            x + width*0.5 - x_gap,
            x + width*1.5 - x_gap,
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
            "male_biased" : "#5E5F7A", # wisteria blue
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

    bottom = [0,0,0,0]*4
    subcat_tissue = ["abdomen","head_thorax"]
    chr_labels = [f"{chr}:{tissue}" for tissue in subcat_tissue for chr in ["A","X"]]
    bar_labels = [f"{conserved_rank}:{sublabel}" for sublabel in chr_labels for conserved_rank in x_coords]

    # stack bars from bottom to top
    for sex_bias in colors_dict["abdomen"].keys():
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
        assert len(x_coords_bars)==len(y_coords)

        ax.bar(x_coords_bars, y_coords, width = width, bottom=bottom, color= colors_list)
        
        bottom = [bottom[i]+y_coords[i] for i in range(len(y_coords))]

    # ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<1 else f'{int(x)}%'))
    ax.set_xticks(ticks=x_coords, labels=bar_labels)
    ax.tick_params(axis='y', labelsize=fs)
    ax.tick_params(axis='x', labelsize=fs, rotation=90)
    ax.set_xticks(ticks=x_coords_bars, labels=bar_labels)
    
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

    for chromosome, path in table_paths_dict.items():
        # if chromosome=="X":
        #     continue
        print(f" --> {chromosome}")
        outfile_path = path.replace(".tsv", "_conservation_rank_analysis.tsv")
        make_rank_summary_table(path, outfile_path=outfile_path, min_LFC=1, p_threshold=0.05)
        summary_table_paths[chromosome] = outfile_path

    plot_outfile_name=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DE_conservation_rank_proportions.png"
    make_category_proportion_plot(rank_summary_path_A=summary_table_paths["A"], 
                                  rank_summary_path_X=summary_table_paths["X"],
                                  outfile=plot_outfile_name)