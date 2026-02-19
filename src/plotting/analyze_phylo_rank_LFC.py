"""
Make a table for the phylogenetic rank analyses that has all the cmac genes in a row with the following columns
* geneID
* most distant conservation rank
* significant sex bias (1: female biased, -1: male biased, 0: unbiased) for
    * abdomen
    * head+thorax tissue 
"""


import pandas as pd

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

        
def make_category_proportion_plot(rank_summary_path_A:str, rank_summary_path_X:str, outfile_name:str):
    """
    plot a stacked bar to show the proportion of male/female/unbiased for X and A, and for abdominal and head+thorax
    """
    sum_dict = {
        "abdomen" : {
            "female_biased" : 0,
            "male_biased" : 0,
            "unbiased" : 0,
        },
        "head_thorax" : {
            "female_biased" : 0,
            "male_biased" : 0,
            "unbiased" : 0,
        },
    }
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
    
    def modify_counts_dict(df, sums_dict):
        unique_ranks = df["conservation_rank"].unique()
        rank_sums_dict = {}
        for rank in unique_ranks:
            rank_sums_dict[rank] = sums_dict

        for i, row in df.iterrows():
            row_rank = row["conservation_rank"]
            SB_int_a = row["abdomen_sex_bias_category"]
            SB_int_h = row["head_thorax_sex_bias_category"]
            rank_sums_dict[row_rank]["abdomen"][categorize_SB_int(SB_int_a)] += 1
            rank_sums_dict[row_rank]["head_thorax"][categorize_SB_int(SB_int_h)] += 1
        print(rank_sums_dict)
        return rank_sums_dict

    print(f"X rank sums dict")
    X_sums_dict = modify_counts_dict(X_df, sums_dict=sum_dict)
    print(f"A rank sums dict")
    A_sums_dict = modify_counts_dict(A_df, sums_dict=sum_dict)



if __name__ == "__main__":
    username = "milena"

    table_paths_dict = get_full_table_path(username=username)
    summary_table_paths = {}

    for chromosome, path in table_paths_dict.items():
        print(f" --> {chromosome}")
        outfile_path = path.replace(".tsv", "_conservation_rank_analysis.tsv")
        make_rank_summary_table(path, outfile_path=outfile_path, min_LFC=1, p_threshold=0.05)
        summary_table_paths[chromosome] = outfile_path

    plot_outfile_name=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DE_conservation_rank_proportions.png"
    make_category_proportion_plot(rank_summary_path_A=summary_table_paths["A"], 
                                  rank_summary_path_X=summary_table_paths["X"],
                                  outfile_name=plot_outfile_name)