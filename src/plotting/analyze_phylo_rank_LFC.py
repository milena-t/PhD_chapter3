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
        elif log2FC < min_log2FC:
            return 0
        elif log2FC > 0: # positive -> female bias
            return 1
        elif min_log2FC < 0 : #negative -> male bias
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
        for i,row in full_summary_df.iterrows():
            
            # classify expression
            abd_sex_bias = classify_expression(log2FC = row["LFC_abdomen"], pvalue_FDR=row["LFC_abdomen"])
            ht_sex_bias = classify_expression(log2FC = row["LFC_head+thorax"], pvalue_FDR=row["FDR_pval_head+thorax"])
            # print(f"sex_bias a: {abd_sex_bias}, h+t: {ht_sex_bias}")
            
            # other vals
            trans_ID = row["focal_transcript"]
            cons_rank = row["level_most_dist_ortholog"]

            out_row = f"{trans_ID}\t{cons_rank}\t{abd_sex_bias}\t{ht_sex_bias}\n"
            outfile.write(out_row)
    print(f"done writing to: {outfile_path}!")

        

if __name__ == "__main__":
    username = "milena"

    table_paths_dict = get_full_table_path(username=username)

    for chromosome, path in table_paths_dict.items():
        print(f" --> {chromosome}")
        outfile_path = path.replace(".tsv", "_conservation_rank_analysis.tsv")
        make_rank_summary_table(path, outfile_path=outfile_path, min_LFC=1, p_threshold=0.05)