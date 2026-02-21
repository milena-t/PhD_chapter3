"""
Investigate if the dNdS for different species pairs with Cmac is correlated with the log2FC in either tissue
use the conservation rank and A/X as fixed factors
"""

import pandas as pd

def get_full_table_path(username="miltr339"):
    out_dict = {
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_A_chr.tsv",
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_X_chr.tsv",
    }
    return out_dict


def reorder_df_to_have_dNdS_cols(full_table_paths_dict, outfile = ""):
    """
    split the single dNdS column into three: one for each Cmac species pair comparison
    This means that there is only one line for each unique Cmac transcripts
    """

    X_df = pd.read_csv(full_table_paths_dict["X"], sep="\t")
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    
    X_df["chromosome"] = [1]*X_df.shape[0]
    A_df["chromosome"] = [0]*A_df.shape[0]

    df = pd.concat([A_df,X_df], ignore_index=True)

    unique_df=df.drop_duplicates(subset=["focal_transcript"])
    unique_df=unique_df.drop(labels=["focal_species","other_species","other_transcript","positive_selection","dN","dS","dN/dS","focal_gene_ID"], axis=1)

    partners_list = list(set(df["other_species"]))

    for partner in partners_list:
        partner_df = df[df["other_species"]==partner]

        unique_df[f"{partner}_dN/dS"] = unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["dN/dS"])))
        unique_df[f"{partner}_dS"] = unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["dS"])))
        unique_df[f"{partner}_dN"] = unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["dN"])))
        unique_df[f"{partner}_positive_selection"] = unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["positive_selection"])))

    if outfile!="":
        unique_df.to_csv(outfile, index=False, sep="\t")
        print(f"rearranged data saved to file: {outfile}")

    return unique_df



def statistical_analysis(full_table_paths_dict, table_outfile=""):  

    df = reorder_df_to_have_dNdS_cols(full_table_paths_dict, outfile=table_outfile)
    ### chatgpt says to transform the dNdS values like this: df["log_dNdS"] = np.log(df["dNdS"] + 0.001)


if __name__ == "__main__":

    username = "milena"
    full_tables_dict = get_full_table_path(username=username)
    statistical_analysis(full_tables_dict, table_outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/paml_stats_outfile_table.tsv")