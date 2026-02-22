"""
Investigate if the dNdS for different species pairs with Cmac is correlated with the log2FC in either tissue
use the conservation rank and A/X as fixed factors
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

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

        unique_df[f"{partner}_dNdS"] = pd.to_numeric(unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["dN/dS"]))), errors='coerce')
        unique_df[f"{partner}_dS"] = pd.to_numeric(unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["dS"]))), errors='coerce')
        unique_df[f"{partner}_dN"] = pd.to_numeric(unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["dN"]))), errors='coerce')
        unique_df[f"{partner}_positive_selection"] = unique_df["focal_transcript"].map(dict(zip(partner_df["focal_transcript"], partner_df["positive_selection"])))

    if outfile!="":
        unique_df.to_csv(outfile, index=False, sep="\t")
        print(f"rearranged data saved to file: {outfile}")

    return unique_df



def statistical_analysis(full_table_paths_dict, table_outfile=""):  
    """
    fit linear regression to see if the response variable (dNdS) is explained by the log2FC in either tissue, with the fixed
    factors of sex chromosome category and conservation rank. NaNs are automatically dropped
    Also do separate tests for male and female-biased genes
    """
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    partners_list = list(set(A_df["other_species"]))
    print(partners_list)
    df = reorder_df_to_have_dNdS_cols(full_table_paths_dict, outfile=table_outfile)
    print(df)
    a_df_m = df[df["LFC_abdomen"] < 0 ]
    a_df_f = df[df["LFC_abdomen"] > 0 ]
    ht_df_m = df[df["LFC_head+thorax"] < 0 ]
    ht_df_f = df[df["LFC_head+thorax"] > 0 ]

    ## test for goodness of fit between the ols with log and glm with gamma distribution models
    ## test with B_siliquastri_dN/dS
    print(f"test goodness of fit, see code for model specification. OLS with log-transformed dNdS vsl GLM with gamma distribution")
    
    for partner in partners_list:
        filt_df = a_df_f[a_df_f[f"{partner}_dNdS"]>0]
        filt_df = filt_df[filt_df[f"{partner}_dNdS"].notna()]
        print(f"* {partner}")

        zeros_test = (a_df_f[f"{partner}_dNdS"] <= 0).sum()
        zeros_test2 = (filt_df[f"{partner}_dNdS"] <= 0).sum()
        print(f"\ttest if all dNdS=0 was removed: {zeros_test} -> {zeros_test2} (should go down to zero)")
        dNdS_nan_test = a_df_f[f"{partner}_dNdS"].isna().sum()
        dNdS_nan_test2 = filt_df[f"{partner}_dNdS"].isna().sum()
        print(f"\ttest if all dNdS=nan was removed: {dNdS_nan_test} -> {dNdS_nan_test2} (should go down to zero)")
        LFC_nan_test = a_df_f[f"LFC_abdomen"].isna().sum()
        LFC_nan_test2 = filt_df[f"LFC_abdomen"].isna().sum()
        print(f"\ttest if all logFC=nan was removed: {LFC_nan_test} -> {LFC_nan_test2} (should go down to zero)")

        ols = smf.ols("np.log(B_siliquastri_dNdS) ~ LFC_abdomen + C(chromosome) + level_most_dist_ortholog", data=filt_df).fit()
        print(f"  - ols : {ols.aic}")
        gamma = smf.glm(
            "B_siliquastri_dNdS ~ LFC_abdomen + C(chromosome) + level_most_dist_ortholog",
            data=filt_df,
            family=sm.families.Gamma(link=sm.families.links.Log())
        ).fit()
        print(f"  - gamma : {gamma.aic}")
        print()
    




if __name__ == "__main__":

    username = "milena"
    full_tables_dict = get_full_table_path(username=username)
    statistical_analysis(full_tables_dict, table_outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/paml_stats_outfile_table.tsv")