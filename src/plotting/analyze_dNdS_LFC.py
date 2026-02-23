"""
Investigate if the dNdS for different species pairs with Cmac is correlated with the log2FC in either tissue
use the conservation rank and A/X as fixed factors
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import patsy

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
    
    X_df["chromosome"] = ["X"]*X_df.shape[0]
    A_df["chromosome"] = ["A"]*A_df.shape[0]

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



def statistical_analysis_dNdS(full_table_paths_dict, table_outfile="", max_dNdS=2):  
    """
    fit linear regression to see if the response variable (dNdS) is explained by the log2FC in either tissue, with the fixed
    factors of sex chromosome category and conservation rank. NaNs are automatically dropped
    Also do separate tests for male and female-biased genes
    """
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    partners_list = list(set(A_df["other_species"]))
    print(partners_list)

    reorder_dNdS_partners = False
    if reorder_dNdS_partners:
        df = reorder_df_to_have_dNdS_cols(full_table_paths_dict, outfile=table_outfile)
    else:
        X_df = pd.read_csv(full_table_paths_dict["X"], sep="\t")
        A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")

        X_df["chromosome"] = ["X"]*X_df.shape[0]
        A_df["chromosome"] = ["A"]*A_df.shape[0]

        df = pd.concat([A_df,X_df], ignore_index=True)

    df = df.rename(columns={'LFC_head+thorax': 'LFC_head_thorax'})
    dfs_dict = { 
        "female_biased_abdomen": df[df["LFC_abdomen"] > 0 ],
        "male_biased_abdomen" : df[df["LFC_abdomen"] < 0  ],
        "female_biased_head_thorax": df[df["LFC_head_thorax"] > 0 ],
        "male_biased_head_thorax" : df[df["LFC_head_thorax"] < 0 ]
    }
    ## test for goodness of fit between the ols with log and glm with gamma distribution models
    ## test with B_siliquastri_dN/dS
    
    for partner in partners_list:
        if reorder_dNdS_partners:
            filt_df = df[df[f"{partner}_dNdS"]>0]
            filt_df = filt_df[filt_df[f"{partner}_dNdS"].notna()]
            filt_df = filt_df[filt_df[f"{partner}_dNdS"] < max_dNdS]
        else:
            df_ren = df.rename(columns={'dN/dS': f"{partner}_dNdS"})
            filt_df = df_ren[df_ren["other_species"]==partner]
            filt_df[f"{partner}_dNdS"] = pd.to_numeric(filt_df[f"{partner}_dNdS"], errors='coerce')
            filt_df = filt_df[filt_df[f"{partner}_dNdS"]>0]
            filt_df = filt_df[filt_df[f"{partner}_dNdS"].notna()]
            filt_df = filt_df[filt_df[f"{partner}_dNdS"] < max_dNdS]
        print(f"\n////////////////// {partner} //////////////////")

        if False:
            ### test filtering stuff and see if there is anything that could make problems
            def test_filtering(colname, nan_only=False):
                print(f"   {colname} filtering:")
                if nan_only:
                    pass
                else:
                    zeros_test = (a_df_f[f"{colname}"] <= 0).sum()
                    zeros_test2 = (filt_df[f"{colname}"] <= 0).sum()
                    print(f"\ttest if all 0 was removed: {zeros_test} -> {zeros_test2} (should go down to zero)")
                nan_test = a_df_f[f"{colname}"].isna().sum()
                nan_test2 = filt_df[f"{colname}"].isna().sum()
                print(f"\ttest if all nan was removed: {nan_test} -> {nan_test2} (should go down to zero)")
                infinity = np.isfinite(filt_df[f"{colname}"]).all()
                print(f"\tall values <inf :  {infinity}")

            test_filtering(f"{partner}_dNdS")
            test_filtering(f"LFC_abdomen")
            test_filtering(f"chromosome", nan_only=True)
            test_filtering(f"level_most_dist_ortholog")

        if True:
            formula = f"{partner}_dNdS ~ (LFC_abdomen + LFC_head_thorax) * C(chromosome) * level_most_dist_ortholog"
                ## the syntax with the parentheses (LFC_abdomen + LFC_head_thorax) * C(chromosome) means this:
                # LFC_abdomen + LFC_head_thorax + C(chromosome) + LFC_abdomen:C(chromosome) + LFC_head_thorax:C(chromosome)
            test = smf.quantreg(formula=formula, data=filt_df).fit(q=0.5) # q=0.5 means we estimate the median
            print(test.summary())
            model="lreg"

        else:
            ### ordinary least-squasres does not have normal residuals!
            formula = f"np.log({partner}_dNdS) ~ (LFC_abdomen + LFC_head_thorax) * C(chromosome) + level_most_dist_ortholog"
            test = smf.ols(formula=formula, data=filt_df).fit(cov_type="HC3")
            model = "OLS"
            print(test.summary())
            res = stats.normaltest(test.resid)
            if False:
                plt.scatter(x=test.fittedvalues, y=test.resid)
                plt.axhline(0, color="red")
                plt.ylabel("residuals")
                plt.xlabel("fitted values")
                plt.suptitle(f"C. maculatus and {partner} ")
                plt.show()

            ## GLM does not work well
            formula = f"{partner}_dNdS ~ (LFC_abdomen + LFC_head_thorax) * C(chromosome) + level_most_dist_ortholog"
            test = smf.glm(formula=formula,data=filt_df,family=sm.families.Gamma(link=sm.families.links.Log())).fit()
            model = "GLM"
            print(test.summary())
            res = stats.normaltest(test.resid_response)

            if False:
                # plot residuals distribution
                # test normal distribution of residuals
                if res.pvalue <0.05:
                    res_nom = f"{model} residuals NOT normally distributed"
                else:
                    res_nom = f"{model} residuals normally distributed!"
                if model=="OLS":
                    plt.hist(res)
                    plt.suptitle(f"C. maculatus and {partner} dNdS {model} regression residuals:\n{res_nom}")
                else:
                    # here is the reason why GLM does not work well
                    plt.scatter(test.fittedvalues, test.resid_deviance)
                    plt.suptitle(f"C. maculatus and {partner} dNdS {model} regression\npredicted vs. observed values")
                plt.show()

            elif model == "OLS":
                # test normal distribution of residuals
                if res.pvalue <0.05:
                    print(f"{model} residuals NOT normally distributed")
                else:
                    print(f"{model} residuals normally distributed!")

        print()
    

def statistical_analysis_pos_sel(full_table_paths_dict):
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    partners_list = list(set(A_df["other_species"]))
    print(partners_list)

    X_df = pd.read_csv(full_table_paths_dict["X"], sep="\t")
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    X_df["chromosome"] = ["X"]*X_df.shape[0]
    A_df["chromosome"] = ["A"]*A_df.shape[0]
    df = pd.concat([A_df,X_df], ignore_index=True)
    df = df.rename(columns={'LFC_head+thorax': 'LFC_head_thorax'})

    for partner in partners_list:
        print(f"\n////////////////// {partner} //////////////////")
        filt_df = df[df["other_species"]==partner]
        filt_df = filt_df[filt_df["positive_selection"] != "NaN"]
        filt_df = filt_df[filt_df["positive_selection"].notna()]
        filt_df["positive_selection"] = filt_df["positive_selection"].map({False: 0, True: 1})
        # print(filt_df["positive_selection"].unique(), filt_df["positive_selection"].dtype)

        formula = f"positive_selection ~ (LFC_abdomen + LFC_head_thorax) * C(chromosome) * level_most_dist_ortholog"

        test = smf.logit(formula=formula, data=filt_df).fit()
        print(test.summary())




if __name__ == "__main__":

    username = "miltr339"
    full_tables_dict = get_full_table_path(username=username)
    reorg_table_outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/paml_stats_outfile_table.tsv"
    
    ## median quantile regression for dNdS
    # statistical_analysis_dNdS(full_tables_dict, table_outfile=f"")

    ## logistic regression for categorical response (positive selection True/False)
    statistical_analysis_pos_sel(full_table_paths_dict=full_tables_dict)


