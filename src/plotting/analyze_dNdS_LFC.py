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




def make_phylogeny_rank_dict(summary_file_df, focal_species, max_dNdS=2):
    """
    make a dict with rank orders like { 1 : [list, of, dNdS, numbers] ,  2 [more, dNdS, numbers] , ... }
    """
    filtered_df = summary_file_df[summary_file_df["other_species"] == focal_species]
    if max_dNdS >0:
        filtered_df = filtered_df[filtered_df["dN/dS"]<max_dNdS]
    filtered_df["dN/dS"] = pd.to_numeric(filtered_df["dN/dS"], errors='coerce')
    filtered_df = filtered_df.dropna(subset="dN/dS")
    
    dNdS_dict = { i : [] for i in range(1,6)}

    for p_rank, log2FC_val  in zip(filtered_df["level_most_dist_ortholog"], filtered_df["dN/dS"]):
        dNdS_dict[p_rank].append(log2FC_val)

    return dNdS_dict


def plot_dNdS_rank_conserved(summary_paths_AX_list:dict, outfile = "", maxdNdS = 0):
    """
    check if genes with a higher phylogeny conservarion rank have higher log2FC values 
    if abs_LFC=True then do abs() around LFC to assess general sex bias and don't differentiate male-female contrast
    """

    summary_data_A = pd.read_csv(summary_paths_AX_list["A"], sep = "\t", index_col=False)
    summary_data_X = pd.read_csv(summary_paths_AX_list["X"], sep = "\t", index_col=False)
    partner_species = list(set(summary_data_X["other_species"]))
    partner_lists_A = {partner : [] for partner in partner_species}
    partner_lists_X = {partner : [] for partner in partner_species}

    for partner in partner_species:
        partner_lists_A[partner] = make_phylogeny_rank_dict(summary_data_A, focal_species=partner, max_dNdS=maxdNdS)
        partner_lists_X[partner] = make_phylogeny_rank_dict(summary_data_X, focal_species=partner, max_dNdS=maxdNdS)

    
    y_label = f"dNdS"
    if maxdNdS>0:
        y_label = f"dNdS (max. < {maxdNdS})"
    
    fs = 30 # font size

    # set figure aspect ratio
    aspect_ratio = 12 / 8
    height_pixels = 1000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels


    def plot_dNdS_subplot(lists_A, lists_X, fs, title, colors_dict, lw=2):
        
        ## make A and X lists alternating to plot
        AX_lists = []
        for rank in lists_A.keys():
            AX_lists.extend([lists_A[rank], lists_X[rank]])
        
        # tick_labels = [f"{i // 2 + 1}\n({len(vals)})" for i,vals in enumerate(sex_biased_lists)]
        tick_labels = [f"({len(vals)})" for i,vals in enumerate(AX_lists)] # plot conservation rank label on separate axis
        pos_adjust=0.2
        tick_pos = [i+pos_adjust if i%2==1 else i-pos_adjust for i in range(1,len(tick_labels)+1)]
        bp = ax.boxplot(AX_lists, positions=tick_pos, patch_artist=True)

        # set axis labels
        tick_fs_factor = 0.75
        ax.set_xticks(ticks = tick_pos, labels = tick_labels, fontsize=fs*tick_fs_factor, rotation=45, ha='right')
        ax.tick_params(axis='x', labelsize=fs*tick_fs_factor)#, rotation = 90)
        ax.tick_params(axis='y', labelsize=fs*0.9)
        ax.set_title(title, fontsize=fs)

        ax2 = ax.secondary_xaxis('bottom')
        ax2.set_xticks([i+0.5 for i in range(1,10,2)])
        ax2.set_xticklabels([i for i in range(1,6)], fontsize=fs)
        ax2.spines['bottom'].set_position(('outward', 60))   
        ax2.xaxis.set_ticks_position('none')
        ax2.spines['bottom'].set_visible(False)
        ax2.tick_params(axis='x', labelsize=fs)

        ## modify boxplot colors
        for i, box in enumerate(bp['boxes']):
            if i%2==0:
                box.set(facecolor=colors_dict["fill"], edgecolor=colors_dict["edge"], linewidth=2)
            else:
                box.set(facecolor=colors_dict["X_fill"], edgecolor=colors_dict["X_edge"], linewidth=2)
        for i, median in enumerate(bp['medians']):
            if i%2==0:
                median.set(color=colors_dict['medians'], linewidth=lw)
            else:
                median.set(color=colors_dict['X_medians'], linewidth=lw)
        for i, whisker in enumerate(bp['whiskers']):
            # print(f"whisker: {i}")
            if i//2 % 2==0:
                whisker.set(color=colors_dict['edge'], linestyle='-',linewidth=lw)
            else:
                whisker.set(color=colors_dict['X_edge'], linestyle='-',linewidth=lw)
        for i, cap in enumerate(bp['caps']):
            if i//2 % 2==0:
                cap.set(color=colors_dict['edge'],linewidth=lw)
            else:
                cap.set(color=colors_dict['X_edge'],linewidth=lw)
        for i, flier in enumerate(bp['fliers']):
            if i%2==0:
                flier.set(marker='.', markerfacecolor=colors_dict['edge'], markeredgecolor=colors_dict['edge'])
            else:
                flier.set(marker='.', markerfacecolor=colors_dict['X_edge'], markeredgecolor=colors_dict['X_edge'])

    colors = {
        "fill" : "#246A73", # stormy teal
        "edge" : "#174C54", # dark teal
        "medians" : "#54A6A2", # tropical teal
        # "lines" : "#7D93B5", # lavender grey
        "X_fill" : "#7E3A7E", # grape soda
        "X_edge" : "#672B67", # velvet purple
        "X_medians" : "#C877C8", # orchid mist
        # "X_lines" : "#7D93B5" # lavender grey
    }


    for partner in partner_species:
        fig, ax = plt.subplots(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
        print(f"/////////////////////// {partner} ///////////////////////")
        partner_title=partner.replace("_", ". ")
        plot_title = f"C. maculatus vs. {partner_title} pairwise comparison"
        plot_dNdS_subplot(lists_A=partner_lists_A[partner], lists_X=partner_lists_X[partner], fs=fs, title=plot_title, colors_dict=colors)
        fig.supxlabel(f"(number of genes)\nconservation rank of C. maculatus ortholog", fontsize = fs)
        fig.supylabel(y_label, fontsize = fs, x=0.0, y=0.625)

        # layout (left, bottom, right, top)
        plt.tight_layout(rect=[0.0, 0.05, 1, 1])

        outfile_partner = outfile.replace(".png", f"{partner}.png")
        # transparent background
        plt.savefig(outfile_partner, dpi = 300, transparent = True)
        # non-transparent background
        filename_tr = outfile_partner.replace(".png", "_white_bg.png")
        plt.savefig(filename_tr, dpi = 300, transparent = False)
        print(f"plot saved in current working directory as: {outfile_partner} and {filename_tr}")





if __name__ == "__main__":

    username = "miltr339"
    full_tables_dict = get_full_table_path(username=username)
    reorg_table_outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/paml_stats_outfile_table.tsv"
    
    if False:
        ###################################################
        ## median quantile regression for dNdS
        statistical_analysis_dNdS(full_tables_dict, table_outfile=f"")
        ###################################################

    if False:
        ###################################################
        ## logistic regression for categorical response (positive selection True/False)
        statistical_analysis_pos_sel(full_table_paths_dict=full_tables_dict)
        ###################################################

    if True:
        ###################################################
        ## plotting 
        filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_vs_conservation_rank.png"
        plot_dNdS_rank_conserved(summary_paths_AX_list= full_tables_dict, outfile = filename, maxdNdS = 0)