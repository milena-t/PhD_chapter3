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
from matplotlib.ticker import FuncFormatter

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


def make_sex_bias_cat_row(row, tissue = "abdomen"):
    if row[f"LFC_{tissue}"] < -1:
        if row[f"FDR_pval_{tissue}"]<0.05:
            return "male"
        else:
            return "unbiased"
    elif row[f"LFC_{tissue}"] >1:
        if row[f"FDR_pval_{tissue}"]<0.05:
            return "female"
        else:
            return "unbiased"
    else:
        return "unbiased"


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
    df = df.rename(columns={'FDR_pval_head+thorax': 'FDR_pval_head_thorax'})

    df["SB_abdomen"] = df.apply(make_sex_bias_cat_row, axis=1, args=("abdomen",))
    df["SB_head_thorax"] = df.apply(make_sex_bias_cat_row, axis=1, args=("head_thorax",))
    df["LFC_abdomen"] = abs(df["LFC_abdomen"])
    df["LFC_head_thorax"] = abs(df["LFC_head_thorax"])
    print(df)
    
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
                
            if "C_chinensis" in partner:
                # only significantly sex-biased genes -> remove LFC
                formula_a = f"{partner}_dNdS ~  C(SB_abdomen)  * C(chromosome) * level_most_dist_ortholog"
                formula_ht = f"{partner}_dNdS ~  C(SB_head_thorax)  * C(chromosome) * level_most_dist_ortholog"
                
                print(f"\n------------> abdomen")
                test = smf.quantreg(formula=formula_a, data=filt_df).fit(q=0.5) # q=0.5 means we estimate the median
                print(test.summary())
                print(f"\n------------> head+thorax")
                test = smf.quantreg(formula=formula_ht, data=filt_df).fit(q=0.5) # q=0.5 means we estimate the median
                print(test.summary())
            else:
                formula = f"{partner}_dNdS ~  C(chromosome) * level_most_dist_ortholog"
                # the syntax with the parentheses (LFC_abdomen + LFC_head_thorax) * C(chromosome) means this:
                # LFC_abdomen + LFC_head_thorax + C(chromosome) + LFC_abdomen:C(chromosome) + LFC_head_thorax:C(chromosome)
 
                test = smf.quantreg(formula=formula, data=filt_df).fit(q=0.5) # q=0.5 means we estimate the median
                print(test.summary())
            
                try:
                    # test three-way interactions relevance with 
                    interactions_test_string = """
                    LFC_head_thorax:C(chromosome)[T.X]:level_most_dist_ortholog = 0,
                    LFC_abdomen:C(chromosome)[T.X]:level_most_dist_ortholog = 0,
                    C(SB_abdomen)[T.male]:C(chromosome)[T.X]:level_most_dist_ortholog = 0,
                    C(SB_head_thorax)[T.male]:C(chromosome)[T.X]:level_most_dist_ortholog = 0
                    """
                    wald_test = test.wald_test(interactions_test_string, scalar = True)
                    print(f"wald test for all three-way interactions: {wald_test}")
                except:
                    print("no Wald test could be performed")

                print(f"\n---------------> test without conservation rank to see if excluding it makes chromosome significant")
                # do one test without age rank to see if chromosome becomes significant to explain the results from the permutation test
                formula = f"{partner}_dNdS ~  C(chromosome)"
                test = smf.quantreg(formula=formula, data=filt_df).fit(q=0.5) # q=0.5 means we estimate the median
                print(test.summary())
            print("\n")

    

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
    df = df.rename(columns={'FDR_pval_head+thorax': 'FDR_pval_head_thorax'})
    df["SB_abdomen"] = df.apply(make_sex_bias_cat_row, axis=1, args=("abdomen",))
    df["SB_head_thorax"] = df.apply(make_sex_bias_cat_row, axis=1, args=("head_thorax",))
    df["LFC_abdomen"] = abs(df["LFC_abdomen"])
    df["LFC_head_thorax"] = abs(df["LFC_head_thorax"])

    for partner in partners_list:
        print(f"\n////////////////// {partner} //////////////////")
        filt_df = df[df["other_species"]==partner]
        filt_df = filt_df[filt_df["positive_selection"] != "NaN"]
        filt_df = filt_df[filt_df["positive_selection"].notna()]
        filt_df["positive_selection"] = filt_df["positive_selection"].map({False: 0, True: 1})
        
        if True:
            formula = f"positive_selection ~ (LFC_abdomen + LFC_head_thorax + C(SB_abdomen) + C(SB_head_thorax)) * C(chromosome) * level_most_dist_ortholog"

            # when adjusting to the same as contninuous dNdS -> remove LFC and merge significant and nonsignificant into one category where there are three factor levels

            # the syntax with the parentheses (LFC_abdomen + LFC_head_thorax) * C(chromosome) means this:
            # LFC_abdomen + LFC_head_thorax + C(chromosome) + LFC_abdomen:C(chromosome) + LFC_head_thorax:C(chromosome)
            if "C_chinensis" in partner:
                # only significantly sex-biased genes -> remove LFC
                formula_a = f"positive_selection ~  C(SB_abdomen)  * C(chromosome) * level_most_dist_ortholog"
                formula_ht = f"positive_selection ~  C(SB_head_thorax)  * (C(chromosome) + level_most_dist_ortholog)"
                formula_no = f"positive_selection ~  C(chromosome) * level_most_dist_ortholog"
                formula_nono = f"positive_selection ~  C(chromosome)"
                
                print(f"\n------------> abdomen")
                test = smf.logit(formula=formula_a, data=filt_df).fit()
                print(test.summary())
                print(f"\n------------> head+thorax")
                test = smf.logit(formula=formula_ht, data=filt_df).fit()
                print(test.summary())
                print(f"\n------------> no sex bias")
                test = smf.logit(formula=formula_no, data=filt_df).fit()
                print(test.summary())
                print(f"\n------------> no conservation rank")
                test = smf.logit(formula=formula_nono, data=filt_df).fit()
                print(test.summary())
            else:
                
                formula = f"positive_selection ~  C(chromosome) * level_most_dist_ortholog"
                test = smf.logit(formula=formula, data=filt_df).fit()
                print(test.summary())
                
                print(f"\n---------------> test without conservation rank to see if excluding it makes chromosome significant")
                # do one test without age rank to see if chromosome becomes significant to explain the results from the permutation test
                formula = f"positive_selection ~  C(chromosome)"
                test = smf.logit(formula=formula, data=filt_df).fit()
                print(test.summary())

            try:
                interactions_test_string = """
                LFC_head_thorax:C(chromosome)[T.X]:level_most_dist_ortholog = 0,
                LFC_abdomen:C(chromosome)[T.X]:level_most_dist_ortholog = 0,
                C(SB_abdomen)[T.male]:C(chromosome)[T.X]:level_most_dist_ortholog = 0,
                C(SB_head_thorax)[T.male]:C(chromosome)[T.X]:level_most_dist_ortholog = 0
                """
                wald_test = test.wald_test(interactions_test_string, scalar = True)
                print(f"wald test: {wald_test}")
            except:
                pass





def make_phylogeny_rank_dict(summary_file_df, focal_species="", max_dNdS=2):
    """
    make a dict with rank orders like { 1 : [list, of, dNdS, numbers] ,  2 [more, dNdS, numbers] , ... }
    """
    if focal_species!="":
        filtered_df = summary_file_df[summary_file_df["other_species"] == focal_species]
        if max_dNdS >0:
            filtered_df = filtered_df[filtered_df["dN/dS"]<max_dNdS]
        filtered_df["dN/dS"] = pd.to_numeric(filtered_df["dN/dS"], errors='coerce')
        filtered_df = filtered_df.dropna(subset="dN/dS")
    else:
        filtered_df = summary_file_df
    
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


def compare_conservation_rank_proportions(summary_paths_AX):

    summary_data_A = pd.read_csv(summary_paths_AX["A"], sep = "\t", index_col=False)
    summary_data_X = pd.read_csv(summary_paths_AX["X"], sep = "\t", index_col=False)
    unique_A=summary_data_A.drop_duplicates(subset=["focal_transcript"])
    unique_X=summary_data_X.drop_duplicates(subset=["focal_transcript"])
    total_A = unique_A.shape[0]
    total_X = unique_X.shape[0]

    counts_A = {rank : f"{100*len(list)/total_A:.3f}%" for rank, list in make_phylogeny_rank_dict(unique_A, focal_species="", max_dNdS=0).items()}
    counts_X = {rank : f"{100*len(list)/total_X:.3f}%" for rank, list in make_phylogeny_rank_dict(unique_X, focal_species="", max_dNdS=0).items()}

    print(f"count proportions of conservation ranks in A and X:")
    for rank in counts_A.keys():
        print(f"{rank} -->\tA = {counts_A[rank]}\t X = {counts_X[rank]}")



def make_phylogeny_rank_nested_dict(summary_file_df, tissue, max_dNdS=2, pos_sel = False):
    """
    make a dict with rank orders like 
    { 1 : {
        "A" : { "male" : [list, of, dNdS, numbers] ,
                "female" : [more, dNdS, numbers] , 
                "unbiased" : [...] },
        "X" : { "male" : [list, of, dNdS, numbers] ,
                "female" : [more, dNdS, numbers] , 
                "unbiased" : [...] },
        },
      2 : { ... },
    ...
    }
    if pos_sel then the list is True/False according to the positive_selection column
    """

    summary_file_df["dN/dS"] = pd.to_numeric(summary_file_df["dN/dS"], errors='coerce')
    filtered_df = summary_file_df.dropna(subset="dN/dS")
    if max_dNdS>0:
        filtered_df = filtered_df[filtered_df["dN/dS"] < max_dNdS]
    
    dNdS_dict = {
    i: { "A": {"male": [], "female": [], "unbiased": []},
         "X": {"male": [], "female": [], "unbiased": []} }
    for i in range(1, 6)}

    if pos_sel:
        for p_rank, SB_cat, possel, chr  in zip(filtered_df["level_most_dist_ortholog"], filtered_df[f"SB_{tissue}"],filtered_df["positive_selection"],filtered_df["chromosome"]):
            dNdS_dict[p_rank][chr][SB_cat].append(possel)
    else:    
        for p_rank, SB_cat, dNdS, chr  in zip(filtered_df["level_most_dist_ortholog"], filtered_df[f"SB_{tissue}"],filtered_df["dN/dS"],filtered_df["chromosome"]):
            dNdS_dict[p_rank][chr][SB_cat].append(dNdS)

    return dNdS_dict


def boxplot_dNdS(full_table_paths_dict, outfile, maxdNdS=2, partner_species="C_chinensis", pos_sel = False):

    X_df = pd.read_csv(full_table_paths_dict["X"], sep="\t")
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    X_df["chromosome"] = ["X"]*X_df.shape[0]
    A_df["chromosome"] = ["A"]*A_df.shape[0]
    df = pd.concat([A_df,X_df], ignore_index=True)

    df = df.rename(columns={'LFC_head+thorax': 'LFC_head_thorax'})
    df = df.rename(columns={'FDR_pval_head+thorax': 'FDR_pval_head_thorax'})
    print(df)
    df["SB_abdomen"] = df.apply(make_sex_bias_cat_row, axis=1, args=("abdomen",))
    df["SB_head_thorax"] = df.apply(make_sex_bias_cat_row, axis=1, args=("head_thorax",))
    df = df[df["other_species"]==partner_species]

    if pos_sel:
        nested_vals_dict_a = make_phylogeny_rank_nested_dict(df, tissue="abdomen", max_dNdS=maxdNdS, pos_sel=True)
        nested_vals_dict_ht = make_phylogeny_rank_nested_dict(df, tissue="head_thorax", max_dNdS=maxdNdS, pos_sel=True)
    else:
        nested_vals_dict_a = make_phylogeny_rank_nested_dict(df, tissue="abdomen", max_dNdS=maxdNdS)
        nested_vals_dict_ht = make_phylogeny_rank_nested_dict(df, tissue="head_thorax", max_dNdS=maxdNdS)
    
    y_label = f"dN/dS"
    if maxdNdS>0:
        y_label = f"dN/dS (max. {maxdNdS})"
    
    fs = 30 # font size

    # set figure aspect ratio
    aspect_ratio = 12 / 8
    height_pixels = 1000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    def plot_dNdS_subplot(AX_dicts, fs, title, colors_dict, lw=2, outfile = "plot.png", pos_sel = False):

        fig, ax = plt.subplots(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
        SB_order_list =["male","unbiased","female"]
        
        color_faces = []
        color_edges = []
        color_medians = []

        for SB_cat in SB_order_list:
            color_faces.append(colors_dict[SB_cat]["fill"])
            color_edges.append(colors_dict[SB_cat]["edge"])
            color_medians.append(colors_dict[SB_cat]["medians"])
        
        if pos_sel:
            
            ## make tick labels in the same order as the box plot lists to be sure that everyting is right
            tick_labels = []
            AX_lists_pos = []
            AX_lists_neg = []
            cols_list_sig = []
            cols_list_unsig = []

            print(f"{title}")
            for i in range(1,6): # hard-code rank order. dicts are ordered so it should be right by default but just to be sure
                chr_dict = AX_dicts[i]
                print(f"rank {i}")
                for chr in ["A","X"]:
                    print(f"\t - {chr}")
                    SB_dict = chr_dict[chr]
                    for SB_cat in SB_order_list:

                        ### to double check plot coloring and stuff right: plot verbose tick labels!
                        # num_genes =f"({len(SB_dict[SB_cat])}):{chr}:{SB_cat}"
                        ### 

                        num_genes =len(SB_dict[SB_cat])
                        print(f"\t\t - {SB_cat} has {num_genes} genes")
                        tick_labels.append(f"({num_genes})")

                        A_nonsig = len([val for val in SB_dict[SB_cat] if val == False])*100.0
                        A_sig = len([val for val in SB_dict[SB_cat] if val == True])*100.0
                        try:
                            AX_lists_pos.append(A_sig/num_genes)
                        except:
                            AX_lists_pos.append(0)
                        try:
                            AX_lists_neg.append(A_nonsig/num_genes)
                        except:
                            AX_lists_neg.append(0)

                        cols_list_sig.append(colors_dict[SB_cat]["sig"])
                        cols_list_unsig.append(colors_dict[SB_cat]["nonsig"])

            ## make tick positions so that the groups of three within each chromosome are closer together
            tick_pos = [i for i in range(1,len(tick_labels)+1)]
            box_width = 1/5
            pos_adjust=box_width*1.3
            for i in range(len(tick_labels)):
                if i%3 == 0:
                    tick_pos[i] =1+ i//3 -pos_adjust
                if i%3 == 1:
                    tick_pos[i] =1+ i//3 
                if i%3 == 2:
                    tick_pos[i] =1+ i//3 +pos_adjust
                print(f"{i} : {tick_pos[i]}")

            ax.bar(tick_pos, AX_lists_neg, width = box_width, label='M1a', color=cols_list_unsig)#, alpha=0.7)
            ax.bar(tick_pos, AX_lists_pos, width = box_width, bottom=AX_lists_neg, label='M2a', color= cols_list_sig)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<1 else f'{int(x)}%'))

        else:

            ## make tick labels in the same order as the box plot lists to be sure that everyting is right
            tick_labels = []
            AX_lists = []
            print(f"{title}")
            for i in range(1,6): # hard-code rank order. dicts are ordered so it should be right by default but just to be sure
                chr_dict = AX_dicts[i]
                print(f"rank {i}")
                for chr in ["A","X"]:
                    print(f"\t - {chr}")
                    SB_dict = chr_dict[chr]
                    for SB_cat in SB_order_list:

                        ### to double check plot coloring and stuff right: plot verbose tick labels!
                        # num_genes =f"({len(SB_dict[SB_cat])}):{chr}:{SB_cat}"
                        ### 

                        num_genes =f"({len(SB_dict[SB_cat])})"
                        print(f"\t\t - {SB_cat} has {num_genes} genes")
                        tick_labels.append(num_genes)
                        AX_lists.append(SB_dict[SB_cat])
            ## make tick positions so that the groups of three within each chromosome are closer together
            tick_pos = [i for i in range(1,len(tick_labels)+1)]
            box_width = 1/5
            pos_adjust=box_width*1.3
            for i in range(len(tick_labels)):
                if i%3 == 0:
                    tick_pos[i] =1+ i//3 -pos_adjust
                if i%3 == 1:
                    tick_pos[i] =1+ i//3 
                if i%3 == 2:
                    tick_pos[i] =1+ i//3 +pos_adjust
                print(f"{i} : {tick_pos[i]}")

            bp = ax.boxplot(AX_lists, positions=tick_pos, widths=box_width, patch_artist=True)
            
            ## make fancy colors for the boxplot
            if True:
                for i, box in enumerate(bp["boxes"]):
                    box.set(facecolor=color_faces[i % 3], edgecolor=color_edges[i % 3], linewidth = lw)
                for i, median in enumerate(bp['medians']):
                    median.set(color=color_medians[i % 3], linewidth=lw)
                for i, whisker in enumerate(bp['whiskers']):
                    whisker.set(color = color_edges[i//2 % 3], linewidth = lw, linestyle='-')
                for i, cap in enumerate(bp['caps']):
                    cap.set(color = color_edges[i//2 % 3], linewidth = lw)
                for i, flier in enumerate(bp['fliers']):
                    flier.set(markerfacecolor=color_edges[i % 3], markeredgecolor=color_edges[i % 3], linewidth = lw, marker='.')



        # set axis labels
        tick_fs_factor = 0.7
        ax.set_xticks(ticks = tick_pos, labels = tick_labels, fontsize=fs*tick_fs_factor)#, rotation=45, ha='right')
        ax.tick_params(axis='x', labelsize=fs*tick_fs_factor, rotation = 90)
        ax.tick_params(axis='y', labelsize=fs*0.9)
        ax.set_title(title, fontsize=fs)
        ax.set_xlim(min(tick_pos)-pos_adjust, max(tick_pos)+pos_adjust)

        ## plot rank labels
        fs_factor=1.5
        ax2 = ax.secondary_xaxis('bottom')
        ax2.set_xticks([i+0.5 for i in range(1,10,2)])
        ax2.set_xticklabels([i for i in range(1,6)], fontsize=fs*fs_factor)
        ax2.spines['bottom'].set_position(('outward', 120))   
        ax2.xaxis.set_ticks_position('none')
        ax2.spines['bottom'].set_visible(False)
        ax2.tick_params(axis='x', labelsize=fs*fs_factor)

        ## plot chromosome categories
        ax3 = ax.secondary_xaxis('bottom')
        ax3.set_xticks([i for i in range(1,11)])
        ax3.set_xticklabels(["A" if i%2==1 else "X" for i in range(1,11)], fontsize=fs)
        ax3.spines['bottom'].set_position(('outward', 80))   
        ax3.xaxis.set_ticks_position('none')
        ax3.spines['bottom'].set_visible(False)
        ax3.tick_params(axis='x', labelsize=fs)

        ## highlight X-chromosomes
        for i in range(1,11):
            if i % 2 == 0:  # only the X-sections
                left = i-box_width-pos_adjust
                right = i+box_width+pos_adjust
                ax.axvspan(left, right, color="#7A6B70", alpha=0.2, linewidth = 0, zorder=0)

        fig.supxlabel(f"conservation rank, chromosome location, and (number of genes)", fontsize = fs)
        fig.supylabel(y_label, fontsize = fs, x=0.0, y=0.625)

        # layout (left, bottom, right, top)
        plt.tight_layout(rect=[0.0, 0.05, 1, 1])

        # transparent background
        plt.savefig(outfile, dpi = 300, transparent = True)
        # non-transparent background
        filename_tr = outfile.replace(".png", "_white_bg.png")
        plt.savefig(filename_tr, dpi = 300, transparent = False)
        print(f"plot saved in current working directory as: {outfile} and {filename_tr}")

    colors = {
        "male" : {
            "fill" : "#495E83", # dusk blue
            "edge" : "#374C6E", # dusk blue darker
            "medians" : "#A7CCED", # icy blue
            "sig" : "#3C5786", # dusk blue
            "nonsig" : "#677CA2", # glaucous
        },
        "female" : {
            "fill" : "#AB354A", # cherry rose
            "edge" : "#771C2C", # dark amaranth
            "medians" : "#EA9AA9", # cotton candy
            "sig" : "#A4273E", # cherry rose
            "nonsig" : "#BC5265" # dusty mauve
        },
        "unbiased" : {
            "fill" : "#816874", # dusty mauve
            "edge" : "#5B414E", # mauve shadow
            "medians" : "#C3A7B5", # lilac ash
            "sig" : "#6B555F", # taupe grey
            "nonsig" : "#967A88" #dusty mauve
        }
    }

    if pos_sel:
        title_suffix = f"pos. sel. genes in"
    else:
        title_suffix = f"dN/dS by sex bias for"

    plot_dNdS_subplot(AX_dicts=nested_vals_dict_a, outfile=outfile.replace(".png", f"_abdomen.png"), 
        title=f"{title_suffix} abdominal tissue", colors_dict=colors, fs=fs, pos_sel=pos_sel)
    plot_dNdS_subplot(AX_dicts=nested_vals_dict_ht,outfile=outfile.replace(".png", f"_head_thorax.png"), 
        title=f"{title_suffix} head+thorax tissue", colors_dict=colors, fs=fs, pos_sel=pos_sel)



def make_phylogeny_rank_merged_dict(summary_file_df, tissue, max_dNdS=2, pos_sel =False):
    """
    make a dict with X/Y like this where the rank orders are merged 
    {
        "A" : { "male" : [list, of, dNdS, numbers] ,
                "female" : [more, dNdS, numbers] , 
                "unbiased" : [...] },
        "X" : { "male" : [list, of, dNdS, numbers] ,
                "female" : [more, dNdS, numbers] , 
                "unbiased" : [...] },
    }
    """

    summary_file_df["dN/dS"] = pd.to_numeric(summary_file_df["dN/dS"], errors='coerce')
    filtered_df = summary_file_df.dropna(subset="dN/dS")
    if max_dNdS>0:
        filtered_df = filtered_df[filtered_df["dN/dS"] < max_dNdS]
    
    dNdS_dict = { "A": {"male": [], "female": [], "unbiased": []},
                  "X": {"male": [], "female": [], "unbiased": []} }

    if pos_sel:
        for SB_cat, possel, chr  in zip(filtered_df[f"SB_{tissue}"],filtered_df["positive_selection"],filtered_df["chromosome"]):
            dNdS_dict[chr][SB_cat].append(possel)
    else:
        for SB_cat, dNdS, chr  in zip(filtered_df[f"SB_{tissue}"],filtered_df["dN/dS"],filtered_df["chromosome"]):
            dNdS_dict[chr][SB_cat].append(dNdS)

    return dNdS_dict


def boxplot_dNdS_merge_rank(full_table_paths_dict, outfile, maxdNdS=2, partner_species="C_chinensis", pos_sel=False):

    X_df = pd.read_csv(full_table_paths_dict["X"], sep="\t")
    A_df = pd.read_csv(full_table_paths_dict["A"], sep="\t")
    X_df["chromosome"] = ["X"]*X_df.shape[0]
    A_df["chromosome"] = ["A"]*A_df.shape[0]
    df = pd.concat([A_df,X_df], ignore_index=True)

    df = df.rename(columns={'LFC_head+thorax': 'LFC_head_thorax'})
    df = df.rename(columns={'FDR_pval_head+thorax': 'FDR_pval_head_thorax'})
    print(df)
    df["SB_abdomen"] = df.apply(make_sex_bias_cat_row, axis=1, args=("abdomen",))
    df["SB_head_thorax"] = df.apply(make_sex_bias_cat_row, axis=1, args=("head_thorax",))
    df = df[df["other_species"]==partner_species]

    if pos_sel:
        nested_vals_dict_a = make_phylogeny_rank_merged_dict(df, tissue="abdomen", max_dNdS=maxdNdS, pos_sel=pos_sel)
        nested_vals_dict_ht = make_phylogeny_rank_merged_dict(df, tissue="head_thorax", max_dNdS=maxdNdS, pos_sel=pos_sel)
    else:
        nested_vals_dict_a = make_phylogeny_rank_merged_dict(df, tissue="abdomen", max_dNdS=maxdNdS)
        nested_vals_dict_ht = make_phylogeny_rank_merged_dict(df, tissue="head_thorax", max_dNdS=maxdNdS)

    
    y_label = f"dN/dS"
    if maxdNdS>0:
        y_label = f"dN/dS (max. {maxdNdS})"
    
    fs = 30 # font size

    # set figure aspect ratio
    if pos_sel:
        aspect_ratio = 4 / 8
    else:
        aspect_ratio = 3 / 8
    height_pixels = 1000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    def plot_dNdS_subplot(AX_dicts, fs, title, colors_dict, lw=2, outfile = "plot.png"):

        fig, ax = plt.subplots(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
        SB_order_list =["male","unbiased","female"]
        
        color_faces = []
        color_edges = []
        color_medians = []

        for SB_cat in SB_order_list:
            color_faces.append(colors_dict[SB_cat]["fill"])
            color_edges.append(colors_dict[SB_cat]["edge"])
            color_medians.append(colors_dict[SB_cat]["medians"])
        

        if pos_sel:
            ## make tick labels in the same order as the box plot lists to be sure that everyting is right
            tick_labels = []
            AX_lists_pos = []
            AX_lists_neg = []
            cols_list_sig = []
            cols_list_unsig = []

            print(f"{title}")
            for chr in ["A","X"]:
                print(f"\t - {chr}")
                SB_dict = AX_dicts[chr]
                for SB_cat in SB_order_list:

                    ### to double check plot coloring and stuff right: plot verbose tick labels!
                    # num_genes =f"({len(SB_dict[SB_cat])}):{chr}:{SB_cat}"
                    ### 

                    num_genes =len(SB_dict[SB_cat])
                    print(f"\t\t - {SB_cat} has {num_genes} genes")
                    tick_labels.append(f"({num_genes})")

                    A_nonsig = len([val for val in SB_dict[SB_cat] if val == False])*100.0
                    A_sig = len([val for val in SB_dict[SB_cat] if val == True])*100.0
                    try:
                        AX_lists_pos.append(A_sig/num_genes)
                    except:
                        AX_lists_pos.append(0)
                    try:
                        AX_lists_neg.append(A_nonsig/num_genes)
                    except:
                        AX_lists_neg.append(0)

                    cols_list_sig.append(colors_dict[SB_cat]["sig"])
                    cols_list_unsig.append(colors_dict[SB_cat]["nonsig"])

            ## make tick positions so that the groups of three within each chromosome are closer together
            tick_pos = [i for i in range(1,len(tick_labels)+1)]
            box_width = 1/5
            pos_adjust=box_width*1.3
            for i in range(len(tick_labels)):
                if i%3 == 0:
                    tick_pos[i] =1+ i//3 -pos_adjust
                if i%3 == 1:
                    tick_pos[i] =1+ i//3 
                if i%3 == 2:
                    tick_pos[i] =1+ i//3 +pos_adjust
                print(f"{i} : {tick_pos[i]}")

            ax.bar(tick_pos, AX_lists_neg, width = box_width, label='M1a', color=cols_list_unsig)#, alpha=0.7)
            ax.bar(tick_pos, AX_lists_pos, width = box_width, bottom=AX_lists_neg, label='M2a', color= cols_list_sig)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<1 else f'{int(x)}%'))



        else:
            ## make tick labels in the same order as the box plot lists to be sure that everyting is right
            tick_labels = []
            AX_lists = []
            print(f"{title}")
            for chr in ["A","X"]:
                print(f"\t - {chr}")
                SB_dict = AX_dicts[chr]
                for SB_cat in SB_order_list:

                    ### to double check plot coloring and stuff right: plot verbose tick labels!
                    # num_genes =f"({len(SB_dict[SB_cat])}):{chr}:{SB_cat}"
                    num_genes =f"({len(SB_dict[SB_cat])})"
                    ### 

                    print(f"\t\t - {SB_cat} has {num_genes} genes")
                    tick_labels.append(num_genes)
                    AX_lists.append(SB_dict[SB_cat])

            ## make tick positions so that the groups of three within each chromosome are closer together
            tick_pos = [i for i in range(1,len(tick_labels)+1)]
            box_width = 1/5
            pos_adjust=box_width*1.3
            for i in range(len(tick_labels)):
                if i%3 == 0:
                    tick_pos[i] =1+ i//3 -pos_adjust
                if i%3 == 1:
                    tick_pos[i] =1+ i//3 
                if i%3 == 2:
                    tick_pos[i] =1+ i//3 +pos_adjust
                print(f"{i} : {tick_pos[i]}")
            
            bp = ax.boxplot(AX_lists, positions=tick_pos, widths=box_width, patch_artist=True)



        # set axis labels
        tick_fs_factor = 0.7
        ax.set_xticks(ticks = tick_pos, labels = tick_labels, fontsize=fs*tick_fs_factor)#, rotation=45, ha='right')
        ax.tick_params(axis='x', labelsize=fs*tick_fs_factor, rotation = 90)
        ax.tick_params(axis='y', labelsize=fs*0.9)
        ax.set_title(title, fontsize=fs*0.8)
        ax.set_xlim(min(tick_pos)-pos_adjust, max(tick_pos)+pos_adjust)

        ## plot chromosome categories
        ax3 = ax.secondary_xaxis('bottom')
        ax3.set_xticks([i+1 for i in range(2)])
        ax3.set_xticklabels(["A" if i%2==1 else "X" for i in range(1,3)], fontsize=fs)
        ax3.spines['bottom'].set_position(('outward', 80))   
        ax3.xaxis.set_ticks_position('none')
        ax3.spines['bottom'].set_visible(False)
        ax3.tick_params(axis='x', labelsize=fs)

        ## highlight X-chromosomes
        for i in range(1,11):
            if i % 2 == 0:  # only the X-sections
                left = i-box_width-pos_adjust
                right = i+box_width+pos_adjust
                ax.axvspan(left, right, color="#7A6B70", alpha=0.2, linewidth = 0, zorder=0)

        if pos_sel==False:
            ## make fancy colors for the boxplot
            if True:
                for i, box in enumerate(bp["boxes"]):
                    box.set(facecolor=color_faces[i % 3], edgecolor=color_edges[i % 3], linewidth = lw)
                for i, median in enumerate(bp['medians']):
                    median.set(color=color_medians[i % 3], linewidth=lw)
                for i, whisker in enumerate(bp['whiskers']):
                    whisker.set(color = color_edges[i//2 % 3], linewidth = lw, linestyle='-')
                for i, cap in enumerate(bp['caps']):
                    cap.set(color = color_edges[i//2 % 3], linewidth = lw)
                for i, flier in enumerate(bp['fliers']):
                    flier.set(markerfacecolor=color_edges[i % 3], markeredgecolor=color_edges[i % 3], linewidth = lw, marker='.')

        fig.supxlabel(f"(number of genes)\nchromosome", fontsize = fs*0.9)
        if pos_sel:
            fig.supylabel(y_label, fontsize = fs*0.8, x=0.1, y=0.625)
        else:
            fig.supylabel(y_label, fontsize = fs*0.8, x=0.0, y=0.625)
        # layout (left, bottom, right, top)
        plt.tight_layout(rect=[0.0, 0.05, 1, 1])

        # transparent background
        plt.savefig(outfile, dpi = 300, transparent = True)
        # non-transparent background
        filename_tr = outfile.replace(".png", "_white_bg.png")
        plt.savefig(filename_tr, dpi = 300, transparent = False)
        print(f"plot saved in current working directory as: {outfile} and {filename_tr}")

    colors = {
        "male" : {
            "fill" : "#495E83", # dusk blue
            "edge" : "#374C6E", # dusk blue darker
            "medians" : "#A7CCED", # icy blue
            "sig" : "#3C5786", # dusk blue
            "nonsig" : "#677CA2", # glaucous
        },
        "female" : {
            "fill" : "#AB354A", # cherry rose
            "edge" : "#771C2C", # dark amaranth
            "medians" : "#EA9AA9", # cotton candy
            "sig" : "#A4273E", # cherry rose
            "nonsig" : "#BC5265" # dusty mauve
        },
        "unbiased" : {
            "fill" : "#816874", # dusty mauve
            "edge" : "#5B414E", # mauve shadow
            "medians" : "#C3A7B5", # lilac ash
            "sig" : "#6B555F", # taupe grey
            "nonsig" : "#967A88" #dusty mauve
        }
    }


    plot_dNdS_subplot(AX_dicts=nested_vals_dict_a, outfile=outfile.replace(".png", f"_abdomen.png"), 
        title="dN/dS abd.", colors_dict=colors, fs=fs)
    plot_dNdS_subplot(AX_dicts=nested_vals_dict_ht,outfile=outfile.replace(".png", f"_head_thorax.png"), 
        title="dN/dS h+t", colors_dict=colors, fs=fs)





if __name__ == "__main__":

    username = "miltr339"
    full_tables_dict = get_full_table_path(username=username)
    reorg_table_outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/paml_stats_outfile_table.tsv"
    
    if False:
        ###################################################
        ## median quantile regression for dNdS as continuous response
        statistical_analysis_dNdS(full_tables_dict, table_outfile=f"")
        # compare_conservation_rank_proportions(full_tables_dict)
        ###################################################

    if True:
        pos_sel = True # if true plot bar charts with proportion of positive selection
        if False:
            ###################################################
            ### boxplot for positive selection or dNdS by rank/chromosome/sex bias
                           # if false plot boxplot of dNdS values
            if pos_sel:
                filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/pos_sel_vs_conservation_rank_boxplot.png"
            else:
                filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_vs_conservation_rank_boxplot.png"
            boxplot_dNdS(full_tables_dict, outfile=filename, pos_sel=pos_sel)
            ###################################################
        else:
            ###################################################
            ### boxplot for dNdS by chromosome/sex bias with rank categories merged
            if pos_sel:
                filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/pos_sel_merged_conservation_rank_boxplot.png"
            else:
                filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_merged_conservation_rank_boxplot.png"
            boxplot_dNdS_merge_rank(full_tables_dict, outfile=filename, pos_sel=pos_sel)
            ###################################################


    if False:
        ###################################################
        ## logistic regression for categorical response (positive selection True/False)
        statistical_analysis_pos_sel(full_table_paths_dict=full_tables_dict)
        ###################################################

    if False:
        ###################################################
        ## scatterplots dNdS vs logFC 
        filename=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_vs_conservation_rank.png"
        plot_dNdS_rank_conserved(summary_paths_AX_list= full_tables_dict, outfile = filename, maxdNdS = 0)