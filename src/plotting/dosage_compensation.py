"""
Assess dosage compensation by comparing mean expression on A and X-linked genes in males and females
"""

import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt
import parse_gff as gff
import scipy.stats as sts
from statannotations.Annotator import Annotator


username="milena"

rpkm_norm_counts = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_gene_counts_edgeR_rpkm_length_normalized.txt"
Cmac_annotation = f"/Users/{username}/work/native_annotations/all_native_annot/C_maculatus_superscaffolded_LomeRNA_braker_isoform_filtered.gff"
sex_chromosomes = { 
    "X" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'],
    "Y" : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164'],
    "all" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83','scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']}


def get_geneID_list(gff_contig_dict, verbose=True):
    """get geneIDs list of all geneIDs included in the output of gff.parse_gff3_by_contig()"""
    out_list = []
    for scaffold, feature_list in gff_contig_dict.items():
        if verbose:
            print(f"\t{scaffold}: {len(feature_list)} genes")
        for feature in feature_list:
            out_list.append(feature.feature_id)
    print(f"{len(out_list)} features in total\n")
    return out_list


def make_boxplot(data_dict:dict, outfile:str, tissue:str):
    fs = 15 # font size
    lw = 2 # line width

    colors_dict = {
        "fill" : "#495E83", # dusk blue
        "edge" : "#374C6E", # dusk blue darker
        "medians" : "#A7CCED", # icy blue
        "lines" : "#7D93B5", # lavender grey
        "FB_fill" : "#AB354A", # cherry rose
        "FB_edge" : "#771C2C", # dark amaranth
        "FB_medians" : "#EA9AA9", # cotton candy
        # "FB_lines" : "#7D93B5" # lavender grey
    }

    # set figure aspect ratio
    aspect_ratio = 18 / 12 # height / width
    height_pixels = 1400  # Height in pixels
    dpi = 300
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(1,1,figsize=(width_pixels/dpi, height_pixels/dpi))
    
    tick_labels = [f"{key}\n({len(lists)})" for key,lists in data_dict.items()]
    lists = [means_list for means_list in data_dict.values()]
    print(len(lists))

    bp = ax.boxplot(lists, patch_artist=True)   
    # set axis labels
    ax.set_ylabel("fpkm normalized counts", fontsize=fs)
    ax.set_yscale('log')
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=fs) 
    ax.set_xticks(ticks = range(1,len(lists)+1), labels = tick_labels, fontsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    title = f"Expression of X and A linked genes\nin {tissue} tissues"
    ax.set_title(title, fontsize=fs)

    ## modify boxplot colors
    if True:
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

    # layout rect=(left, bottom, right, top)
    plt.tight_layout()# rect=[0.0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = dpi, transparent = True)
    print(f"plot saved in current working directory as: {outfile}")


def test_median_diff(data_dict):
    """mood's median test as implemented in scipy.stats.median_Test"""
    for sex in ["male","female"]:
        print(f" ......... {sex} .........")
        statistic,pvalue,median,table = sts.median_test(data_dict[f"{sex} X"], data_dict[f"{sex} A"])
        print(f"\tMood's : \t Chi-sq statistic: {statistic}, p-value: {pvalue}")

        stat, pval = sts.mannwhitneyu(data_dict[f"{sex} X"], data_dict[f"{sex} A"], alternative='two-sided')
        print(f"\tMann-Whitney : \t U statistic: {stat}, p-value: {pval}")
        
    print()


def make_boxplot_both_tissues(data_dict:dict, outfile:str, tissue:str):
    fs = 18 # font size
    lw = 2 # line width

    colors_dict = {
        "fill" : "#495E83", # dusk blue
        "edge" : "#374C6E", # dusk blue darker
        "medians" : "#A7CCED", # icy blue
        "lines" : "#7D93B5", # lavender grey
        "FB_fill" : "#AB354A", # cherry rose
        "FB_edge" : "#771C2C", # dark amaranth
        "FB_medians" : "#EA9AA9", # cotton candy
        # "FB_lines" : "#7D93B5" # lavender grey
    }

    # set figure aspect ratio
    aspect_ratio = 32 / 12 # width / height
    height_pixels = 1400  # Height in pixels
    dpi = 300
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(1,1,figsize=(width_pixels/dpi, height_pixels/dpi))

    title = f"Expression of X and A linked genes"
    # samples_sorted_keys = sorted(data_dict.keys())
    ## hard code the order of the boxplots:
    samples_sorted_keys = [
        'reproductive_male_X', 
        'reproductive_male_A', 
        'reproductive_female_X', 
        'reproductive_female_A', 
        'somatic_male_X',
        'somatic_male_A', 
        'somatic_female_X', 
        'somatic_female_A', 
    ]

    tick_labels = []
    for key in samples_sorted_keys:
        means_list = data_dict[key]
        lablist = key.split("_")
        tissue = lablist[0]
        sex = lablist[1]
        chr = lablist[2]
        label_name = chr
        # tick_labels.append(f"{tissue}\n{label_name}\n({len(means_list)})")
        tick_labels.append(f"({len(means_list)})\n{label_name}")
    lists = [data_dict[key] for key in samples_sorted_keys]
    offset=0.1
    shrink=1
    tickpos = [i*shrink+offset if i%2==1 else i*shrink-offset for i in range(1,len(lists)+1)]
    print(len(lists))

    bp = ax.boxplot(lists, positions = tickpos, widths=0.6, patch_artist=True)   

    # set axis labels
    ax.set_ylabel("fpkm normalized counts", fontsize=fs)
    ax.set_yscale('log')
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=fs) 
    ax.set_xticks(ticks = tickpos, labels = tick_labels, fontsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title(title, fontsize=fs)
    # label reproductive and somatics
    half = len(tickpos)//2
    ax2 = ax.secondary_xaxis('bottom')
    ax2.set_xticks([mean(tickpos[:half]), mean(tickpos[half:])])
    ax2.set_xticklabels(["abdomen","head+thorax"], fontsize=fs)
    ax2.spines['bottom'].set_position(('outward', 40))   
    ax2.xaxis.set_ticks_position('none')
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(axis='x', labelsize=fs)

    ## modify boxplot colors
    ## (i//1) and (i//2) for alternating colors
    ## (i//2) and (i//4) for pairs of two
    if True:
        for i, box in enumerate(bp['boxes']):
            if (i//2) % 2==0:
                box.set(facecolor=colors_dict["fill"], edgecolor=colors_dict["edge"], linewidth=2)
            else:
                box.set(facecolor=colors_dict["FB_fill"], edgecolor=colors_dict["FB_edge"], linewidth=2)
        for i, median in enumerate(bp['medians']):
            if (i//2) % 2==0:
                median.set(color=colors_dict['medians'], linewidth=lw)
            else:
                median.set(color=colors_dict['FB_medians'], linewidth=lw)
        for i, whisker in enumerate(bp['whiskers']):
            # print(f"whisker: {i}")
            if (i//4) % 2==0:
                whisker.set(color=colors_dict['edge'], linestyle='-',linewidth=lw)
            else:
                whisker.set(color=colors_dict['FB_edge'], linestyle='-',linewidth=lw)
        for i, cap in enumerate(bp['caps']):
            if (i//4) % 2==0:
                cap.set(color=colors_dict['edge'],linewidth=lw)
            else:
                cap.set(color=colors_dict['FB_edge'],linewidth=lw)
        for i, flier in enumerate(bp['fliers']):
            if (i//2) % 2==0:
                flier.set(marker='.', markerfacecolor=colors_dict['edge'], markeredgecolor=colors_dict['edge'])
            else:
                flier.set(marker='.', markerfacecolor=colors_dict['FB_edge'], markeredgecolor=colors_dict['FB_edge'])


    ## make statistical annotation
    if True:
        def add_significance_bar_log(ax, x1, x2, sample1, sample2, data, y, factor=1.5, color='black', lw=lw, fs=fs):
            """
            This function was modified from one created by claude code
            """
            # run the test
            data1 = data[sample1]
            data2 = data[sample2]
            stat, p = sts.mannwhitneyu(data1, data2, alternative='two-sided')

            # convert p-value to stars
            text = 'ns'
            color="#8e8e8e" #light grey
            tick_top = y * factor
            if p < 0.05:
                text = '*'
                color="#343434"# darker grey
                ax.text((x1+x2)/2, tick_top*0.7, text, ha='center', va='bottom', color=color, fontsize=fs)
                print(f"  * Mann-Whitney {sample1}-{sample2} : \t U statistic: {stat}, p-value: {p}")
            else:
                print(f"    Mann-Whitney {sample1}-{sample2} : \t U statistic: {stat}, p-value: {p}")
                ax.text((x1+x2)/2, tick_top, text, ha='center', va='bottom', color=color, fontsize=fs)
            ax.plot([x1, x1, x2, x2], [y, tick_top, tick_top, y], lw=lw, color=color)
                
        y0 = max(data_dict['reproductive_male_X'] + data_dict['reproductive_male_A'] + data_dict['reproductive_female_X'] + data_dict['reproductive_female_A'] + 
                data_dict['somatic_male_X'] + data_dict['somatic_male_A'] + data_dict['somatic_female_X'] + data_dict['somatic_female_A']) * 1.1

        add_significance_bar_log(ax=ax, x1=tickpos[0], x2=tickpos[1], sample1=samples_sorted_keys[0], sample2=samples_sorted_keys[1], data=data_dict, y=y0 * 2)
        add_significance_bar_log(ax=ax, x1=tickpos[2], x2=tickpos[3], sample1=samples_sorted_keys[2], sample2=samples_sorted_keys[3], data=data_dict, y=y0 * 2)
        add_significance_bar_log(ax=ax, x1=tickpos[4], x2=tickpos[5], sample1=samples_sorted_keys[4], sample2=samples_sorted_keys[5], data=data_dict, y=y0 * 2)
        add_significance_bar_log(ax=ax, x1=tickpos[6], x2=tickpos[7], sample1=samples_sorted_keys[6], sample2=samples_sorted_keys[7], data=data_dict, y=y0 * 2)
        
        add_significance_bar_log(ax=ax, x1=tickpos[0], x2=tickpos[2], sample1=samples_sorted_keys[0], sample2=samples_sorted_keys[2], data=data_dict, y=y0 * 50)
        add_significance_bar_log(ax=ax, x1=tickpos[4], x2=tickpos[6], sample1=samples_sorted_keys[4], sample2=samples_sorted_keys[6], data=data_dict, y=y0 * 50)

        ax.set_ylim(0,y0*1000)
        


    # layout rect=(left, bottom, right, top)
    plt.tight_layout()# rect=[0.0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = dpi, transparent = True)
    print(f"plot saved in current working directory as: {outfile}")




if __name__ == "__main__":

    # get A and X gene lists
    gff_contig = gff.parse_gff3_by_contig(Cmac_annotation, featurecategory=[gff.FeatureCategory.Gene])
    gff_X = { scaffold : genes_list for scaffold,genes_list in gff_contig.items() if scaffold in sex_chromosomes["X"]}
    gff_A = { scaffold : genes_list for scaffold,genes_list in gff_contig.items() if scaffold not in sex_chromosomes["all"]}
    X_genes = get_geneID_list(gff_X)
    A_genes = get_geneID_list(gff_A, verbose=False)

    # read normalized counts
    column_categories = {
        "reproductive_F" : ["AVf1","AVf2","AVf3"],
        "reproductive_M" : ["AVm1","AVm2","AVm3"],
        "somatic_F" : ["HtVf1","HtVf2","HtVf3"],
        "somatic_M" : ["HtVm1","HtVm2","HtVm3"]
    }

    tissues=["somatic","reproductive"]

    if False:
        for tissue in tissues:
            print(f"tissue : {tissue}")
            df = pd.read_csv(rpkm_norm_counts, sep="\t")
            X_df = df.loc[X_genes]
            A_df = df.loc[A_genes]
            ## make a column of gene-wise means
            X_F_df = X_df[column_categories[f"{tissue}_F"]]
            X_F_df["gene_mean"] = X_F_df.mean(axis=1)
            X_M_df = X_df[column_categories[f"{tissue}_M"]]
            X_M_df["gene_mean"] = X_M_df.mean(axis=1)
            A_F_df = A_df[column_categories[f"{tissue}_F"]]
            A_F_df["gene_mean"] = A_F_df.mean(axis=1)
            A_M_df = A_df[column_categories[f"{tissue}_M"]]
            A_M_df["gene_mean"] = A_M_df.mean(axis=1)
            
            filt_expr =True
            if not filt_expr:
                # no expression filtering to remove 0-counts genes
                sample_mean_expressions = {
                    "male X" : X_M_df["gene_mean"].tolist(),
                    "female X" : X_F_df["gene_mean"].tolist(),
                    "male A" : A_M_df["gene_mean"].tolist(),
                    "female A" : A_F_df["gene_mean"].tolist(),
                }
                outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DC_boxplot_{tissue}.png"    
            else:
                # remove genes that have mean expression 0 across all samples from the list
                min_expr = 0
                sample_mean_expressions = {
                    "male X" : [expr_val for expr_val in X_M_df["gene_mean"].tolist() if expr_val>min_expr],
                    "female X" : [expr_val for expr_val in X_F_df["gene_mean"].tolist() if expr_val>min_expr],
                    "male A" : [expr_val for expr_val in A_M_df["gene_mean"].tolist() if expr_val>min_expr],
                    "female A" : [expr_val for expr_val in A_F_df["gene_mean"].tolist() if expr_val>min_expr],
                }
                outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DC_boxplot_{tissue}_filtered.png"    


            for key,val in sample_mean_expressions.items():
                print(f"{key} : {len(val)}")

            make_boxplot(data_dict=sample_mean_expressions, outfile = outfile, tissue=tissue)

            ## test for statistical differences: 
            print(f"\n -------------- {tissue} stats (filtered: {filt_expr}) --------------")
            test_median_diff(data_dict=sample_mean_expressions)

    ## make shared plot of somatic and reproductive for publication
    if True:
        min_expr=0

        df = pd.read_csv(rpkm_norm_counts, sep="\t")
        X_df = df.loc[X_genes]
        A_df = df.loc[A_genes]

        sample_mean_expressions = {}
        tissues=["somatic","reproductive"]

        for tissue in tissues:
            X_F_df = X_df[column_categories[f"{tissue}_F"]]
            X_F_df["gene_mean"] = X_F_df.mean(axis=1)
            X_M_df = X_df[column_categories[f"{tissue}_M"]]
            X_M_df["gene_mean"] = X_M_df.mean(axis=1)
            A_F_df = A_df[column_categories[f"{tissue}_F"]]
            A_F_df["gene_mean"] = A_F_df.mean(axis=1)
            A_M_df = A_df[column_categories[f"{tissue}_M"]]
            A_M_df["gene_mean"] = A_M_df.mean(axis=1)
            sample_mean_expressions[f"{tissue}_male_X"] = [expr_val for expr_val in X_M_df["gene_mean"].tolist() if expr_val>min_expr]
            sample_mean_expressions[f"{tissue}_female_X"] = [expr_val for expr_val in X_F_df["gene_mean"].tolist() if expr_val>min_expr]
            sample_mean_expressions[f"{tissue}_male_A"] = [expr_val for expr_val in A_M_df["gene_mean"].tolist() if expr_val>min_expr]
            sample_mean_expressions[f"{tissue}_female_A"] = [expr_val for expr_val in A_F_df["gene_mean"].tolist() if expr_val>min_expr]
            
        outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DC_boxplot_both_tissues_filtered.png"    
        make_boxplot_both_tissues(data_dict=sample_mean_expressions, outfile = outfile, tissue=tissue)