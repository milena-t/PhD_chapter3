"""
Assess dosage compensation by comparing mean expression on A and X-linked genes in males and females
"""

import pandas as pd
import matplotlib.pyplot as plt
import parse_gff as gff

username="miltr339"
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
    ax.set_ylabel("rpkm normalized counts", fontsize=fs)
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
        sample_mean_expressions = {
            "male X" : X_M_df["gene_mean"].tolist(),
            "female X" : X_F_df["gene_mean"].tolist(),
            "male A" : A_M_df["gene_mean"].tolist(),
            "female A" : A_F_df["gene_mean"].tolist(),
        }
        for key,val in sample_mean_expressions.items():
            print(f"{key} : {len(val)}")

        make_boxplot(data_dict=sample_mean_expressions, outfile = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/DC_boxplot_{tissue}.png", tissue=tissue)
    

