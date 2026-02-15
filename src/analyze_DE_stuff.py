"""
analyze the info from the big summary tables generated with PhD_chapter3/src/make_DE_analysis_table.py
"""

import parse_gff as gff
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
import scipy.stats as sts


def get_summary_paths(username = "miltr339"):
    summary_dict = {
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_A_chr.tsv",
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_X_chr.tsv",
    }
    return summary_dict

def get_DE_paths(username = "miltr339"):
    paths = {
        "head_thorax" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_sex_DE_head_thorax_edgeR.txt",
        "abdomen" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_sex_DE_abdomen_edgeR.txt",
    }
    return paths

def get_Cmac_superscaffolded_XY_contigs():
    return { "X" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'], "Y" : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']}


def get_contig_lengths(contigs_list, assembly_index):
    """
    get the length of all contigs in question from the assembly
    returns a dict with { contigID : length }
    """
    out_dict = {contig : 0 for contig in contigs_list}  
    with open(assembly_index) as f:
        for line in f:
            fields = line.split("\t")
            if fields[0] in contigs_list:
                out_dict[fields[0]] = int(fields[1])
    return out_dict





def plot_dosage_compensation(summary_paths, annotation, assembly_index, X_list, outfile = ""):
    """
    plot the log2FC of all genes on the X in both head+thorax and abdominal tissues on the y,
    the X-axis is the position of the gene on the X
    """

    summary_data_abdomen = pd.read_csv(summary_paths["abdomen"], sep = "\t", index_col=False)
    summary_data_head_thorax = pd.read_csv(summary_paths["head_thorax"], sep = "\t", index_col=False)
    filtered_abdomen = summary_data_abdomen[pd.notna(summary_data_abdomen["logFC"])]
    filtered_head_thorax = summary_data_head_thorax[pd.notna(summary_data_head_thorax["logFC"])]

    ## get LFC information for every gene for plotting { transcriptID : LFC_float }
    LFC_dict_abdomen = dict(zip(filtered_abdomen["geneID"], filtered_abdomen["logFC"]))
    LFC_dict_head_thorax = dict(zip(filtered_head_thorax["geneID"], filtered_head_thorax["logFC"]))

    assert sorted(list(LFC_dict_abdomen.keys())) == sorted(list(LFC_dict_head_thorax.keys())) # check that there is data for all transcripts in both

    ## parse gene position data by contig for plotting
    annotation_dict = gff.parse_gff3_general(annotation)
    coordinates_dict = { contig : {} for contig in X_list} # dict with { contig : { transcriptID : [start,end], ... }  }
    for transcriptID in LFC_dict_abdomen.keys():
        transcript = annotation_dict[transcriptID]
        if transcript.contig not in X_list:
            # print(transcript)
            continue
        # coordinates_dict[transcript.contig][transcript.feature_id] = [transcript.start, transcript.end]
        coordinates_dict[transcript.contig][transcript.feature_id] = np.mean([transcript.start, transcript.end])

    ## get the X contig lenghts for plotting the X axis correctly
    X_contig_lengths_dict = get_contig_lengths(contigs_list=X_list,assembly_index=assembly_index)
    total_x_length = sum(list(X_contig_lengths_dict.values()))
    # sort by size for plotting
    lengths = [lengths for lengths in X_contig_lengths_dict.values()]
    lengths.sort(reverse=True)
    contig_lengths_keys = {length : contig for contig, length in X_contig_lengths_dict.items()}
    contig_names_sorted_by_length = [contig_lengths_keys[length] for length in lengths]
    
    ################
    ### plotting ###
    ################
    colors_dict = {
        "abdomen" : "#E24D28", # fiery terracotta
        "head+thorax" : "#39676A", # stormy teal
        "separators" : "#ABABAB" # silver
    }
    tissues = ["abdomen","head+thorax"]

    fs = 25 # font size
    ps = 20 # point size

    # set figure aspect ratio
    aspect_ratio = 20 / 12
    height_pixels = 1200  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    
    ## setup contig labels, do the same as ReVis stacked histogram
    curr_contig_start = 0
    x_contig_coords = []
    x_contig_labels = []

    ax2 = ax.twinx() # axis that plots the contig labels
    ax2.spines['right'].set_visible(False)

    ax2.yaxis.set_ticks_position('none')  # Hide the default ticks
    ax2.yaxis.set_ticklabels([])  # Hide the default tick labels 
    
    plt.axvline(x=curr_contig_start, color=colors_dict["separators"], linestyle="--")
    for contig in contig_names_sorted_by_length:
        
        num_genes = len(coordinates_dict[contig])
        x_coord = [0.0]*num_genes
        abdomen = [0.0]*num_genes
        head_thorax = [0.0]*num_genes

        for i, transcript in enumerate(coordinates_dict[contig].keys()):
            x_coord[i] = coordinates_dict[contig][transcript]+curr_contig_start
            abdomen[i] = LFC_dict_abdomen[transcript]
            head_thorax[i] = LFC_dict_head_thorax[transcript]
        
        ## sort by X coordinate
        combined = sorted(zip(x_coord, abdomen, head_thorax))
        x_coord, abdomen, head_thorax = map(list, zip(*combined))
        
        ax.plot(x_coord, abdomen, color = colors_dict["abdomen"], alpha=0.5)
        ax.scatter(x_coord, abdomen, color = colors_dict["abdomen"], s=ps)
        ax.plot(x_coord, head_thorax, color = colors_dict["head+thorax"], alpha=0.5)
        ax.scatter(x_coord, head_thorax, color = colors_dict["head+thorax"], s=ps)
        
        x_contig_coords.append(curr_contig_start + X_contig_lengths_dict[contig]*0.5)
        x_contig_labels.append(contig)
        
        mean_min = curr_contig_start
        curr_contig_start += X_contig_lengths_dict[contig]
        mean_max = curr_contig_start

        plt.axvline(x=curr_contig_start, color=colors_dict["separators"], linestyle="--")

        print(f"{contig} : {num_genes} transcripts")
        ## plot mean and standard error of the mean
        a_mean = np.mean(abdomen)
        a_sem = sts.sem(abdomen)
        ht_mean = np.mean(head_thorax)
        ht_sem = sts.sem(head_thorax)

        if np.isnan(a_sem)==False:
            print(f"\t * abdomen mean log2FC: {a_mean:.3f}, SEM: {a_sem:.3f}")
            print(f"\t * head+thorax mean log2FC: {ht_mean:.3f}, SEM: {ht_sem:.3f}")
            ax.axhline(y=a_mean, xmin=mean_min/total_x_length, xmax=mean_max/total_x_length, color=colors_dict["abdomen"], linestyle="-") 
            ax.axhline(y=ht_mean, xmin=mean_min/total_x_length, xmax=mean_max/total_x_length, color=colors_dict["head+thorax"], linestyle="-") 
            ax.axhspan(ymin=a_mean-a_sem, ymax=a_mean+a_sem, xmin=mean_min/total_x_length, xmax=mean_max/total_x_length,color=colors_dict["abdomen"],  alpha=0.3, linewidth = 0)
            ax.axhspan(ymin=ht_mean-ht_sem, ymax=ht_mean+ht_sem, xmin=mean_min/total_x_length, xmax=mean_max/total_x_length, color=colors_dict["head+thorax"], alpha=0.3, linewidth = 0)


    ## set x axis labels and ticks correctly with contig names and stuff
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 0.01 else f'{x / 1e6:.0f} Mb'))
    ax2 = ax.secondary_xaxis('bottom')
    ax2.set_xticks(x_contig_coords)
    rot = 90
    if len(x_contig_labels[0])<2:
        rot = 0
    ax2.set_xticklabels(x_contig_labels, rotation=rot, fontsize=fs*0.7)
    # Adjust the position of the secondary x-axis
    ax2.spines['bottom'].set_position(('outward', 40))    
    ax2.xaxis.set_ticks_position('none')
    ax2.spines['bottom'].set_visible(False)
    # change tick fontsizes
    ax.tick_params(axis='x', labelsize=fs) 
    ax.tick_params(axis='y', labelsize=fs)
    ax2.tick_params(axis='x', labelsize=fs*0.8, rotation = rot)

    # ax.set_xlim(0-0.05*curr_contig_start, curr_contig_start+0.05*curr_contig_start)
    ax.set_xlim(0, curr_contig_start)

    ## plot for legend
    ax.plot(-1e3, 0, color = colors_dict["abdomen"], marker=".", markersize=ps, label = "abdomen")
    ax.plot(-1e3, 0, color = colors_dict["head+thorax"], marker=".", markersize=ps, label = "head+thorax")
    ax.legend(fontsize=fs)

    ax.set_ylabel(f"log2FC for female-male", fontsize = fs)
    ax2.set_xlabel(f"C. maculatus X chromosome", fontsize = fs, labelpad = 20)

    # layout (left, bottom, right, top)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = outfile.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile} and {filename_tr}")


def plot_sex_bias_bar_chart(summary_paths, annotation, X_list, sig_p_level = 0.001, minLFC = 1,outfile = ""):
    """
0 3  plot a bar chart to show proportion of significantly male or female biased genes on X or A chromosome
    """

    summary_data_abdomen = pd.read_csv(summary_paths["abdomen"], sep = "\t", index_col=False)
    summary_data_head_thorax = pd.read_csv(summary_paths["head_thorax"], sep = "\t", index_col=False)
    filtered_abdomen = summary_data_abdomen[pd.notna(summary_data_abdomen["logFC"])]
    filtered_head_thorax = summary_data_head_thorax[pd.notna(summary_data_head_thorax["logFC"])]

    ## get LFC information for every gene for plotting { transcriptID : LFC_float }
    ## the data is LFC female-male, so positive values are female-biased
    LFC_dict_abdomen = dict(zip(filtered_abdomen["geneID"], filtered_abdomen["logFC"]))
    LFC_dict_head_thorax = dict(zip(filtered_head_thorax["geneID"], filtered_head_thorax["logFC"]))
    assert sorted(list(LFC_dict_abdomen.keys())) == sorted(list(LFC_dict_head_thorax.keys())) # check that there is data for all transcripts in both

    ## get p-values
    pval_dict_abdomen = dict(zip(filtered_abdomen["geneID"], filtered_abdomen["FDR"]))
    pval_dict_head_thorax = dict(zip(filtered_head_thorax["geneID"], filtered_head_thorax["FDR"]))

    ## parse gene position data
    annotation_dict = gff.parse_gff3_general(annotation)
    head_thorax_summary_A = {
        "sig_male" : 0,
        "sig_female" : 0,
        "unbiased" : 0
    }
    head_thorax_summary_X = {
        "sig_male" : 0,
        "sig_female" : 0,
        "unbiased" : 0
    }
    abdomen_summary_A = {
        "sig_male" : 0,
        "sig_female" : 0,
        "unbiased" : 0
    }
    abdomen_summary_X = {
        "sig_male" : 0,
        "sig_female" : 0,
        "unbiased" : 0
    }

    # { contig : {} for contig in X_list} # dict with { contig : { transcriptID : [start,end], ... }  }
    for geneID in LFC_dict_abdomen.keys():
        gene = annotation_dict[geneID]

        if gene.contig in X_list:
            ## X-linked
            if pval_dict_abdomen[geneID]<sig_p_level:
                if LFC_dict_abdomen[geneID] > minLFC:
                    abdomen_summary_X["sig_female"]+=1 #female biased
                elif LFC_dict_abdomen[geneID] < -1*minLFC:
                    abdomen_summary_X["sig_male"]+=1
                else:
                    abdomen_summary_X["unbiased"]+=1
            else:
                abdomen_summary_X["unbiased"]+=1

            if pval_dict_head_thorax[geneID]<sig_p_level:
                if LFC_dict_head_thorax[geneID] > minLFC:
                    head_thorax_summary_X["sig_female"]+=1
                elif LFC_dict_head_thorax[geneID] < -1*minLFC:
                    head_thorax_summary_X["sig_male"]+=1
                else:
                    head_thorax_summary_X["unbiased"]+=1
            else:
                head_thorax_summary_X["unbiased"]+=1
        
        else:
            ## A-linked
            if pval_dict_abdomen[geneID]<sig_p_level:
                if LFC_dict_abdomen[geneID] > minLFC:
                    abdomen_summary_A["sig_female"]+=1
                elif LFC_dict_abdomen[geneID] < -1*minLFC:
                    abdomen_summary_A["sig_male"]+=1
                else:
                    abdomen_summary_A["unbiased"]+=1
            else:
                abdomen_summary_A["unbiased"]+=1

            if pval_dict_head_thorax[geneID]<sig_p_level:
                if LFC_dict_head_thorax[geneID] > minLFC:
                    head_thorax_summary_A["sig_female"]+=1
                elif LFC_dict_head_thorax[geneID] < -1*minLFC:
                    head_thorax_summary_A["sig_male"]+=1
                else:
                    head_thorax_summary_A["unbiased"]+=1
            else:
                head_thorax_summary_A["unbiased"]+=1

    print(f"head_thorax_summary_A: {head_thorax_summary_A}")
    print(f"head_thorax_summary_X: {head_thorax_summary_X}")
    print(f"abdomen_summary_A: {abdomen_summary_A}")
    print(f"abdomen_summary_X: {abdomen_summary_X}")

    ## make percentages for plot 
    sum_head_thorax_summary_A=sum(list(head_thorax_summary_A.values()))
    sum_head_thorax_summary_X=sum(list(head_thorax_summary_X.values()))
    sum_abdomen_summary_A=sum(list(abdomen_summary_A.values()))
    sum_abdomen_summary_X=sum(list(abdomen_summary_X.values()))
    head_thorax_percentage_A = {cat : val*100.0/sum_head_thorax_summary_A for cat,val in head_thorax_summary_A.items()}
    head_thorax_percentage_X = {cat : val*100.0/sum_head_thorax_summary_X for cat,val in head_thorax_summary_X.items()}
    abdomen_percentage_A = {cat : val*100.0/sum_abdomen_summary_A for cat,val in abdomen_summary_A.items()}
    abdomen_percentage_X = {cat : val*100.0/sum_abdomen_summary_X for cat,val in abdomen_summary_X.items()}

    ### plot stacked bar chart
    fs = 25
    width=0.4
    x_subtr=width*1.1/2.0
    x_coords = [1-x_subtr,1+x_subtr, 2-x_subtr,2+x_subtr]
    fig, ax = plt.subplots(1, 1, figsize=(16, 11)) 

    colors_dict ={
        "sig_female" : "#DC4141", # scarlet rush
        "sig_male" : "#7B8CE0", # wisteria blue
        "unbiased" : "#44455F", # vintage grape
    }
    labels_dict ={
        "sig_female" : "female biased",
        "sig_male" : "male biased", 
        "unbiased" : "unbiased",
    }
    bottom = [0,0,0,0]
    for DE_category in reversed(head_thorax_summary_A.keys()):

        y_coords = [head_thorax_percentage_A[DE_category],
                    head_thorax_percentage_X[DE_category],
                    abdomen_percentage_A[DE_category],
                    abdomen_percentage_X[DE_category]]
        ax.bar(x_coords, y_coords, width = width, bottom=bottom, label=labels_dict[DE_category], color= colors_dict[DE_category])
        
        bottom = [bottom[i]+y_coords[i] for i in range(len(y_coords))]

    # add numbers to bar chart
    vals=[]
    for DE_category in reversed(head_thorax_summary_A.keys()):
        # values in order left to right and bottom bar to top, so A_nonsig, X_nonsig, A_sig, X_sig
        vals.extend([f"{head_thorax_summary_A[DE_category]}", f"{head_thorax_summary_X[DE_category]}", f"{abdomen_summary_A[DE_category]}", f"{abdomen_summary_X[DE_category]}"])
    for i,bar in enumerate(ax.patches):
        ax.text(
            # Put the text in the middle of each bar. get_x returns the start
            # so we add half the width to get to the middle.
            bar.get_x() + bar.get_width() / 2,
            # Vertically, add the height of the bar to the start of the bar,
            # along with the offset.
            bar.get_height() + bar.get_y() -7,
            # This is actual value we'll show.
            vals[i],
            # Center the labels and style them a bit.
            ha='center',
            color='w',
            weight='bold',
            size=fs*0.9
        )

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<1 else f'{int(x)}%'))
    ax.set_xticks(ticks=x_coords, labels=["A: head+thorax","X: head+thorax","A: abdomen","X: abdomen"])
    ax.tick_params(axis='y', labelsize=fs)
    ax.tick_params(axis='x', labelsize=fs, rotation=90)

    ax.set_xlabel('')
    # ax.set_ylabel(f"Percentage of all genes", fontsize = fs)
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=fs)
    ax.tick_params(axis='x', labelsize=fs) 

    if minLFC == 0:
        leg_title = f"sig. threshold\np<{sig_p_level}"
    else:
        leg_title = f"sig. thresholds:\np<{sig_p_level} and\n|log2FC|>{minLFC}"
    plt.legend(reverse=True, fontsize=fs*0.8, loc="lower left", 
               title=leg_title, title_fontsize=fs*0.8)
    plt.title(f"Percent of male/female/unbiased genes", fontsize=fs*1.2)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = outfile.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile} and {filename_tr}")


def make_phylogeny_rank_dict(summary_file_df, min_p):
    """
    make a dict with rank orders like { 1 : [list, of, Log2FC, numbers] ,  2 [more, log2FC, numbers] , ... }
    """
    filtered_abdomen = summary_file_df[pd.notna(summary_file_df["LFC_abdomen"])]
    filtered_head_thorax = summary_file_df[pd.notna(summary_file_df["LFC_head+thorax"])]
    
    LFC_dict_abdomen = { i : [] for i in range(1,6)}
    LFC_dict_head_thorax = { i : [] for i in range(1,6)}

    if min_p == 0:
        for p_rank, log2FC_val  in zip(filtered_abdomen["level_most_dist_ortholog"], filtered_abdomen["LFC_abdomen"]):
            LFC_dict_abdomen[p_rank].append(log2FC_val)
        for p_rank, log2FC_val  in zip(filtered_head_thorax["level_most_dist_ortholog"], filtered_head_thorax["LFC_head+thorax"]):
            LFC_dict_head_thorax[p_rank].append(log2FC_val)
    else:
        for p_rank, log2FC_val,pval  in zip(filtered_abdomen["level_most_dist_ortholog"], filtered_abdomen["LFC_abdomen"], filtered_abdomen["FDR_pval_abdomen"]):
            if pval<min_p:
                LFC_dict_abdomen[p_rank].append(log2FC_val)
        for p_rank, log2FC_val,pval  in zip(filtered_head_thorax["level_most_dist_ortholog"], filtered_head_thorax["LFC_head+thorax"], filtered_head_thorax["FDR_pval_head+thorax"]):
            if pval<min_p:
                LFC_dict_head_thorax[p_rank].append(log2FC_val)

    return LFC_dict_abdomen, LFC_dict_head_thorax



def check_DE_phylogeny_rank_conserved(summary_paths_AX_dict:dict, outfile = "", abs_LFC=False, sig_p_threshold = 0, sep_MF=True):
    """
    check if genes with a higher phylogeny conservarion rank have higher log2FC values 
    if abs_LFC=True then do abs() around LFC to assess general sex bias and don't differentiate male-female contrast
    """

    summary_data_A = pd.read_csv(summary_paths_AX_dict["A"], sep = "\t", index_col=False)
    summary_data_X = pd.read_csv(summary_paths_AX_dict["X"], sep = "\t", index_col=False)

    LFC_dict_abdomen_A, LFC_dict_head_thorax_A = make_phylogeny_rank_dict(summary_data_A, min_p = sig_p_threshold)
    LFC_dict_abdomen_X, LFC_dict_head_thorax_X = make_phylogeny_rank_dict(summary_data_X, min_p = sig_p_threshold)

    if abs_LFC:
        LFC_lists_abdomen_A = [[abs(val) for val in vals_list] for vals_list in LFC_dict_abdomen_A.values()]
        LFC_lists_head_thorax_A = [[abs(val) for val in vals_list] for vals_list in LFC_dict_head_thorax_A.values()]
        LFC_lists_abdomen_X = [[abs(val) for val in vals_list] for vals_list in LFC_dict_abdomen_X.values()]
        LFC_lists_head_thorax_X = [[abs(val) for val in vals_list] for vals_list in LFC_dict_head_thorax_X.values()]
        y_label = f"|log2FC|"
    else:
        LFC_lists_abdomen_A = [vals_list for vals_list in LFC_dict_abdomen_A.values()]
        LFC_lists_head_thorax_A = [vals_list for vals_list in LFC_dict_head_thorax_A.values()]
        LFC_lists_abdomen_X = [vals_list for vals_list in LFC_dict_abdomen_X.values()]
        LFC_lists_head_thorax_X = [vals_list for vals_list in LFC_dict_head_thorax_X.values()]
        y_label = f"log2FC for female-male"
    
    if sig_p_threshold>0:
        y_label = f"{y_label}, p<{sig_p_threshold}"


    #### plot 2x2 boxplots
    
    fs = 25 # font size

    # set figure aspect ratio
    aspect_ratio = 18 / 12
    height_pixels = 1200  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(2,2,figsize=(width_pixels / 100, height_pixels / 100), dpi=100)

    def plot_DE_subplot(lists, row,col, fs, title, colors_dict, abs_logFC, lw=2):
        tick_labels = [f"{i+1} ({len(vals)})" for i,vals in enumerate(lists)]
        bp = ax[row,col].boxplot(lists, patch_artist=True)   

        # set axis labels
        ax[row,col].tick_params(axis='x', labelsize=fs) 
        ax[row,col].set_xticks(ticks = range(1,6), labels = tick_labels, fontsize=fs*0.8)
        ax[row,col].tick_params(axis='y', labelsize=fs)
        ax[row,col].set_title(title, fontsize=fs)

        ## modify boxplot colors
        for box in bp['boxes']:
            box.set(facecolor=colors_dict["fill"], edgecolor=colors_dict["edge"], linewidth=2)
        for median in bp['medians']:
            median.set(color=colors_dict['medians'], linewidth=lw)
        for whisker in bp['whiskers']:
            whisker.set(color=colors_dict['edge'], linestyle='-',linewidth=lw)
        for cap in bp['caps']:
            cap.set(color=colors_dict['edge'],linewidth=lw)
        for flier in bp['fliers']:
            flier.set(marker='.', markerfacecolor=colors_dict['edge'], markeredgecolor=colors_dict['edge'])
    
        if abs_logFC == False:
            ax[row,col].axhline(y=0, linestyle = "--", color = "black")


    def plot_DE_subplot_sep_MF(lists, row,col, fs, title, colors_dict, abs_logFC, lw=2):
        
        ## split into male- and female biased for separate boxplots
        sex_biased_lists = []
        for lst in lists:
            male_biased = [x for x in lst if x < 0]
            female_biased = [x for x in lst if x > 0]
            sex_biased_lists.extend([male_biased, female_biased])

        tick_labels = [f"{i // 2 + 1}\n({len(vals)})" for i,vals in enumerate(sex_biased_lists)]
        tick_pos = range(1,len(tick_labels)+1)
        bp = ax[row,col].boxplot(sex_biased_lists, patch_artist=True)

        # set axis labels
        tick_fs_factor = 0.6
        ax[row,col].set_xticks(ticks = tick_pos, labels = tick_labels, fontsize=fs*tick_fs_factor)
        ax[row,col].tick_params(axis='x', labelsize=fs*tick_fs_factor)
        ax[row,col].tick_params(axis='y', labelsize=fs*0.9)
        ax[row,col].set_title(title, fontsize=fs)

        ## modify boxplot colors
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
    
        if abs_logFC == False:
            ax[row,col].axhline(y=0, linestyle = "--", color = "black")

    colors = {
        "fill" : "#495E83", # dusk blue
        "edge" : "#374C6E", # dusk blue darker
        "medians" : "#A7CCED", # icy blue
        "lines" : "#7D93B5", # lavender grey
        "FB_fill" : "#AB354A", # cherry rose
        "FB_edge" : "#771C2C", # dark amaranth
        "FB_medians" : "#EA9AA9", # cotton candy
        # "FB_lines" : "#7D93B5" # lavender grey
    }

    if sep_MF:
        print(f" separate male and female sex biased genes")
        plot_DE_subplot_sep_MF(LFC_lists_abdomen_A, row=0, col=0, fs=fs, title=f"Autosomes: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot_sep_MF(LFC_lists_head_thorax_A, row=0, col=1, fs=fs, title=f"Autosomes: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot_sep_MF(LFC_lists_abdomen_X, row=1, col=0, fs=fs, title=f"X-chromosome: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot_sep_MF(LFC_lists_head_thorax_X, row=1, col=1, fs=fs, title=f"X-chromosome: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        fig.supxlabel(f"conservation rank\n(number of genes)", fontsize = fs)
    else:
        plot_DE_subplot(LFC_lists_abdomen_A, row=0, col=0, fs=fs, title=f"Autosomes: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot(LFC_lists_head_thorax_A, row=0, col=1, fs=fs, title=f"Autosomes: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot(LFC_lists_abdomen_X, row=1, col=0, fs=fs, title=f"X-chromosome: abdomen", colors_dict=colors, abs_logFC=abs_LFC)
        plot_DE_subplot(LFC_lists_head_thorax_X, row=1, col=1, fs=fs, title=f"X-chromosome: head+thorax", colors_dict=colors, abs_logFC=abs_LFC)
        fig.supxlabel(f"conservation rank, and (number of genes)", fontsize = fs)
    fig.supylabel(y_label, fontsize = fs)

    # layout (left, bottom, right, top)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # transparent background
    plt.savefig(outfile, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = outfile.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile} and {filename_tr}")



if __name__ == "__main__":
    
    username = "milena"
    Cmac_annotation = f"/Users/{username}/work/native_annotations/all_native_annot/C_maculatus_superscaffolded_LomeRNA_braker_isoform_filtered.gff"
    # milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse/Cmac_Lome_diverse/braker/braker_isoform_filtered.gff C_maculatus_superscaffolded_LomeRNA_braker_isoform_filtered.gff
    Cmac_assembly = f"/Users/{username}/work/assemblies_masked_uniform/C_maculatus_superscaffolded_genomic_fasta.masked.fai" # I only need contig lengths so use assembly index
    summary_paths = get_summary_paths(username=username)
    DE_paths = get_DE_paths(username=username)
    Cmac_X_contigs_list = get_Cmac_superscaffolded_XY_contigs()["X"]

    if False:
        plot_dosage_compensation(summary_paths=DE_paths, annotation=Cmac_annotation, assembly_index=Cmac_assembly, 
            X_list=Cmac_X_contigs_list, 
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/X_sex_bias.png")

    if False:
        plot_sex_bias_bar_chart(summary_paths=DE_paths, annotation=Cmac_annotation, X_list=Cmac_X_contigs_list, 
            sig_p_level = 0.001, minLFC = 0.5,
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/all_sex_bias_proportion.png")

    if True:
        abs_logFC = False
        check_DE_phylogeny_rank_conserved(summary_paths_AX_dict=summary_paths,
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/conservation_rank_all_sex_bias_proportion.png",
            abs_LFC=abs_logFC, sig_p_threshold=0)
        check_DE_phylogeny_rank_conserved(summary_paths_AX_dict=summary_paths,
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/conservation_rank_sig_sex_bias_proportion.png",
            abs_LFC=abs_logFC, sig_p_threshold=0.05)
