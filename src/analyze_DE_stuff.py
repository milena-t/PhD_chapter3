"""
analyze the info from the big summary tables generated with PhD_chapter3/src/make_DE_analysis_table.py
"""

import parse_gff as gff
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter


def get_summary_paths(username = "miltr339"):
    summary_dict = {
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_A_chr.tsv",
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_X_chr.tsv",
    }
    return summary_dict

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



def plot_dosage_compensation(summary_path, annotation, assembly_index, X_list, outfile = ""):
    """
    plot the log2FC of all genes on the X in both head+thorax and abdominal tissues on the y,
    the X-axis is the position of the gene on the X
    """


    summary_data = pd.read_csv(summary_path, sep = "\t")

    
    ## get LFC information for every gene for plotting { transcriptID : LFC_float }
    filtered = summary_data[pd.notna(summary_data["LFC_abdomen"])]
    LFC_dict_abdomen = dict(zip(filtered["focal_transcript"], filtered["LFC_abdomen"]))
    LFC_dict_head_thorax = dict(zip(filtered["focal_transcript"], filtered["LFC_head+thorax"]))

    assert sorted(list(LFC_dict_abdomen.keys())) == sorted(list(LFC_dict_head_thorax.keys())) # check that there is data for all transcripts in both

    ## parse gene position data by contig for plotting
    annotation_dict = gff.parse_gff3_general(annotation)
    coordinates_dict = { contig : {} for contig in X_list} # dict with { contig : { transcriptID : [start,end], ... }  }
    for transcriptID in LFC_dict_abdomen.keys():
        transcript = annotation_dict[transcriptID]
        if transcript.contig not in X_list:
            print(transcript)
            continue
        # coordinates_dict[transcript.contig][transcript.feature_id] = [transcript.start, transcript.end]
        coordinates_dict[transcript.contig][transcript.feature_id] = np.mean([transcript.start, transcript.end])

    ## get the X contig lenghts for plotting the X axis correctly
    print(f"get X contig lengths ...")
    X_contig_lengths_dict = get_contig_lengths(contigs_list=X_list,assembly_index=assembly_index)
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
        print(f"{contig} : {num_genes} transcripts")
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
        curr_contig_start += X_contig_lengths_dict[contig]
        plt.axvline(x=curr_contig_start, color=colors_dict["separators"], linestyle="--")

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

    ax.set_xlim(0-0.05*curr_contig_start, curr_contig_start+0.05*curr_contig_start)

    ## plot for legend
    ax.plot(-1e3, 0, color = colors_dict["abdomen"], marker=".", markersize=ps, label = "abdomen")
    ax.plot(-1e3, 0, color = colors_dict["head+thorax"], marker=".", markersize=ps, label = "head+thorax")
    ax.legend(fontsize=fs)

    ax.set_ylabel(f"log2FC for female-male", fontsize = fs)
    ax2.set_xlabel(f"C. maculatus X chromosome", fontsize = fs, labelpad = 20)

    plt.tight_layout(rect=[0, 0.05, 1, 1])

    plt.savefig(outfile, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile}")




if __name__ == "__main__":
    
    username = "miltr339"
    Cmac_annotation = f"/Users/{username}/work/native_annotations/all_native_annot/C_maculatus_superscaffolded_LomeRNA_braker_isoform_filtered.gff"
    # milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse/Cmac_Lome_diverse/braker/braker_isoform_filtered.gff
    Cmac_assembly = f"/Users/{username}/work/assemblies_masked_uniform/C_maculatus_superscaffolded_genomic_fasta.masked.fai" # I only need contig lengths so use assembly index
    summary_paths = get_summary_paths()
    Cmac_X_contigs_list = get_Cmac_superscaffolded_XY_contigs()["X"]

    if True:
        plot_dosage_compensation(summary_path=summary_paths["X"], annotation=Cmac_annotation, assembly_index=Cmac_assembly, X_list=Cmac_X_contigs_list, 
            outfile=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/X_sex_bias.png")

