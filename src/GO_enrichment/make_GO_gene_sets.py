"""
filter lists of gene sets from the summary tables
and do some plots to see which ones to select
"""

import pandas as pd

from matplotlib import pyplot as plt
from matplotlib_venn import venn3


def get_full_table_path(username="miltr339"):
    out_dict = {
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_A_chr.tsv",
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/DE_summary_table_X_chr.tsv",
    }
    return out_dict


def test_geneID_overlap(table_path:str, venn_name = f"venn", venn_title = ""):
    """
    Make a venn diagramm of the three species to see how the geneIDs overlap between them
    """

    df = pd.read_csv(table_path, sep="\t")
    species_list = sorted(list(set(list(df["other_species"]))))

    geneIDs_lists = { species : [] for species in species_list}
    geneIDs_pos_sel_lists = { species : [] for species in species_list}
    
    for species in species_list:
        species_df = df[df["other_species"]==species]
        geneIDs_lists[species] = set(species_df["focal_transcript"])

        pos_sel_species_df = species_df[species_df["positive_selection"] == True]
        geneIDs_pos_sel_lists[species] = set(pos_sel_species_df["focal_transcript"])
    
        print(f"{species} ({len(geneIDs_lists[species])})")

    ### plot all genes overlap
    venn3(subsets=(geneIDs_lists[species_list[0]], geneIDs_lists[species_list[1]], geneIDs_lists[species_list[2]]), 
    set_labels=(species_list[0], species_list[1], species_list[2]),
    set_colors=("orange", "blue", "red"), alpha=0.7)

    plt.tight_layout()
    plt.title(f"{venn_title} : all genes")
    plt.savefig(venn_name, dpi = 300, transparent = False)
    plt.clf()
    print(f"venn diagramm saved here: {venn_name}")
    ## all list is the intersection
    all_intersection = list(geneIDs_lists[species_list[0]].intersection(geneIDs_lists[species_list[1]], geneIDs_lists[species_list[2]]))
    ## all list union
    all_union = list(set([geneID for species_list in geneIDs_lists.values() for geneID in species_list]))

    ### plot positively selected genes overlap
    venn3(subsets=(geneIDs_pos_sel_lists[species_list[0]], geneIDs_pos_sel_lists[species_list[1]], geneIDs_pos_sel_lists[species_list[2]]), 
    set_labels=(species_list[0], species_list[1], species_list[2]),
    set_colors=("orange", "blue", "red"), alpha=0.7)
    plt.tight_layout()
    plt.title(f"{venn_title} : positively selected")
    venn_name_sig = venn_name.replace(".png", "_pos_sel_genes.png")
    plt.savefig(venn_name_sig, dpi = 300, transparent = False)
    plt.clf()
    print(f"venn diagramm saved here: {venn_name_sig}")
    ## sig list is the intersection
    sig_intersection = list(geneIDs_pos_sel_lists[species_list[0]].intersection(geneIDs_pos_sel_lists[species_list[1]], geneIDs_pos_sel_lists[species_list[2]]))
    ## all list union
    sig_union = list(set([geneID for species_list in geneIDs_pos_sel_lists.values() for geneID in species_list]))



if __name__ == "__main__":

    username = "miltr339"
    tables_dict = get_full_table_path(username=username)
    outfiles_dir = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/GO_enrichment"

    for chr in ["X", "A"]:
        print(f"///////////////// {chr} /////////////////")
        test_geneID_overlap(table_path= tables_dict[chr], 
        venn_name=f"{outfiles_dir}/geneID_overlap_Venn_{chr}.png", 
        venn_title = f"geneID overlap {chr}")