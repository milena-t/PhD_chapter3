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


def test_geneID_overlap(table_path:str,chr = "", outdir = f"", venn_name = f"venn", venn_title = ""):
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
        all_IDs_outstr = ",".join(list(geneIDs_lists[species]))
        all_outname=f"{outdir}geneIDs_{chr}_{species}_all.txt"
        with open(all_outname, "w") as outfile:
            outfile.write(f"geneID,pos_sel\n")

            pseudo_pval_dict = {True : 1 , False: 0} # make numeric for topGO in R
            for geneID,pos_sel in zip(species_df["focal_transcript"], species_df["positive_selection"]):
                if pos_sel == True or pos_sel ==False:
                    outfile.write(f"{geneID},{pseudo_pval_dict[pos_sel]}\n")
                else:
                    pass
            print(f"list written to: {all_outname}")

        pos_sel_species_df = species_df[species_df["positive_selection"] == True]
        geneIDs_pos_sel_lists[species] = set(pos_sel_species_df["focal_transcript"])
        sig_IDs_outstr = ",".join(list(geneIDs_pos_sel_lists[species]))
        sig_outname=f"{outdir}geneIDs_{chr}_{species}_sig.txt"
        with open(sig_outname, "w") as outfile:
            outfile.write(sig_IDs_outstr)
            print(f"list written to: {sig_outname}")

        print(f"{species} ({len(geneIDs_lists[species])})")

    ### plot all genes overlap
    venn3(subsets=(geneIDs_lists[species_list[0]], geneIDs_lists[species_list[1]], geneIDs_lists[species_list[2]]), 
    set_labels=(species_list[0], species_list[1], species_list[2]),
    set_colors=("orange", "blue", "red"), alpha=0.7)

    plt.title(f"{venn_title} : all genes")
    plt.savefig(venn_name, dpi = 300, transparent = False)
    plt.clf()
    print(f"venn diagramm saved here: {venn_name}")
    ## all list is the intersection
    all_intersection = list(geneIDs_lists[species_list[0]].intersection(geneIDs_lists[species_list[1]], geneIDs_lists[species_list[2]]))
    ## all list union
    all_union = list(set([geneID for species_list in geneIDs_lists.values() for geneID in species_list]))
    all_outname = f"{outdir}all_geneIDs_{chr}_all_species_union.txt"
    with open(all_outname, "w") as outfile:
        outfile.write(",".join(all_union))
        print(f"list written to: {all_outname}")


    ### plot positively selected genes overlap
    venn3(subsets=(geneIDs_pos_sel_lists[species_list[0]], geneIDs_pos_sel_lists[species_list[1]], geneIDs_pos_sel_lists[species_list[2]]), 
    set_labels=(species_list[0], species_list[1], species_list[2]),
    set_colors=("orange", "blue", "red"), alpha=0.7)
    plt.title(f"{venn_title} : positively selected")
    venn_name_sig = venn_name.replace(".png", "_pos_sel_genes.png")
    plt.savefig(venn_name_sig, dpi = 300, transparent = False)
    plt.clf()
    print(f"venn diagramm saved here: {venn_name_sig}")
    ## sig list is the intersection
    sig_intersection = list(geneIDs_pos_sel_lists[species_list[0]].intersection(geneIDs_pos_sel_lists[species_list[1]], geneIDs_pos_sel_lists[species_list[2]]))
    ## all list union
    sig_union = list(set([geneID for species_list in geneIDs_pos_sel_lists.values() for geneID in species_list]))
    sig_outname = f"{outdir}sig_geneIDs_{chr}_all_species_union.txt"
    with open(sig_outname, "w") as outfile:
        outfile.write(",".join(sig_union))
        print(f"list written to: {sig_outname}")


if __name__ == "__main__":

    username = "miltr339"
    tables_dict = get_full_table_path(username=username)
    outfiles_dir = f"/Users/{username}/work/PhD_code/PhD_chapter3/data/GO_enrichment"

    for chr in ["X", "A"]:
        print(f"///////////////// {chr} /////////////////")
        test_geneID_overlap(table_path= tables_dict[chr], chr=chr, outdir = f"{outfiles_dir}/",
        venn_name=f"{outfiles_dir}/geneID_overlap_Venn_{chr}.png", 
        venn_title = f"geneID overlap {chr}")