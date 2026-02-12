"""
make a summary table for every ortholog from a C. maculatus perspective. includes 
* all Cmac-orthologs 
    * cmac transcript ID
    * partner transcript ID
* dN value
* dS value
* dN/dS
* positive selection? (binary)
* Cmac gene ID
    * normalized read count (log-transformed)
    * LFC abdomen
    * LFC head+thorax
* other orthologs the cmac gene may have in species outside bruchids
"""

from tqdm import tqdm
from analyze_site_classes import read_site_classes
import os

def get_lookup_tables(username="miltr339"):
    """
    tables that contain the information about which orthologs have which transcripts in them
    """
    tables={
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/ortholog_IDs_X_transcript_IDs_association.txt",
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/ortholog_IDs_A_transcript_IDs_association.txt",
    }
    return tables

def get_DE_paths(username="miltr339"):
    """
    all output files for the differential expression
    """
    tables = {
       "counts" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_gene_counts_edgeR_normalized.txt",
       "DE_abdomen" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_sex_DE_abdomen_edgeR.txt",
       "DE_head_thorax" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_sex_DE_head_thorax_edgeR.txt",
    }
    return tables

def get_paml_paths(username="miltr339"):
    """
    all outputs from the paml branch and site model analyses

    rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_A/site_classes_summary_A-linked.txt /Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/site_classes_summary_A-linked.txt
    rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_X/site_classes_summary_X-linked.txt /Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/site_classes_summary_X-linked.txt
    rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_A_branch_model/dNdS_by_ortholog_A-linked_updated_species.txt /Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_by_ortholog_A-linked_updated_species.txt
    rsync -azP milenatr@pelle.uppmax.uu.se:/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_results_X_branch_model/dNdS_by_ortholog_X-linked_updated_species.txt /Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_by_ortholog_X-linked_updated_species.txt
    """
    tables = {
        "A" : {
            "sites" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/site_classes_summary_A-linked.txt",
            "dNdS" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_by_ortholog_A-linked_updated_species.txt",
            },
        "X" : {
            "sites" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/site_classes_summary_X-linked.txt",
            "dNdS" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/dNdS_by_ortholog_X-linked_updated_species.txt",
            }
    }
    return tables



def make_summary_table(lookup_table_path, DE_paths_dict, paml_paths_dict, chromosome, focal_species = "C_maculatus", excl_species=["T_castaneum", "C_septempunctata"], outfile_name = "orthologs_DE_summary.txt"):
    """
    make the summary table for either A or X chromosome with the columns
    * all Cmac-orthologs 
        * cmac transcript ID
        * partner transcript ID
    * dN value
    * dS value
    * dN/dS
    * positive selection? (binary)
    * Cmac gene ID
        * normalized read count (log-transformed)
        * LFC abdomen
        * LFC head+thorax
    """

    ## read input data
    sites_path = paml_paths_dict["sites"]
    if not os.path.isfile(sites_path):
        raise RuntimeError(f"FILE: {sites_path} does not exist!")
    site_classes_dict,_ = read_site_classes(sites_path, excl_list=excl_species, read_for_table = True)
    print(list(site_classes_dict.keys())[:10])
    
    dNdS_path = paml_paths_dict["dNdS"]
    if not os.path.isfile(dNdS_path):
        raise RuntimeError(f"FILE: {dNdS_path} does not exist!")
    pair_dir_split = f"_{chromosome}-linked_ortholog_" # parse the dNdS directory analysis naming scheme which is like 'B_siliquastri_C_maculatus_X-linked_ortholog_68'

    with open(lookup_table_path, "r") as ortholog_lookup_table, open(outfile_name, "w") as outfile:

        ortholog_lines = [line.strip() for line in ortholog_lookup_table.readlines() if not any(e in line for e in excl_species)]

        ## loop through lookup table and add all the stuff from the different input files to each ortholog line
        for line in ortholog_lines:
            pair_dir,transcripts = line.split(":")
            if any(e in pair_dir for e in excl_species):
                continue
            try:
                pair,ortholog_num = pair_dir.split(pair_dir_split)
                g1,s1,g2,s2=pair.split("_")
                species1 = f"{g1}_{s1}"
                species2 = f"{g2}_{s2}"
            except:
                continue

            print(f"{pair}: ortholog {ortholog_num}")
            print(f"site class")

            ortholog = site_classes_dict[pair][int(ortholog_num)]
            print(ortholog)
            print(ortholog.sig_pos_selection)

            break


if __name__=="__main__":
    
    ### load data
    username="miltr339"
    outdir_tables=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables"
    lookup_tables_dict = get_lookup_tables(username=username)
    DE_paths_dict = get_DE_paths(username=username)
    paml_paths_dict = get_paml_paths(username=username)

    # chromosomes = ["A", "X"]
    chromosomes = ["X"]

    for chromosome in chromosomes:
        print(f"\n//////////////////////// {chromosome} ////////////////////////")
        make_summary_table(lookup_table_path=lookup_tables_dict[chromosome], 
            DE_paths_dict=DE_paths_dict, 
            paml_paths_dict=paml_paths_dict[chromosome], 
            chromosome=chromosome, 
            focal_species = "C_maculatus", excl_species=["T_castaneum", "C_septempunctata"],
            outfile_name=f"{outdir_tables}/DE_summary_table.tsv")
        
