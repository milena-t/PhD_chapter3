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
from analyze_site_classes import read_site_classes,get_pairs_from_summary,species_names_from_pair
import os
import pandas as pd

def get_lookup_tables(username="miltr339"):
    """
    tables that contain the information about which orthologs have which transcripts in them
    """
    tables={
        "X" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/ortholog_IDs_X_transcript_IDs_association.txt",
        "A" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/ortholog_IDs_A_transcript_IDs_association.txt",
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
            "sites" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_A_BH_corrected.txt",
            "dNdS" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/dNdS_by_ortholog_A-linked_updated_species.txt",

            },
        "X" : {
            "sites" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/site_classes_summary_X_BH_corrected.txt",
            "dNdS" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/dNdS_by_ortholog_X-linked_updated_species.txt",

            }
    }
    return tables


def get_BRH_paths(username="miltr339"):
    """
    blast best reciprocal hits of Cmac and all other species, to determine how far back the gene of interest is conserved
    """
    tables = {
        "D_melanogaster" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_D_melanogaster_BRHs.tsv",
        "T_castaneum" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_T_castaneum_BRHs.tsv",
        "T_freemani" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_T_freemani_BRHs.tsv",
        "C_magnifica" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_C_magnifica_BRHs.tsv",
        "C_septempunctata" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_C_septempunctata_BRHs.tsv",
        "A_obtectus" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_A_obtectus_BRHs.tsv",
        "B_siliquastri" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_B_siliquastri_BRHs.tsv",
        "C_chinensis" : f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables/C_maculatus_C_chinensis_BRHs.tsv",
    }
    dist_dict = {
        "D_melanogaster" : 5,
        "T_castaneum" : 4,
        "T_freemani" : 4,
        "C_magnifica" : 4,
        "C_septempunctata" : 4,
        "A_obtectus" : 3,
        "B_siliquastri" : 2,
        "C_chinensis" : 1,
    }
    return tables, dist_dict


def read_ortholog_dNdS(summary_path, excl_list = [], pair_split_string="_"):
    """
    read the dNdS by-ortholog summary file into a dictionary like 
    dict = {
        pair : {
            ortholog_num : {
                dN : float,
                dS : float,
                dNdS : float
            }
        }
    }
    """
    if not os.path.isfile(summary_path):
        raise RuntimeError(f"FILE: {summary_path} does not exist!")

    pairs_list = get_pairs_from_summary(summary_path, excl_list=excl_list)
    
    out_dict = {pair : {} for pair in pairs_list}
    no_dNdS = {pair : [] for pair in pairs_list}
    print(f" *  {len(out_dict)} unique pairs")


    with open(summary_path, "r") as summary_file:
        for line in summary_file.readlines():
            try:
                ortholog_dir, values = line.strip().split(":")
            except:
                raise RuntimeError(f"line cannot be parsed! \n{line}")
                # continue
            try:
                pair_name,num_dir= ortholog_dir.split(pair_split_string)
                pair_number = int(num_dir.split("_dNdS")[0])
            except:
                if ortholog_dir[-4:] == ".out" or ortholog_dir[-4:] == ".log" or len(ortholog_dir)<5:
                    continue
                else:
                    raise RuntimeError(f"cannot parse pair number in {pair_name}, from filepath: {ortholog_dir}")

            species1,species2=species_names_from_pair(pair_name)
            pair = f"{species1}_{species2}"
            if pair not in pairs_list:
                continue # skip species that are excluded by excl_list above

            try:
                values_dict = { val.split("=")[0] : val.split("=")[1]  for val in values.split(",")}
            except:
                no_dNdS[pair].append(ortholog_dir)
                continue

            out_dict[pair][pair_number]=values_dict

    return out_dict, no_dNdS



def read_LFC(summary_path):
    """
    read the dNdS by-ortholog summary file into a dictionary like 
    dict = {
        geneID : {
            LFC : float,
            FDR : float # FDR corrected p-value
        }
    }
    """
    if not os.path.isfile(summary_path):
        raise RuntimeError(f"FILE: {summary_path} does not exist!")

    out_dict = {}

    with open(summary_path, "r") as summary_file:
        for line in summary_file.readlines()[1:]: # ignore header line
            try:
                geneID,logFC,logCPM,LR_num,PValue,FDR_p = line.strip().split("\t")
            except:
                raise RuntimeError(f"line cannot be parsed! \n{line}")
                # continue
            
            out_dict[geneID]={"LFC" : logFC, "FDR" : FDR_p}

    return out_dict


def read_BRHs(brhs_paths_dict):
    """
    read all Cmac BRHs into a dict like
    {
        species1 : [Cmac,transcript,IDs,that,have,a,BRH,with,this,species,...],
        species2 : [more,Cmac,transcript,IDs,...]
    }
    also mind that the IDs are transcripts and have the tailing '_1'
    """
    out_dict = {}
    for pair_species,path_to_table in brhs_paths_dict.items():
        brh_df = pd.read_csv(path_to_table, sep="\t", names=["C_maculatus", pair_species])
        out_dict[pair_species] = [match[:-2] for match in brh_df['C_maculatus'].values.tolist()]
        
    return out_dict


def make_summary_table(lookup_table_path, DE_paths_dict, paml_paths_dict, chromosome, brhs_lists, distances_dict, focal_species = "C_maculatus", excl_species=["T_castaneum", "C_septempunctata"], outfile_name = "orthologs_DE_summary.txt"):
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

    header = [
        "focal_species",
        "other_species",
        "ortholog_number",
        "focal_transcript",
        "other_transcript",
        "positive_selection",
        "dN",
        "dS",
        "dN/dS",
        "focal_gene_ID",
        "LFC_abdomen",
        "FDR_pval_abdomen",
        "LFC_head+thorax",
        "FDR_pval_head+thorax",
        "level_most_dist_ortholog"
    ]
    header_line = "\t".join(header)
    header_line = f"{header_line}\n"

    ## read input data
    ## positively selected sites
    sites_path = paml_paths_dict["sites"]
    if not os.path.isfile(sites_path):
        raise RuntimeError(f"FILE: {sites_path} does not exist!")
    site_classes_dict,_ = read_site_classes(sites_path, excl_list=excl_species, read_for_table = True)
    print(list(site_classes_dict.keys()))

    ## dNdS
    dNdS_path = paml_paths_dict["dNdS"]
    if not os.path.isfile(dNdS_path):
        raise RuntimeError(f"FILE: {dNdS_path} does not exist!")
    pair_dir_split = f"_{chromosome}-linked_ortholog_" # parse the dNdS directory analysis naming scheme which is like 'B_siliquastri_C_maculatus_X-linked_ortholog_68'
    dNdS_values_dict,_ = read_ortholog_dNdS(dNdS_path, excl_list=excl_species,pair_split_string=pair_dir_split)

    ## Log-fold change
    # abdomen
    DE_abdomen = read_LFC(DE_paths_dict["DE_abdomen"])
    # head+thorax
    DE_head_thorax = read_LFC(DE_paths_dict["DE_head_thorax"])

    with open(lookup_table_path, "r") as ortholog_lookup_table, open(outfile_name, "w") as outfile:
        
        outfile.write(header_line)
        ortholog_lines = [line.strip() for line in ortholog_lookup_table.readlines() if not any(e in line for e in excl_species)]

        ## loop through lookup table and add all the stuff from the different input files to each ortholog line
        for line in tqdm(ortholog_lines):
            pair_dir,transcripts = line.split(":")
            if any(e in pair_dir for e in excl_species):
                continue
            try:
                pair,ortholog_num = pair_dir.split(pair_dir_split)
                species1,species2 = species_names_from_pair(pair)
                transcript1,transcript2 = transcripts.split(",")
            except:
                print(f"ERROR: {line}")
                continue
            if species1 != focal_species and species2 != focal_species:
                continue
            
            if species1 == focal_species:
                line_focal_species = species1
                other_species = species2
                focal_transcript = transcript1
                other_transcript = transcript2
            elif species2 == focal_species:
                line_focal_species = species2
                other_species = species1
                focal_transcript = transcript2
                other_transcript = transcript1
            
            ### Differential expression
            focal_gene = focal_transcript.split(".t")[0] ## !!! this only works for thge BRAKER annotation, where the format is g1.t1 for the transcript and then g1 for the parent gene
            try:
                abd_LFC = DE_abdomen[focal_gene]["LFC"]
                abd_FDR = DE_abdomen[focal_gene]["FDR"]
                ht_LFC = DE_head_thorax[focal_gene]["LFC"]
                ht_FDR = DE_head_thorax[focal_gene]["FDR"]
            except:
                # not expressed
                abd_LFC = "NaN"
                abd_FDR = 1
                ht_LFC = "NaN"
                ht_FDR = 1

            ### positive selection
            try: 
                ortholog = site_classes_dict[pair][int(ortholog_num)]
                pos_sel = ortholog.sig_pos_selection
            except:
                # if this is absent, then likely paml didn't work because the gene does not have a correct start/stop codon or some other structural issue
                # also gives "not_found" for dNdS
                pos_sel = "NaN"

            ### dNdS
            dNdS_ortholog = dNdS_values_dict[pair][int(ortholog_num)]
            try:
                dN = dNdS_ortholog["dN"]
                dS = dNdS_ortholog["dS"]
                dNdS = dNdS_ortholog["dNdS"]
            except:
                raise RuntimeError(f"could not parse dNdS from {pair} ortholog {ortholog_num}! \n\tthis is what was read: {dNdS_ortholog}")

            ### most distant ortholog of the focal gene
            ortholog_species = [species for species,ortholog_list in brhs_lists.items() if focal_transcript in ortholog_list]
            ortholog_levels = [distances_dict[species] for species in ortholog_species]
            max_dist_level = max(ortholog_levels)

            ## make outfile line
            outfile_line_orthologs = f"{line_focal_species}\t{other_species}\t{ortholog_num}\t{focal_transcript}\t{other_transcript}"
            outfile_line_paml = f"{pos_sel}\t{dN}\t{dS}\t{dNdS}"
            outfile_line_LFC = f"{focal_gene}\t{abd_LFC}\t{abd_FDR}\t{ht_LFC}\t{ht_FDR}"
            outfile_line_all = f"{outfile_line_orthologs}\t{outfile_line_paml}\t{outfile_line_LFC}\t{max_dist_level}\n"

            outfile.write(outfile_line_all)
            # break

    print(f"outfile written to: {outfile_name}")     

if __name__=="__main__":
    
    ### load data
    username="milena"
    outdir_tables=f"/Users/{username}/work/PhD_code/PhD_chapter3/data/DE_analysis/paml_summary_tables"
    lookup_tables_dict = get_lookup_tables(username=username)
    DE_paths_dict = get_DE_paths(username=username)
    paml_paths_dict = get_paml_paths(username=username)
    brh_tables_paths, distances_dict = get_BRH_paths(username=username)
    brhs_list = read_BRHs(brh_tables_paths)
    
    chromosomes = ["A", "X"]

    if True:
        for chromosome in chromosomes:
            print(f"\n//////////////////////// {chromosome} ////////////////////////")
            make_summary_table(lookup_table_path=lookup_tables_dict[chromosome], 
                DE_paths_dict=DE_paths_dict, 
                paml_paths_dict=paml_paths_dict[chromosome], 
                chromosome=chromosome, 
                focal_species = "C_maculatus", excl_species=["T_castaneum", "C_septempunctata"],
                brhs_lists = brhs_list, distances_dict = distances_dict,
                outfile_name=f"{outdir_tables}/DE_summary_table_{chromosome}_chr.tsv")
            
