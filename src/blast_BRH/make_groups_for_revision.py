"""
For the revision I will test paml site models with a four-species phylogeny from the perspective of C. maculatus.
I will find all other Bruchini orhtologs of any C. maculatus gene, and then check if they are also each other's partners to identify groups of four 1-to-1 orthologs
I will then use them to run paml site models and identify positively selected orthologs
"""

import pandas as pd 

def get_paths(username = "miltr339"):
    dirname=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/"
    files_dict = {
        ("A_obtectus","B_siliquastri") : f"{dirname}A_obtectus_B_siliquastri_BRH.tsv",
        ("A_obtectus","C_maculatus") : f"{dirname}A_obtectus_C_maculatus_BRH.tsv",
        ("A_obtectus","C_chinensis") : f"{dirname}A_obtectus_C_chinensis_BRH.tsv",
        ("B_siliquastri","C_maculatus") : f"{dirname}B_siliquastri_C_maculatus_BRH.tsv",
        ("B_siliquastri","C_chinensis") : f"{dirname}B_siliquastri_C_chinensis_BRH.tsv",
        ("C_chinensis","C_maculatus") : f"{dirname}C_chinensis_C_maculatus_BRH.tsv",
    }
    return files_dict


class SpeciesOrtholog:
    """
    A species ortholog with a key of one geneID, species ident and a dict of orthologs { geneID : species }
    """
    def __init__(self, geneID:str, species:str, orthologs:dict, chromosome:str) -> None:
        self.geneID = geneID
        self.species=species
        self.orthologs=orthologs
        self.chromosome=chromosome

    def add_ortholog(self, key_value_tuple:tuple):
        orthologID,species=key_value_tuple
        self.orthologs[orthologID]= species

    def test_ortholog_presence(self, species:str, orthologID:str):
        if orthologID in self.orthologs:
            if self.orthologs[orthologID] == species:
                return True
            else:
                return False
        else:
            return False
    
    def __str__(self):
        return(
            f"""
gene {self.geneID} in species: {self.species} (on chromosome {self.chromosome})
\thas {len(self.orthologs)} orthologs: {self.orthologs}
            """
        )

class OrthoGroup:
    """
    group of orthologs that is 1-to-1 in all species as a dictionary like { species : geneID }
    """
    def __init__(self, groupID:int, species_dict:dict, chromosome:str) -> None:
        self.groupID=groupID
        self.species_dict=species_dict
        self.chromosome=chromosome
    

def get_pair_partner(pair_tuple, species):
    pair_list = list(pair_tuple)
    pair_list.remove(species)
    assert len(pair_list)==1
    return pair_list[0]

def get_all_orthologs(pair_paths_dir, species_list, outfile):
    """
    write a list of 1-to-1 orthologs across all species in species_list into an outfile
    """

    ## read all orthologs into class
    orthologs_dict = { species : {} for species in species_list}
    for species in species_list:
        species_orthologs = {}
        print(f"================= {species} =================")
        for pair_tuple, pair_filepath in pair_paths_dir.items():
            if species not in pair_tuple:
                continue

            species_partner = get_pair_partner(pair_tuple=pair_tuple, species=species)
            df=pd.read_csv(pair_filepath, sep="\t")

            chr_colname = "chromosome"
            if pair_tuple[1]==species:
                chr_colname = "chromosome.1"
            
            print(f"\t* {species_partner} : \t{pair_filepath}")
            for i,ortholog_line in df.iterrows():
                species_geneID = ortholog_line[species]
                partner_geneID = ortholog_line[species_partner]
                chromosome = ortholog_line[chr_colname]
                if species_geneID not in species_orthologs:
                    species_orthologs[species_geneID] = SpeciesOrtholog(geneID=species_geneID, species=species, orthologs={partner_geneID:species_partner}, chromosome=chromosome)
                else:
                    species_orthologs[species_geneID].add_ortholog(key_value_tuple = (partner_geneID,species_partner))
        if False:
            if species=="A_obtectus":
                print(species_orthologs["rna-AOBTE_LOCUS3-2"])
            elif species=="B_siliquastri":
                print(species_orthologs["BRAKERILHT00000000982"])
            elif species=="C_chinensis":
                print(species_orthologs["CALCHIM00000001954"])
            elif species=="C_maculatus":
                print(species_orthologs["g1087.t1"])

        orthologs_dict[species]=species_orthologs
        
    print("\n")
    for species, orthologs_dict_sp in orthologs_dict.items():
        print(f"{species} has {len(orthologs_dict_sp)} orthologs (Type: {type(orthologs_dict_sp)})")
    
    # loop through all orthologs in species_list[0] and check their orthologs with the full ortholog dictionary. 
    # Make list of all of these gene IDs and check if list(set()) gives the same length as species list (only four unique ones)
    # immediately skip past geneIDs that aren't represented in one other species
    # also skip past ones where the chromosome assignment is not uniquely A or X but mixed

    focal_species=species_list[0]
    num_species = len(species_list)
    skipped_focal_geneID_list = []
    wrong_chr = 0
    all_OGs = len(orthologs_dict[focal_species])

    one_to_ones = {}
    OG_id = 1
    for focal_geneID, ortholog_instance in orthologs_dict[focal_species].items():
        focal_chr = ortholog_instance.chromosome
        orthologs_list =[(focal_geneID,focal_species)] # make list of tuples (geneID,species) to avoid duplication overwrite stuff that will happen with dicts 

        # skip focal orthologs that don't have the correct number of partners
        if len(ortholog_instance.orthologs) != len(species_list)-1:
            skipped_focal_geneID_list.append(focal_geneID)
            continue
        
        # check the orthologs for all partners
        for partnerID,partner_species in ortholog_instance.orthologs.items():
            # this should not be possible, at least one reciprocal ortholog of focal_geneID should be in the partner species dict under partnerID
            if partnerID not in orthologs_dict[partner_species]:
                print(orthologs_dict[partner_species], len(orthologs_dict[partner_species]))
                print(f"focal species : {focal_species}")
                print(orthologs_dict[focal_species][focal_geneID])
                print(f"partner species : {partner_species}")
                print(orthologs_dict[partner_species][partnerID])
                raise RuntimeError(f"""
{focal_geneID} in {focal_species} does not have a reciprocal ortholog in {partner_species} for partner ortholog {partnerID}!
Something is wrong when making the dictionaries
                """)
            
            # skip partner orthologs that don't have the correct number of partners
            partner_orthologs_dict = orthologs_dict[partner_species][partnerID].orthologs
            if len(partner_orthologs_dict) != num_species-1:
                skipped_focal_geneID_list.append(focal_geneID)
                break  
            # skip partner orthologs that are not on the same chromosome as the focal ortholog
            partner_chr = orthologs_dict[partner_species][partnerID].chromosome
            if partner_chr != focal_chr:
                skipped_focal_geneID_list.append(focal_geneID)
                wrong_chr+=1
                break
            
            for geneID,species in partner_orthologs_dict.items():
                partner_tupel = (geneID,species)
                orthologs_list.append(partner_tupel)
        
        orthologs_list = list(set(orthologs_list))
        orthologs_species_list = list(set([species for gID,species in orthologs_list]))
        # skip when unique set of orthologs and species is not correct
        if len(orthologs_list) != num_species or len(orthologs_species_list) != num_species:
            # print(orthologs_list)
            skipped_focal_geneID_list.append(focal_geneID)
            continue
        
        # write to dict if all requirements are fulfilled
        onetoone_orthologs_dict = {species : geneID for geneID,species in orthologs_list}
        one_to_ones[OG_id] = OrthoGroup(groupID=OG_id,species_dict=onetoone_orthologs_dict, chromosome=focal_chr)
        OG_id+=1

    skipped_focal_geneID_list = list(set(skipped_focal_geneID_list))
    intersection = list(set(skipped_focal_geneID_list) & set(one_to_ones))
    skipped_list_len = len(skipped_focal_geneID_list)
    sum_OGs = len(one_to_ones) + skipped_list_len
    print(f"""
{skipped_list_len} out of {all_OGs} A. obtectus orthologs do not fulfill the requirements for 1-to-1 in all species.
({wrong_chr} were skipped because they were not exclusive to X or A in all species)
{len(one_to_ones)} 1-to-1 orthologs with all {len(species_list)} species detected
{len(one_to_ones)} + {skipped_list_len } = {sum_OGs} (should be {all_OGs}, difference is {sum_OGs - all_OGs}), length of intersection: {len(intersection)})
""")

    with open(outfile, "w") as output_file:
        header = ",".join(species_list)
        output_file.write(f"orthogroup_ID,chromosome,{header}\n") 
        for OG_id, orthogroup in one_to_ones.items():
            geneIDs = ",".join([orthogroup.species_dict[species] for species in species_list])
            outline = f"{orthogroup.groupID},{orthogroup.chromosome},{geneIDs}\n"
            output_file.write(outline)

    print(f"\n-- written to output file: {outfile}")

            

                
            
            


if __name__ == "__main__":

    username="miltr339"
    brhs_tables = get_paths(username=username)

    get_all_orthologs(pair_paths_dir=brhs_tables, 
        species_list = ["A_obtectus","B_siliquastri","C_chinensis","C_maculatus"], 
        outfile=f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_orthologs.csv")
