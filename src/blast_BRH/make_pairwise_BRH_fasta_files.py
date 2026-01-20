"""
make fasta files with the two BRH sequences from the blast results
to calculate the dNdS ratio
"""

from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
import os 

def nucleotides_paths(username = "miltr339"):
    nuc_dir = f"/Users/{username}/work/native_nucleotides/"
    filepaths_dict = {
        "A_obtectus" : f"{nuc_dir}A_obtectus_transcripts.fna",
        "A_verrucosus" : f"{nuc_dir}A_verrucosus_transcripts.fna",
        "B_siliquastri" : f"{nuc_dir}B_siliquastri_transcripts.fna",
        "C_chinensis" : f"{nuc_dir}C_chinensis_transcripts.fna",
        "C_maculatus" : f"{nuc_dir}C_maculatus_transcripts.fna",
        "C_septempunctata" : f"{nuc_dir}C_septempunctata_transcripts.fna",
        "C_magnifica" : f"{nuc_dir}C_magnifica_transcripts.fna",
        "D_carinulata" : f"{nuc_dir}D_carinulata_transcripts.fna",
        "D_melanogaster" : f"{nuc_dir}D_melanogaster_transcripts.fna",
        "D_ponderosae" : f"{nuc_dir}D_ponderosae_transcripts.fna",
        "D_sublineata" : f"{nuc_dir}D_sublineata_transcripts.fna",
        "I_luminosus" : f"{nuc_dir}I_luminosus_transcripts.fna",
        "P_pyralis" : f"{nuc_dir}P_pyralis_transcripts.fna",
        "R_ferrugineus" : f"{nuc_dir}R_ferrugineus_transcripts.fna",
        "T_castaneum" : f"{nuc_dir}T_castaneum_transcripts.fna",
        "T_molitor" : f"{nuc_dir}T_molitor_transcripts.fna",
        "Z_morio" : f"{nuc_dir}Z_morio_transcripts.fna",
    }
    return filepaths_dict


def brh_results(username = "miltr339", X_syntenic = False):
    brh_dir = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/"
    brh_tables_X_syntenic = {
        "A_obtectus" : {
            "A_obtectus" : f"{brh_dir}A_obtectus_A_obtectus_BRH.tsv",
            "B_siliquastri" : f"{brh_dir}A_obtectus_B_siliquastri_BRH.tsv",
            "C_chinensis" : f"{brh_dir}A_obtectus_C_chinensis_BRH.tsv",
            "C_maculatus" : f"{brh_dir}A_obtectus_C_maculatus_BRH.tsv",
            # "D_carinulata" : f"{brh_dir}A_obtectus_D_carinulata_BRH_X_syntenic.tsv",
            # "D_sublineata" : f"{brh_dir}A_obtectus_D_sublineata_BRH.tsv",
        },
        "B_siliquastri" : {
            "B_siliquastri" : f"{brh_dir}B_siliquastri_B_siliquastri_BRH.tsv",
            "C_chinensis" : f"{brh_dir}B_siliquastri_C_chinensis_BRH.tsv",
            "C_maculatus" : f"{brh_dir}B_siliquastri_C_maculatus_BRH.tsv",
            # "D_carinulata" : f"{brh_dir}B_siliquastri_D_carinulata_BRH_X_syntenic.tsv",
            # "D_sublineata" : f"{brh_dir}B_siliquastri_D_sublineata_BRH.tsv",
        },
        "C_chinensis" : {
            "C_chinensis" : f"{brh_dir}C_chinensis_C_chinensis_BRH.tsv",
            "C_maculatus" : f"{brh_dir}C_chinensis_C_maculatus_BRH.tsv",
            # "D_carinulata" : f"{brh_dir}C_chinensis_D_carinulata_BRH_X_syntenic.tsv",
            # "D_sublineata" : f"{brh_dir}C_chinensis_D_sublineata_BRH.tsv",
        },
        "C_maculatus" : {
            "C_maculatus" : f"{brh_dir}C_maculatus_C_maculatus_BRH.tsv",
            # "D_carinulata" : f"{brh_dir}C_maculatus_D_carinulata_BRH_X_syntenic.tsv",
            # "D_sublineata" : f"{brh_dir}C_maculatus_D_sublineata_BRH.tsv",
        },
        "D_carinulata" : {
            "D_carinulata" : f"{brh_dir}D_carinulata_D_carinulata_BRH.tsv",
            "D_sublineata" : f"{brh_dir}D_carinulata_D_sublineata_BRH_X_syntenic.tsv",
        },
        "D_sublineata" : {
            "D_sublineata" : f"{brh_dir}D_sublineata_D_sublineata_BRH.tsv",
        }
    }
    brh_tables = {
        # "A_obtectus" : {
        #     "A_obtectus" : f"{brh_dir}A_obtectus_A_obtectus_BRH.tsv",
        #     "B_siliquastri" : f"{brh_dir}A_obtectus_B_siliquastri_BRH.tsv",
        #     "C_chinensis" : f"{brh_dir}A_obtectus_C_chinensis_BRH.tsv",
        #     "C_maculatus" : f"{brh_dir}A_obtectus_C_maculatus_BRH.tsv",
        #     # "D_carinulata" : f"{brh_dir}A_obtectus_D_carinulata_BRH.tsv",
        #     # "D_sublineata" : f"{brh_dir}A_obtectus_D_sublineata_BRH.tsv",
        # },
        # "B_siliquastri" : {
        #     "B_siliquastri" : f"{brh_dir}B_siliquastri_B_siliquastri_BRH.tsv",
        #     "C_chinensis" : f"{brh_dir}B_siliquastri_C_chinensis_BRH.tsv",
        #     "C_maculatus" : f"{brh_dir}B_siliquastri_C_maculatus_BRH.tsv",
        #     # "D_carinulata" : f"{brh_dir}B_siliquastri_D_carinulata_BRH.tsv",
        #     # "D_sublineata" : f"{brh_dir}B_siliquastri_D_sublineata_BRH.tsv",
        # },
        # "C_chinensis" : {
        #     "C_chinensis" : f"{brh_dir}C_chinensis_C_chinensis_BRH.tsv",
        #     "C_maculatus" : f"{brh_dir}C_chinensis_C_maculatus_BRH.tsv",
        #     # "D_carinulata" : f"{brh_dir}C_chinensis_D_carinulata_BRH.tsv",
        #     # "D_sublineata" : f"{brh_dir}C_chinensis_D_sublineata_BRH.tsv",
        # },
        # "C_maculatus" : {
        #     "C_maculatus" : f"{brh_dir}C_maculatus_C_maculatus_BRH.tsv",
        #     # "D_carinulata" : f"{brh_dir}C_maculatus_D_carinulata_BRH.tsv",
        #     # "D_sublineata" : f"{brh_dir}C_maculatus_D_sublineata_BRH.tsv",
        # },
        # "D_carinulata" : {
        #     "D_carinulata" : f"{brh_dir}D_carinulata_D_carinulata_BRH.tsv",
        #     "D_sublineata" : f"{brh_dir}D_carinulata_D_sublineata_BRH.tsv",
        # },
        # "D_sublineata" : {
        #     "D_sublineata" : f"{brh_dir}D_sublineata_D_sublineata_BRH.tsv",
        # },
        "T_castaneum" : {
            "T_freemani" : f"{brh_dir}T_castaneum_T_freemani_BRH.tsv",
        },
        "C_magnifica" : {
            "C_septempunctata" : f"{brh_dir}C_magnifica_C_septempunctata_BRH.tsv",
        }

    }
    if X_syntenic:
        return brh_tables_X_syntenic
    else:
        return brh_tables


def make_ortholog_fasta_files(brh_tables, nucleotides_dict, chr_type = "X", outdir = ""):
    """
    make fasta files for dNdS from pairwise blast results
    with chromosome type you can restrict the orthologs to X-exclusive or A-exclusive by setting "X" or "A"
    """

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    print(f"\nmaking nucleotide fasta files for {chr_type}-linked orthologs")
    nucleotides_seqs_dict = {}
    print(f"reading nucleotide files...")
    for species in brh_tables.keys():
        print(f"> {species} ({nucleotides_dict[species]})")
        # make dict from the nucleotide files
        nucleotides_seqs_dict[species] = {}
        for cds_record in SeqIO.parse(nucleotides_dict[species], "fasta"):
            nucleotides_seqs_dict[species][cds_record.id] = cds_record
    
    for species in brh_tables.keys():
        print(f">>>>>> {species}")
        not_found_partners_dict = {sp:0 for sp in brh_tables[species].keys()}
        for species_partner, brh_path in brh_tables[species].items(): # tqdm(brh_tables[species].items()):
            count_pairs = 0
            ## read and filter to only include hits on chr_type
            brh_table = pd.read_csv(brh_path, sep="\t")
            print(brh_table)
            brh_filtered = brh_table[brh_table["chromosome"] == chr_type]
            print(f"after filtering {species} for {chr_type}...")
            print(brh_filtered)
            brh_filtered = brh_filtered[brh_filtered["chromosome.1"] == chr_type]
            print(f"\t{species_partner} : {len(brh_filtered)} orthologs on {chr_type} ({brh_path})")
            
            pair_dirname = f"{outdir}{species}_{species_partner}/"
            if not os.path.isdir(pair_dirname):
                os.mkdir(pair_dirname)

            for brh_pair in brh_filtered.itertuples():
                # print(f"ID1 = {brh_pair[1]} , ID2 = {brh_pair[3]}")
                try:
                    record1 = nucleotides_seqs_dict[species][brh_pair[1]]
                    record1.id=f">{brh_pair[1]}_{species}_{chr_type}"
                except:
                    not_found_partners_dict[species]+=1
                    # raise RuntimeError(f"{brh_pair[1]} not found in {nucleotides_dict[species]}")
                try:
                    record2 = nucleotides_seqs_dict[species_partner][brh_pair[3]]
                    record2.id=f">{brh_pair[3]}_{species_partner}_{chr_type}"
                except:
                    not_found_partners_dict[species_partner]+=1
                    # raise RuntimeError(f"{brh_pair[3]} not found in {nucleotides_dict[species_partner]}")
                
                brh_seq_records = [record1,record2]
                fasta_name = f"{pair_dirname}{species}_{species_partner}_{chr_type}-linked_ortholog_{count_pairs}.fasta"
                SeqIO.write(sequences=brh_seq_records, handle=fasta_name, format="fasta")
                # print(f"{record1}\n{record2}")
                # print()
                # break
                count_pairs += 1

            # print(f"\t* {species_partner}: {count_pairs} orthologs")
            print(f"\t\t{count_pairs} counted, {not_found_partners_dict[species_partner]} transcripts not found")
    print(f"\ndone! all outfiles written to \n{outdir}")






if __name__ == "__main__":
    username = "miltr339"
    Dcar_X_syntenic = False
    brh_tables = brh_results(username, X_syntenic=Dcar_X_syntenic)
    nucleotides_dict = nucleotides_paths(username)
    
    if Dcar_X_syntenic:
        outdir_X = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_X_Dcar_X_syntenic/"
        outdir_A = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_A_Dcar_X_syntenic/"
    else:
        outdir_X = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_X/"
        outdir_A = f"/Users/{username}/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_A/"
    
    make_ortholog_fasta_files(brh_tables, nucleotides_dict, chr_type="X", outdir= outdir_X)
    # make_ortholog_fasta_files(brh_tables, nucleotides_dict, chr_type="A", outdir= outdir_A)